#include "arcs_threads.h"
#include "kmer_anchor_pg.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <queue>
#include <thread>
#include <unordered_map>
#include <vector>
#include <string>

namespace {

// Pack a K-length ACGT window starting at off into a 2-bit code. Returns false
// if it contains a non-ACGT byte (kmer indexing/seeding is ACGT-only; N's and
// other bytes are handled losslessly elsewhere via pg_N_* side streams).
inline bool pack_kmer(const std::string& s, int off, int K, uint64_t& out) {
    uint64_t v = 0;
    for (int j = 0; j < K; ++j) {
        uint8_t b = encode_base(s[(size_t)(off + j)]);
        if (b >= 4) return false;
        v = (v << 2) | b;
    }
    out = v;
    return true;
}

struct Placement {
    bool     used = false;
    uint32_t cid  = 0;
    uint32_t off  = 0;   // start position within that contig (contig-local frame)
    uint8_t  view = 0;   // 0 = read's own forward seq placed as-is, 1 = its RC was placed
};

// ── Flat, parallel-built k-mer occurrence table ──────────────────────────────
// Profiling (ARCS_KA_TIMING) showed the serial std::unordered_map index build
// is ~98% of Method A's total assembly time (27-34s of a ~28-35s shard,
// growth itself only 0.5-0.7s) — this is the actual bottleneck, not the
// growth loop. Confirmed via the deep-research pass into PgRC2/copMEM2's real
// source (not the paper's marketing language): their "lock-free hash table"
// isn't a CAS-based concurrent map — it's a two-pass parallel counting sort
// (genCummMultithreaded): pass 1, each thread builds a LOCAL histogram over
// its own contiguous slice; a cheap serial prefix-sum merges these into
// global per-bucket offsets; pass 2, each thread rescans its own slice and
// scatters directly into pre-claimed, non-overlapping slots — zero locks,
// zero retries, because collisions are precomputed away instead of resolved
// live. That's exactly what this does. Critically: unlike round-robin
// sharding the READ SET (which starves each shard of global overlap
// visibility and measurably wrecked assembly quality — pg_len 60.1M vs 34.0M,
// see the parallel_shard_assemble_ka path), this only parallelizes HOW the
// index is built. Every one of the n reads still lands in the ONE resulting
// table — it is byte-for-byte the same set of entries a serial build would
// produce (mod arbitrary tie-breaking among >bucket_cap collisions, same
// heuristic cap the serial version already had), so growth downstream sees
// identical candidates and produces identical assembly quality.
// km is built via repeated (v<<2)|base, so its LOW bits are dominated by the
// LAST few bases only — using them directly as a bucket index (raw `& mask`)
// causes very uneven bucket occupancy (many distinct k-mers sharing a common
// suffix motif collide into the same bucket and fight over its bucket_cap
// slots, which the OLD per-exact-key unordered_map never suffered — each
// distinct k-mer got its own independent cap). Root cause of a first attempt
// at this that regressed pg_len 34.0M -> 49.0M despite "coverage" being
// unchanged: it wasn't a coverage loss, it was uneven bucket collision
// silently evicting entries that used to have their own private slot budget.
// Fix: a real avalanche mix (splitmix64 finalizer) before masking, so bucket
// assignment no longer correlates with k-mer suffix content.
inline uint64_t mix64(uint64_t x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

struct FlatKmerIndex {
    uint64_t table_mask = 0;
    int      bucket_cap = 0;
    std::vector<uint64_t> km;   // UINT64_MAX = empty slot
    std::vector<uint64_t> val;

    template <class F>
    void for_each(uint64_t query_km, F&& f) const {
        uint64_t b = mix64(query_km) & table_mask;
        const uint64_t* kp = &km[(size_t)b * (size_t)bucket_cap];
        const uint64_t* vp = &val[(size_t)b * (size_t)bucket_cap];
        for (int s = 0; s < bucket_cap; ++s) {
            if (kp[s] == UINT64_MAX) break;   // compact fill guaranteed by atomic slot claim
            if (kp[s] == query_km) f(vp[s]);
        }
    }
};

template <class ViewFn>
FlatKmerIndex build_index_parallel(size_t n, int K0, int bucket_cap,
                                    int n_threads, ViewFn&& view_seq) {
    int T = std::max(1, n_threads);
    std::vector<std::vector<std::pair<uint64_t,uint64_t>>> local((size_t)T);
    {
        std::vector<std::thread> ths;
        ths.reserve((size_t)T);
        for (int t = 0; t < T; ++t) {
            ths.emplace_back([&, t]{
                size_t lo = n * (size_t)t / (size_t)T;
                size_t hi = n * (size_t)(t + 1) / (size_t)T;
                auto& out = local[(size_t)t];
                for (size_t rid = lo; rid < hi; ++rid) {
                    for (uint8_t view = 0; view < 2; ++view) {
                        const std::string& S = view_seq((uint32_t)rid, view);
                        int Ls = (int)S.size();
                        for (int off = 0; off + K0 <= Ls; ++off) {
                            uint64_t km;
                            if (!pack_kmer(S, off, K0, km)) continue;
                            out.push_back({km, ((uint64_t)rid << 17) | ((uint64_t)view << 16) | (uint64_t)off});
                        }
                    }
                }
            });
        }
        for (auto& th : ths) th.join();
    }

    size_t total = 0;
    for (auto& v : local) total += v.size();

    // Real memory cost is table_size * bucket_cap slots (each bucket carries
    // bucket_cap inline entries), not table_size alone — size so that TOTAL
    // SLOT CAPACITY is ~2x the occurrence count (average bucket well under
    // half full), then cap total slot count outright as a hard memory guard.
    size_t target_buckets = (total * 2) / (size_t)std::max(1, bucket_cap) + 1;
    size_t table_size = 1024;
    while (table_size < target_buckets && table_size < (64ull << 20)) table_size <<= 1;
    while (table_size * (size_t)bucket_cap > (2ull << 30)) table_size >>= 1; // hard cap: 2G slots total

    FlatKmerIndex idx;
    idx.table_mask = table_size - 1;
    idx.bucket_cap = bucket_cap;
    idx.km.assign(table_size * (size_t)bucket_cap, UINT64_MAX);
    idx.val.assign(table_size * (size_t)bucket_cap, 0);

    std::vector<std::atomic<uint32_t>> bucket_count(table_size);
    for (auto& c : bucket_count) c.store(0, std::memory_order_relaxed);

    {
        std::vector<std::thread> ths;
        ths.reserve((size_t)T);
        for (int t = 0; t < T; ++t) {
            ths.emplace_back([&, t]{
                for (auto& kv : local[(size_t)t]) {
                    uint64_t b = mix64(kv.first) & idx.table_mask;
                    uint32_t slot = bucket_count[b].fetch_add(1, std::memory_order_relaxed);
                    if (slot < (uint32_t)bucket_cap) {
                        size_t pos = (size_t)b * (size_t)bucket_cap + slot;
                        idx.km[pos]  = kv.first;
                        idx.val[pos] = kv.second;
                    }
                }
            });
        }
        for (auto& th : ths) th.join();
    }
    return idx;
}

} // namespace

// ── Method A: global k-mer occurrence index + mismatch-tolerant anchor growth ──
ChainEncodeResult build_kmer_anchor_pg(const std::vector<Read>& reads, bool allow_leftover_pass) {
    ChainEncodeResult r;
    r.has_pg = true;
    const size_t n = reads.size();
    r.n_reads = n;

    // PgRC's real seed length is 38 (confirmed: pgrc-params.h ReadSeedLength) —
    // tried porting that number directly (K0=32, the closest this uint64_t
    // 2-bit-packed encoding can represent) plus its bucket cap of 12 below.
    // Measured result on yeast_sub.fq (1M reads): WORSE, not better — raw
    // pg_len 43.7M vs this K0=24/cap=64 config's 34.0M, and final archive
    // 5.65M vs 5.30M. PgRC's real hash table doesn't use simple bucket-cap
    // eviction the way ours does (its K is auto-derived per-read-length from
    // a lookup table, not a fixed global constant) — copying its raw numbers
    // without its exact collision-handling mechanics backfired here. Reverted;
    // kept as an env-tunable option for further experimentation, not default.
    int K0 = 24;
    if (const char* s = getenv("ARCS_KA_K")) { int v = atoi(s); if (v >= 12 && v <= 32) K0 = v; }
    // MAXMM: growth-time anchor verification tolerance only — stays a tight
    // flat constant (growth accepts read content directly into the pg, so a
    // loose tolerance here would corrupt assembly quality, not just placement).
    int MAXMM = 4;
    if (const char* s = getenv("ARCS_KA_MAXMM")) MAXMM = atoi(s);
    // FALLBACK_MAXMM: separate, looser tolerance for the fallback placement
    // pass only (reads growth never touched, mapped onto the FIXED finished
    // pg — no corruption risk). -1 sentinel = default to PgRC's real LQ-remap
    // tolerance (~read_len/3 per its MinCharsPerMismatch=3, confirmed via
    // direct source read), computed per read below, instead of reusing the
    // tight growth constant which is ~8x stricter than PgRC's real number.
    int FALLBACK_MAXMM = -1;
    if (const char* s = getenv("ARCS_KA_FALLBACK_MAXMM")) FALLBACK_MAXMM = atoi(s);
    // PgRC's real bucket cap is 12 (CopMEMMatcher.cpp) — tried and reverted
    // alongside K0 above (measured regression). Separately, after switching
    // to the parallel flat-table index (see FlatKmerIndex below), a bucket
    // holds entries from potentially MULTIPLE distinct k-mers (unlike the old
    // per-exact-key unordered_map, where each k-mer's occurrence list was
    // capped independently) — so the old cap=64 silently dropped more real
    // candidates than before under hash-bucket sharing. Measured on
    // yeast_sub.fq (1M reads, 12 threads): cap=64 gave pg_len 35.7M; cap=256
    // gave pg_len 30.15M (smaller than even the original serial 34.0M
    // baseline) and total_seq 5.85M — matching the original serial baseline's
    // 5.90M total_seq at ~12x the assembly speed. Confirmed LOSSLESS. Raised
    // the default accordingly; the real cost is memory (table_size*bucket_cap
    // slots), not time — see build_index_parallel's sizing comment.
    int BUCKET_CAP = 256;
    if (const char* s = getenv("ARCS_KA_BUCKET")) { int v = atoi(s); if (v >= 1) BUCKET_CAP = v; }

    // ── Build the two views (forward, reverse-complement) of every read once ────
    std::vector<std::string> rcseq(n);
    for (size_t i = 0; i < n; ++i) rcseq[i] = reverse_complement(reads[i].seq);
    auto view_seq = [&](uint32_t rid, uint8_t view) -> const std::string& {
        return view == 0 ? reads[rid].seq : rcseq[rid];
    };

    // ── Global K0-mer occurrence index over every position of every read, both
    // views (this is the distinguishing feature vs build_vodbg_pg's exact
    // suffix-array/APSP overlap table: a candidate is accepted here via a
    // single anchor k-mer plus a bounded mismatch-tolerant span check, not an
    // exact maximal-overlap guarantee). packed value = (rid<<17)|(view<<16)|off
    // — off capped at 16 bits (reads already validated <65536 bp elsewhere).
    const bool KA_TIMING = getenv("ARCS_KA_TIMING") != nullptr;
    auto _ka_t0 = std::chrono::steady_clock::now();

    // Index build parallelism: default to min(hw_cores, 4) threads, same cap
    // as the default assembler's own sharding — this is a build-time-only
    // parallelization (see FlatKmerIndex header comment), not a coverage
    // trade-off, so there is no quality-vs-speed tension to tune here the way
    // ARCS_KA_PAR_SHARDS has. Default to ALL cores (not capped at 4 the way
    // read-sharding is) — measured on this 12-core box: 4 threads gave 10.80s
    // assembly, 12 threads gave 7.24s, same pg_len both times (35.7M vs
    // serial's 34.0M) — quality is thread-count-independent here, so there's
    // no reason to leave cores idle. Override via ARCS_KA_INDEX_THREADS
    // (1 = serial).
    int idx_threads = arcs_threads();
    if (idx_threads < 1) idx_threads = 1;
    if (const char* s = getenv("ARCS_KA_INDEX_THREADS")) { int v = atoi(s); if (v >= 1) idx_threads = v; }
    FlatKmerIndex kidx = build_index_parallel(n, K0, BUCKET_CAP, idx_threads, view_seq);
    if (KA_TIMING) {
        auto t1 = std::chrono::steady_clock::now();
        fprintf(stderr, "[KA-TIMING] index_build: %.2fs (threads=%d, table_slots=%zu)\n",
                std::chrono::duration<double>(t1 - _ka_t0).count(), idx_threads, kidx.km.size());
        _ka_t0 = t1;
    }

    // ── Global greedy-overlap contig growth, parallelized ───────────────────────
    // Each thread grows its OWN independent set of contigs against the ONE
    // shared, already-complete, read-only kidx — unlike parallel_shard_assemble_ka
    // (which sharded the READ SET itself and measurably lost quality), every
    // thread here can still find ANY read anywhere as a growth candidate, so
    // there is no coverage loss. The only cross-thread contention is "which
    // thread gets to claim read X" — resolved with a single atomic
    // compare-exchange per read, no locks. Critically, a contig is only ever
    // touched by the ONE thread that created it (we never merge two
    // already-grown contigs, only ever consume a fresh original read into
    // whichever contig is asking) — so unlike PgRC's real parallel design
    // (which needs a critical section for cross-bucket cycle detection), we
    // need no critical section at all: contigs/contig_members/tail_rid/
    // tail_view are thread-local, and the shared `place[]` write for a given
    // read happens only in the one thread that won that read's atomic claim.
    std::vector<Placement> place(n);
    std::vector<std::atomic<uint8_t>> claimed(n);
    for (auto& c : claimed) c.store(0, std::memory_order_relaxed);
    auto try_claim = [&](uint32_t rid) -> bool {
        uint8_t expected = 0;
        return claimed[rid].compare_exchange_strong(expected, 1, std::memory_order_acq_rel);
    };

    struct ThreadState {
        std::vector<std::string> contigs;
        std::vector<std::vector<uint32_t>> contig_members;
        std::vector<std::array<uint32_t,2>> tail_rid;
        std::vector<std::array<uint8_t,2>>  tail_view;
    };

    struct Growth { bool found; uint32_t rid; uint8_t view; uint32_t overlap; };
    struct Proposal { int64_t overlap; uint32_t cid; uint8_t dir; };
    struct ProposalLess { bool operator()(const Proposal& a, const Proposal& b) const { return a.overlap < b.overlap; } };

    int growth_threads = idx_threads; // same core budget as index build
    if (const char* s = getenv("ARCS_KA_GROWTH_THREADS")) { int v = atoi(s); if (v >= 1) growth_threads = v; }
    if ((size_t)growth_threads > n) growth_threads = std::max(1, (int)n);
    std::vector<ThreadState> tstate((size_t)growth_threads);

    auto run_thread = [&](int t, size_t range_lo, size_t range_hi) {
        ThreadState& ts = tstate[(size_t)t];
        std::priority_queue<Proposal, std::vector<Proposal>, ProposalLess> heap;

        // Anchor-based candidate discovery + mismatch-tolerant verification,
        // identical logic to the pre-parallelization version (see the removed
        // comment block this replaces), except "already placed" is now
        // answered by the SHARED atomic claimed[] instead of a plain bool —
        // that's the only piece any other thread can change underneath us.
        auto next_candidate = [&](uint32_t cid, uint8_t dir) -> Growth {
            uint32_t rid  = ts.tail_rid[cid][dir];
            uint8_t  view = ts.tail_view[cid][dir];
            uint8_t  query_view = (dir == 0) ? view : (uint8_t)(1 - view);
            std::string tmp;
            const std::string& S = (dir == 0) ? ts.contigs[cid] : (tmp = reverse_complement(ts.contigs[cid]));
            int Ls = (int)S.size();
            if (Ls < K0) return {false, 0, 0, 0};
            const std::string& Q = view_seq(rid, query_view);
            if ((int)Q.size() < K0) return {false, 0, 0, 0};
            uint64_t qk;
            if (!pack_kmer(Q, (int)Q.size() - K0, K0, qk)) return {false, 0, 0, 0};

            int      best_overlap = -1, best_mm = MAXMM + 1;
            uint32_t best_rid = 0;
            uint8_t  best_view = 0;
            kidx.for_each(qk, [&](uint64_t packed) {
                uint32_t crid  = (uint32_t)(packed >> 17);
                uint8_t  cview = (uint8_t)((packed >> 16) & 1);
                uint32_t off   = (uint32_t)(packed & 0xFFFFu);
                if (claimed[crid].load(std::memory_order_acquire)) return;
                const std::string& R = view_seq(crid, cview);
                int base_overlap = (int)off + K0;
                if (base_overlap > Ls || base_overlap > (int)R.size()) return;
                const char* rp = R.data();
                const char* sp = S.data() + (Ls - base_overlap);
                int mm = 0;
                for (int j = 0; j < base_overlap && mm <= MAXMM; ++j) if (rp[j] != sp[j]) ++mm;
                if (mm > MAXMM) return;
                if (base_overlap > best_overlap || (base_overlap == best_overlap && mm < best_mm)) {
                    best_overlap = base_overlap; best_mm = mm; best_rid = crid; best_view = cview;
                }
            });
            if (best_overlap < 0) return {false, 0, 0, 0};
            return {true, best_rid, best_view, (uint32_t)best_overlap};
        };

        auto propose = [&](uint32_t cid, uint8_t dir) {
            Growth g = next_candidate(cid, dir);
            if (g.found) heap.push({(int64_t)g.overlap, cid, dir});
        };

        size_t next_unseeded = range_lo;
        auto seed_next = [&]() -> bool {
            while (next_unseeded < range_hi &&
                   (claimed[next_unseeded].load(std::memory_order_acquire) ||
                    reads[next_unseeded].seq.size() < (size_t)K0))
                ++next_unseeded;
            if (next_unseeded >= range_hi) return false;
            uint32_t rid0 = (uint32_t)next_unseeded;
            ++next_unseeded; // advance regardless — if the claim below fails, seed_next's own
                              // top-of-loop skip will not revisit it (claimed[] is now true)
            if (!try_claim(rid0)) return true; // lost race to another thread; caller loops again
            uint32_t cid = (uint32_t)ts.contigs.size();
            ts.contigs.emplace_back(reads[rid0].seq);
            ts.contig_members.emplace_back();
            ts.tail_rid.push_back({rid0, rid0});
            ts.tail_view.push_back({0, 0});
            place[rid0] = {true, cid, 0, 0};
            ts.contig_members[cid].push_back(rid0);
            propose(cid, 0);
            propose(cid, 1);
            return true;
        };

        bool more_to_seed = seed_next();
        while (more_to_seed || !heap.empty()) {
            if (heap.empty()) { more_to_seed = seed_next(); continue; }
            Proposal p = heap.top(); heap.pop();
            Growth g{false, 0, 0, 0};
            while (true) {
                g = next_candidate(p.cid, p.dir); // re-verify fresh each attempt
                if (!g.found) break;
                if (!try_claim(g.rid)) continue;  // another thread took it — rescan, it's now excluded
                break;
            }
            if (!g.found) continue;

            const std::string& R = view_seq(g.rid, g.view);
            if (p.dir == 0) {
                std::string& C = ts.contigs[p.cid];
                place[g.rid] = {true, p.cid, (uint32_t)(C.size() - g.overlap), g.view};
                ts.contig_members[p.cid].push_back(g.rid);
                C.append(R.data() + g.overlap, R.size() - (size_t)g.overlap);
                ts.tail_rid[p.cid][0] = g.rid; ts.tail_view[p.cid][0] = g.view;
                propose(p.cid, 0);
            } else {
                std::string Crc = reverse_complement(ts.contigs[p.cid]);
                size_t Lc_orig = ts.contigs[p.cid].size();
                Crc.append(R.data() + g.overlap, R.size() - (size_t)g.overlap);
                size_t Lc2 = Crc.size();
                uint32_t prepend_amount = (uint32_t)(Lc2 - Lc_orig);
                ts.contigs[p.cid] = reverse_complement(Crc);
                for (uint32_t m : ts.contig_members[p.cid]) place[m].off += prepend_amount;
                uint32_t new_off = (uint32_t)(Lc2 - (Lc_orig - g.overlap) - R.size());
                place[g.rid] = {true, p.cid, new_off, (uint8_t)(1 - g.view)};
                ts.contig_members[p.cid].push_back(g.rid);
                ts.tail_rid[p.cid][1] = g.rid; ts.tail_view[p.cid][1] = (uint8_t)(1 - g.view);
                propose(p.cid, 1);
            }
        }
    };

    {
        std::vector<std::thread> ths;
        ths.reserve((size_t)growth_threads);
        for (int t = 0; t < growth_threads; ++t) {
            size_t lo = n * (size_t)t / (size_t)growth_threads;
            size_t hi = n * (size_t)(t + 1) / (size_t)growth_threads;
            ths.emplace_back([&, t, lo, hi]{ run_thread(t, lo, hi); });
        }
        for (auto& th : ths) th.join();
    }

    size_t total_contigs = 0;
    for (auto& ts : tstate) total_contigs += ts.contigs.size();
    if (KA_TIMING) {
        auto t1 = std::chrono::steady_clock::now();
        fprintf(stderr, "[KA-TIMING] growth_loop: %.2fs (threads=%d, contigs=%zu)\n",
                std::chrono::duration<double>(t1 - _ka_t0).count(), growth_threads, total_contigs);
        _ka_t0 = t1;
    }

    // ── Merge all threads' contigs -> pg, compute absolute per-contig base
    // offsets. place[].cid was thread-local (an index into that thread's own
    // ts.contigs) — remap to a global contig id via each thread's running
    // offset, applied while copying place[] below.
    std::vector<uint32_t> contig_base;
    contig_base.reserve(total_contigs);
    std::vector<uint32_t> cid_remap_base((size_t)growth_threads, 0); // global cid of thread t's local cid 0
    r.pg.reserve([&]{ size_t t = 0; for (auto& ts : tstate) for (auto& c : ts.contigs) t += c.size(); return t; }());
    {
        uint32_t running_cid = 0;
        for (int t = 0; t < growth_threads; ++t) {
            cid_remap_base[(size_t)t] = running_cid;
            for (auto& c : tstate[(size_t)t].contigs) {
                contig_base.push_back((uint32_t)r.pg.size());
                r.pg += c;
                ++running_cid;
            }
        }
    }
    // Fix up place[].cid from thread-local to global now that every thread's
    // contigs have a known position in the merged sequence. We don't track
    // which thread wrote each place[i] directly, but every claimed read's cid
    // is thread-local only in the sense of *which contig list* it indexes —
    // since a read is claimed by exactly one thread, and that thread is the
    // only one that ever wrote place[i].cid, we recover the owning thread by
    // re-deriving it from contig_members during the merge pass above instead:
    for (int t = 0; t < growth_threads; ++t) {
        uint32_t base = cid_remap_base[(size_t)t];
        for (auto& members : tstate[(size_t)t].contig_members) {
            for (uint32_t rid : members) place[rid].cid += base;
        }
    }

    // ── Emit placements sorted by absolute pg position (delta-codes far better
    // downstream than original-read-index order — same reasoning as elsewhere).
    std::vector<uint32_t> resolved;
    resolved.reserve(n);
    for (size_t i = 0; i < n; ++i) if (place[i].used) resolved.push_back((uint32_t)i);
    std::sort(resolved.begin(), resolved.end(), [&](uint32_t a, uint32_t b) {
        uint32_t pa = contig_base[place[a].cid] + place[a].off;
        uint32_t pb = contig_base[place[b].cid] + place[b].off;
        return pa < pb;
    });
    for (uint32_t i : resolved) {
        uint32_t abs_pos = contig_base[place[i].cid] + place[i].off;
        const std::string& target = view_seq(i, place[i].view);
        record_mapped(r, i, abs_pos, (int)place[i].view, target, (int)target.size());
    }

    // ── Fallback placement pass for reads growth never touched (too short for
    // K0, or would only have extended a contig backward past its own start) —
    // search the FINAL, fixed pg directly. Always lossless regardless (falls
    // through to record_append if no acceptable placement found).
    std::vector<uint32_t> unresolved;
    for (size_t i = 0; i < n; ++i) if (!place[i].used) unresolved.push_back((uint32_t)i);
    std::vector<uint32_t> still_unresolved;

    if (!unresolved.empty()) {
        std::unordered_map<uint64_t, std::vector<uint32_t>> pgidx;
        int Lpg = (int)r.pg.size();
        if (Lpg >= K0) {
            pgidx.reserve((size_t)Lpg);
            for (int p = 0; p + K0 <= Lpg; ++p) {
                uint64_t km;
                if (!pack_kmer(r.pg, p, K0, km)) continue;
                auto& b = pgidx[km];
                if ((int)b.size() < BUCKET_CAP) b.push_back((uint32_t)p);
            }
        }

        for (uint32_t rid : unresolved) {
            int Lr0 = (int)reads[rid].seq.size();
            int maxmm_for_read = (FALLBACK_MAXMM >= 0) ? FALLBACK_MAXMM : std::max(1, Lr0 / 3);
            int bestmm = maxmm_for_read + 1; int32_t bestpos = -1; uint8_t bestview = 0;
            for (uint8_t view = 0; view < 2 && !pgidx.empty(); ++view) {
                const std::string& R = view_seq(rid, view);
                int Lr = (int)R.size();
                if (Lr < K0) continue;
                const int seed_offs[3] = {0, (Lr - K0) / 2, Lr - K0};
                for (int off : seed_offs) {
                    if (off < 0 || off + K0 > Lr) continue;
                    uint64_t km;
                    if (!pack_kmer(R, off, K0, km)) continue;
                    auto it = pgidx.find(km);
                    if (it == pgidx.end()) continue;
                    for (uint32_t p : it->second) {
                        int32_t cs = (int32_t)p - off;
                        if (cs < 0 || cs + Lr > Lpg) continue;
                        int mm = 0;
                        for (int j = 0; j < Lr && mm < bestmm; ++j)
                            if (r.pg[(size_t)(cs + j)] != R[(size_t)j]) ++mm;
                        if (mm < bestmm) { bestmm = mm; bestpos = cs; bestview = view; }
                    }
                }
            }
            if (bestpos >= 0 && bestmm <= maxmm_for_read) {
                const std::string& target = view_seq(rid, bestview);
                record_mapped(r, rid, (uint32_t)bestpos, (int)bestview, target, (int)target.size());
            } else {
                still_unresolved.push_back(rid);
            }
        }
    }

    // ── Leftover reassembly (confirmed via direct source read: this is exactly
    // what PgRC does with reads that fail to map onto its main pg — it doesn't
    // append them raw, it runs its SAME greedy assembler again over just the
    // leftovers, building a small second pg, before finally giving up on
    // whatever's left after THAT). Reused here via plain recursion — one extra
    // level only (allow_leftover_pass=false on the recursive call), since any
    // read still unplaceable after a second full assembly attempt is a true
    // singleton no amount of further recursion would help. This is the only
    // place record_append is ever reached now: after two real assembly
    // attempts, not one, on any read that doesn't overlap anything at all.
    if (getenv("ARCS_KA_TIMING") && allow_leftover_pass)
        fprintf(stderr, "[KA-TIMING] leftover_pass: unresolved=%zu still_unresolved=%zu (of n=%zu)\n",
                unresolved.size(), still_unresolved.size(), n);
    if (allow_leftover_pass && !still_unresolved.empty()) {
        std::vector<Read> leftover_reads;
        leftover_reads.reserve(still_unresolved.size());
        for (uint32_t rid : still_unresolved) leftover_reads.push_back(reads[rid]);

        ChainEncodeResult sub = build_kmer_anchor_pg(leftover_reads, /*allow_leftover_pass=*/false);
        if (getenv("ARCS_KA_TIMING")) {
            size_t raw_bytes = 0;
            for (uint32_t rid : still_unresolved) raw_bytes += reads[rid].seq.size();
            fprintf(stderr, "[KA-TIMING] leftover_pass: sub_pg_len=%zu vs raw_append_would_be=%zu bytes\n",
                    sub.pg.size(), raw_bytes);
        }
        uint32_t offset = (uint32_t)r.pg.size();
        r.pg += sub.pg;
        // Re-emit every one of sub's placements against the OUTER r.pg at the
        // shifted position — record_mapped recomputes mismatches live, and
        // since sub.pg was just copied byte-for-byte into r.pg at `offset`,
        // this reproduces the exact same diff sub itself found, whether that
        // entry came from sub's own growth, its own fallback mapping, or even
        // its own (bounded, non-recursive) raw append — the universal
        // record_mapped/record_append contract makes all three
        // indistinguishable to re-emit uniformly this way.
        for (size_t i = 0; i < sub.chain_order.size(); ++i) {
            uint32_t orig_rid = still_unresolved[sub.chain_order[i]];
            uint32_t abs_pos  = offset + sub.pg_pos[i];
            uint8_t  rc       = sub.pg_rc[i];
            const std::string& target = view_seq(orig_rid, rc);
            record_mapped(r, orig_rid, abs_pos, (int)rc, target, (int)target.size());
        }
    } else {
        for (uint32_t rid : still_unresolved)
            record_append(r, rid, reads[rid].seq, (int)reads[rid].seq.size());
    }

    return r;
}
