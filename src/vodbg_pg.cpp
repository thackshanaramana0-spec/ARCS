#include "vodbg_pg.h"
#include "sa_apsp.h"
#include "flat_kmer_index.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
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
// (bailing the caller to skip this position) if it contains a non-ACGT byte —
// same convention as the rest of the codebase (kmer seeding is ACGT-only; N's
// and other bytes are handled losslessly elsewhere via pg_N_* side streams).
// Only used by the fallback placement pass now (see below) — the main
// discovery/growth path uses exact suffix-array overlaps (sa_apsp.h), not
// k-mer anchors.
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

} // namespace

ChainEncodeResult build_vodbg_pg(const std::vector<Read>& reads, CallData* call_out) {
    ChainEncodeResult r;
    r.has_pg = true;
    const size_t n = reads.size();
    r.n_reads = n;

    int K0 = 24;
    if (const char* s = getenv("ARCS_VODBG_K")) { int v = atoi(s); if (v >= 12 && v <= 31) K0 = v; }
    // MAXMM: fallback-pass mismatch tolerance (the only place it's used — growth
    // is zero-tolerance by design, see HQ/LQ split above). -1 sentinel = default
    // to PgRC's real LQ-remap tolerance (~read_len/3 per its MinCharsPerMismatch=3
    // config, confirmed via direct source read) computed per read below, instead
    // of a flat constant far stricter than that.
    int    MAXMM  = -1;
    double OV_ERR = 0.08;
    if (const char* s = getenv("ARCS_VODBG_MAXMM")) MAXMM  = atoi(s);
    if (const char* s = getenv("ARCS_VODBG_OVERR")) OV_ERR = atof(s);
    constexpr int BUCKET_CAP = 64;

    // ── Build the two views (forward, reverse-complement) of every read once ────
    std::vector<std::string> rcseq(n);
    for (size_t i = 0; i < n; ++i) rcseq[i] = reverse_complement(reads[i].seq);
    auto view_seq = [&](uint32_t rid, uint8_t view) -> const std::string& {
        return view == 0 ? reads[rid].seq : rcseq[rid];
    };

    // ── APSP candidate table (exact suffix-array overlap discovery) ────────────
    // Replaces the old k-mer-anchor occurrence index entirely: that method
    // could only find an overlap if an exact K0-length window happened to
    // land at the right spot — a real, measured sensitivity ceiling (see
    // build_vodbg_pg header notes). A suffix array finds the TRUE maximal
    // suffix-prefix overlap between every pair of reads directly, no anchor
    // window to miss. This is also why growth below needs no per-step
    // mismatch scanning any more: APSP candidates are exact-match verified by
    // construction (Gusfield's classic all-pairs-suffix-prefix technique).
    int max_cands = 8;
    uint32_t search_cap = 2000;
    if (const char* s = getenv("ARCS_VODBG_MAXCAND")) { int v = atoi(s); if (v >= 1 && v <= 64) max_cands = v; }
    if (const char* s = getenv("ARCS_VODBG_SEARCHCAP")) { long v = atol(s); if (v >= 10) search_cap = (uint32_t)v; }

    const bool VB_TIMING = getenv("ARCS_VODBG_TIMING") != nullptr;
    auto _vb_t0 = std::chrono::steady_clock::now();

    std::vector<std::string> both_views((size_t)n * 2);
    for (size_t i = 0; i < n; ++i) { both_views[2*i] = reads[i].seq; both_views[2*i+1] = rcseq[i]; }
    auto apsp = build_apsp_candidates(both_views, (uint32_t)n, max_cands, (uint32_t)K0, search_cap);
    both_views.clear(); both_views.shrink_to_fit(); // no longer needed; view_seq() reconstructs on demand
    if (VB_TIMING) {
        auto t1 = std::chrono::steady_clock::now();
        fprintf(stderr, "[VB-TIMING] sa_apsp_build: %.2fs\n", std::chrono::duration<double>(t1 - _vb_t0).count());
        _vb_t0 = t1;
    }

    // ── HQ/LQ read division (PgRC2's actual architecture) ───────────────────────
    // Pure exact-match assembly is brittle to real sequencing errors — even one
    // mismatch anywhere breaks an otherwise-long true overlap, since exact
    // matching has zero tolerance by definition. The fix isn't to tolerate
    // mismatches *inside* assembly (tried that — it's exactly what caused the
    // cascading-staleness bug found and fixed above). PgRC2's real answer is
    // architectural: keep assembly purely exact, and simply never feed it
    // error-prone reads. Only HQ reads participate in exact-match growth;
    // every other read maps onto the finished assembly afterward via the
    // mismatch-tolerant fallback pass, which already exists and was correct.
    //
    // Real PgRC2 criterion (confirmed via direct source read, pgrc-params.h /
    // its overlap-search termination logic): a read is HQ iff a trial overlap
    // search finds overlaps >= 65% of read length on BOTH sides (left AND
    // right) — an ACHIEVED-overlap-length test, not a standalone quality-score
    // proxy. We already have the exact answer to that test for free: apsp[e]
    // (e = 2*rid+view) already lists every valid suffix-prefix overlap for
    // that read/view, sorted descending by length — apsp[2*rid+0]'s best
    // entry IS the read's best "left" (predecessor) overlap; apsp[2*rid+1]'s
    // best entry, by the same RC-mirror argument used throughout this file, IS
    // the best "right" (successor) overlap for the read's forward view. No
    // separate trial search needed.
    //
    // The earlier quality-score expected-error criterion (this session's first,
    // weaker-proxy attempt) is kept as an explicit opt-in fallback via
    // ARCS_VODBG_HQ_MAXERR, for A/B comparison against the real PgRC2 rule.
    double hq_ov_frac = 0.65;
    if (const char* s = getenv("ARCS_VODBG_HQ_OVFRAC")) hq_ov_frac = atof(s);
    bool use_quality_criterion = (getenv("ARCS_VODBG_HQ_MAXERR") != nullptr);
    double hq_max_expected_errors = 1.0;
    if (use_quality_criterion) hq_max_expected_errors = atof(getenv("ARCS_VODBG_HQ_MAXERR"));
    // ARCS_VODBG_HQ_ADAPTIVE=1: replace the fixed overlap-fraction threshold with
    // one solved, per-dataset, to hit a TARGET HQ percentage. Real, measured
    // reason this exists: a fixed hq_ov_frac does not generalize across
    // datasets with different coverage depth — 0.45 cut yeast_sub.fq's
    // archive by 7% but GREW E. coli's by 6.3% (same fixed value, opposite
    // effect), because at higher coverage the same absolute overlap-length
    // fraction admits a very different, non-comparable slice of reads into
    // "HQ" (more contention over the same candidate pool, more of them
    // failing to actually merge and falling through to costly approximate
    // fallback anyway). Solving for a target HQ FRACTION instead of a target
    // overlap LENGTH fraction is a coverage-depth-normalized proxy: since
    // best_left/best_right are already fully computed for every read (free,
    // reused from the apsp[] table already built above), re-scanning them at
    // a candidate frac is O(n) and cheap — affordable to binary-search over
    // ~24 iterations for a frac that lands within HQ_ADAPTIVE_TOL of the
    // target, entirely before growth starts.
    bool hq_adaptive = getenv("ARCS_VODBG_HQ_ADAPTIVE") != nullptr;
    double hq_adaptive_target = 0.75;
    if (const char* s = getenv("ARCS_VODBG_HQ_ADAPTIVE_TARGET")) hq_adaptive_target = atof(s);

    std::vector<uint32_t> best_left_v, best_right_v;
    if (!use_quality_criterion) {
        best_left_v.assign(n, 0);
        best_right_v.assign(n, 0);
        for (size_t i = 0; i < n; ++i) {
            for (const auto& c : apsp[2*i + 0]) best_left_v[i]  = std::max(best_left_v[i],  c.overlap);
            for (const auto& c : apsp[2*i + 1]) best_right_v[i] = std::max(best_right_v[i], c.overlap);
        }
    }

    auto hq_fraction_at = [&](double frac) -> double {
        size_t n_hq = 0;
        for (size_t i = 0; i < n; ++i) {
            int thresh = (int)(frac * (double)reads[i].seq.size());
            if ((int)best_left_v[i] >= thresh && (int)best_right_v[i] >= thresh) ++n_hq;
        }
        return (double)n_hq / (double)n;
    };

    if (hq_adaptive && !use_quality_criterion) {
        // Monotonic in frac (higher frac -> fewer or equal HQ reads), so a
        // plain bisection over [0,1] converges directly to whichever frac
        // yields the target HQ fraction for THIS dataset.
        double lo = 0.0, hi = 1.0;
        for (int iter = 0; iter < 24; ++iter) {
            double mid = 0.5 * (lo + hi);
            double got = hq_fraction_at(mid);
            if (got > hq_adaptive_target) lo = mid; else hi = mid; // more HQ than target -> need a stricter (higher) frac
        }
        hq_ov_frac = 0.5 * (lo + hi);
        if (getenv("ARCS_VODBG_EXT_DEBUG"))
            fprintf(stderr, "[HQ-ADAPTIVE] target=%.2f solved_frac=%.4f\n", hq_adaptive_target, hq_ov_frac);
    }

    std::vector<bool> is_hq(n, true);
    if (use_quality_criterion) {
        for (size_t i = 0; i < n; ++i) {
            const std::string& q = reads[i].qual;
            if (q.empty()) continue; // no quality info to classify with — don't penalize
            double expected_errors = 0.0;
            for (char c : q) {
                int qv = (int)(unsigned char)c - 33;
                if (qv < 0) qv = 0;
                expected_errors += std::pow(10.0, -qv / 10.0);
            }
            is_hq[i] = expected_errors <= hq_max_expected_errors;
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            int Lr = (int)reads[i].seq.size();
            int thresh = (int)(hq_ov_frac * Lr);
            is_hq[i] = ((int)best_left_v[i] >= thresh) && ((int)best_right_v[i] >= thresh);
        }
    }
    if (getenv("ARCS_VODBG_EXT_DEBUG")) {
        size_t n_hq = 0; for (size_t i = 0; i < n; ++i) if (is_hq[i]) ++n_hq;
        fprintf(stderr, "[HQ-DBG] hq=%zu (%.1f%%) lq=%zu (%.1f%%) mode=%s frac=%.4f\n",
                n_hq, 100.0*n_hq/n, n-n_hq, 100.0*(n-n_hq)/n,
                use_quality_criterion ? "quality-expected-errors" : (hq_adaptive ? "overlap-length(adaptive)" : "overlap-length(PgRC2-real)"),
                hq_ov_frac);
    }

    // ── Global greedy-overlap contig growth ─────────────────────────────────────
    // True global-greedy (not locally-greedy): a plain "grow one contig fully,
    // then move to the next" walk is provably weaker (a 4-approximation to
    // optimal — Locally Greedy Algorithm) than always taking the single best
    // overlap available *anywhere among all currently-growing contigs* (proven
    // >=2.5-approximation for the real Greedy Algorithm; this is the algorithm
    // class PgRC2 actually uses). Concretely: sequential-per-contig growth lets
    // whichever contig happens to reach a shared candidate *first in file
    // order* claim it, even if a different, still-growing contig would have
    // used it in a much longer, more correct overlap. A global priority queue
    // across all active contigs' current best candidates fixes exactly that.
    //
    // What's deliberately NOT implemented: arbitrary contig-to-contig
    // absorption (two already-multi-read contigs merging into one another).
    // That requires composing position/orientation transforms through
    // potentially many levels of prior absorption — real, correctness-
    // sensitive bookkeeping that's easy to get subtly wrong. Instead, every
    // merge here always consumes a single not-yet-used *original* read (never
    // an already-grown contig), which keeps every position/view computation
    // identical to the already-verified single-contig logic below — global
    // contention is fixed with no added risk to losslessness.
    // Parallelized with the exact same design proven on Method A
    // (kmer_anchor_pg.cpp): each thread grows its OWN independent set of
    // contigs against the ONE shared, already-complete, read-only apsp[]
    // table — no read-set sharding, so no coverage loss. The only shared,
    // contested resource is "which read gets claimed", resolved with a single
    // atomic compare-exchange per read; a contig is only ever touched by the
    // thread that created it, so no per-contig lock is needed either.
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

    int growth_threads = (int)std::thread::hardware_concurrency();
    if (growth_threads < 1) growth_threads = 1;
    if (const char* s = getenv("ARCS_VODBG_GROWTH_THREADS")) { int v = atoi(s); if (v >= 1) growth_threads = v; }
    if ((size_t)growth_threads > n) growth_threads = std::max(1, (int)n);
    std::vector<ThreadState> tstate((size_t)growth_threads);

    // apsp[e] is indexed by e as the *prefix side*: its candidates are reads
    // whose SUFFIX matches e's prefix (valid predecessors of e). Backward
    // growth wants exactly that for the current head — direct lookup.
    // Forward growth wants the opposite (successors of the tail), which
    // isn't what apsp[tail] gives directly. Fix: reverse-complementing both
    // sides of a "D's prefix == T's suffix" relationship turns it into
    // "reverse_complement(D)'s suffix == reverse_complement(T)'s prefix" —
    // i.e. reverse_complement(D) is a valid predecessor of reverse_complement(T).
    // So looking up apsp[reverse_complement(T)] and flipping the candidate's
    // view back gives the true successor D. Same flip undoes the RC framing
    // in the returned view either way, which is also exactly the view the
    // Crc-space append (backward branch below) needs — see its comments.
    // Re-verify APSP's claimed base overlap against the ACTUAL current contig
    // content, not just trust it — this is essential, not optional (see the
    // cascading-staleness bug this fixed, described at length in the
    // session's git history / prior comments here).
    // ── Measurement-only diagnostic (ARCS_VODBG_CONTENTION_DEBUG=1): counts
    // every case where a candidate's own theoretical-best entry (apsp[...][0],
    // sorted descending by overlap) is already claimed, forcing a fallback to
    // a shorter candidate. IMPORTANT, measured caveat: this captures BOTH (a)
    // genuine cross-thread races (fixable, in principle, via cross-thread
    // arbitration) AND (b) the much larger, INHERENT case where two DIFFERENT
    // contigs simply, legitimately want the SAME read and only one can have
    // it — true regardless of threading, present even fully serial. Measured
    // on yeast_sub.fq: this diagnostic reports ~3.48M "overlap bytes lost"
    // (22.6% of pg_len) at 12 threads, but a single-threaded run (zero
    // cross-thread races by construction) still lands within ~12KB (0.08%) of
    // the 12-thread pg_len — meaning almost all of what this counter reports
    // is case (b), not a fixable threading inefficiency. Kept as a diagnostic
    // for future reference, but do not read its raw total as "recoverable
    // bytes" without re-running the single-thread A/B check above.
    // Never changes growth behavior; pure side-channel.
    const bool CONTENTION_DEBUG = getenv("ARCS_VODBG_CONTENTION_DEBUG") != nullptr;
    static std::atomic<uint64_t> g_contention_steps{0};
    static std::atomic<uint64_t> g_total_steps{0};
    static std::atomic<uint64_t> g_overlap_lost{0};

    auto run_thread = [&](int t, size_t range_lo, size_t range_hi) {
        ThreadState& ts = tstate[(size_t)t];
        std::priority_queue<Proposal, std::vector<Proposal>, ProposalLess> heap;

        auto next_candidate = [&](uint32_t cid, uint8_t dir) -> Growth {
            uint32_t rid  = ts.tail_rid[cid][dir];
            uint8_t  view = ts.tail_view[cid][dir];
            uint8_t  query_view = (dir == 0) ? (uint8_t)(1 - view) : view;
            std::string tmp;
            const std::string& S = (dir == 0) ? ts.contigs[cid] : (tmp = reverse_complement(ts.contigs[cid]));
            int Ls = (int)S.size();
            const auto& cand_list = apsp[2*rid + query_view];
            if (CONTENTION_DEBUG && !cand_list.empty()) g_total_steps.fetch_add(1, std::memory_order_relaxed);
            bool first_was_claimed = false;
            uint32_t first_overlap = 0;
            if (CONTENTION_DEBUG && !cand_list.empty()) {
                first_overlap = cand_list[0].overlap;
                first_was_claimed = claimed[cand_list[0].rid].load(std::memory_order_acquire);
            }
            for (const auto& c : cand_list) {
                if (claimed[c.rid].load(std::memory_order_acquire) || !is_hq[c.rid]) continue;
                uint8_t cand_view = (uint8_t)(1 - c.view);
                const std::string& R = view_seq(c.rid, cand_view);
                int L = (int)c.overlap;
                if (L > Ls || L > (int)R.size()) continue;
                if (memcmp(R.data(), S.data() + (Ls - L), (size_t)L) != 0) continue; // stale anchor — skip, try next
                if (CONTENTION_DEBUG && first_was_claimed && c.overlap < first_overlap) {
                    g_contention_steps.fetch_add(1, std::memory_order_relaxed);
                    g_overlap_lost.fetch_add((uint64_t)(first_overlap - c.overlap), std::memory_order_relaxed);
                }
                return {true, c.rid, cand_view, c.overlap};
            }
            return {false, 0, 0, 0};
        };

        auto propose = [&](uint32_t cid, uint8_t dir) {
            Growth g = next_candidate(cid, dir);
            if (g.found) heap.push({(int64_t)g.overlap, cid, dir});
        };

        size_t next_unseeded = range_lo;
        auto seed_next = [&]() -> bool {
            while (next_unseeded < range_hi &&
                   (claimed[next_unseeded].load(std::memory_order_acquire) || !is_hq[next_unseeded] ||
                    reads[next_unseeded].seq.size() < (size_t)K0))
                ++next_unseeded;
            if (next_unseeded >= range_hi) return false;
            uint32_t rid0 = (uint32_t)next_unseeded;
            ++next_unseeded;
            if (!try_claim(rid0)) return true; // lost race; caller loops again, skip logic catches up next call
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

        // No tolerant "extension" here, deliberately — pure exact-match growth
        // only (see build_vodbg_pg header notes for why).
        bool more_to_seed = seed_next();
        while (more_to_seed || !heap.empty()) {
            if (heap.empty()) { more_to_seed = seed_next(); continue; }
            Proposal p = heap.top(); heap.pop();
            Growth g{false, 0, 0, 0};
            while (true) {
                g = next_candidate(p.cid, p.dir); // re-verify fresh: state may have changed since queued
                if (!g.found) break;
                if (!try_claim(g.rid)) continue;  // another thread claimed it — rescan, now excluded
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
                // Backward: mirror via reverse_complement, exactly as the original
                // single-contig bidirectional fix — see build_vodbg_pg header notes.
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

    size_t total_contigs_dbg = 0;
    for (auto& ts : tstate) total_contigs_dbg += ts.contigs.size();
    if (VB_TIMING) {
        auto t1 = std::chrono::steady_clock::now();
        fprintf(stderr, "[VB-TIMING] growth_loop: %.2fs (threads=%d, contigs=%zu)\n",
                std::chrono::duration<double>(t1 - _vb_t0).count(), growth_threads, total_contigs_dbg);
        _vb_t0 = t1;
    }
    if (CONTENTION_DEBUG) {
        uint64_t steps = g_total_steps.load(), csteps = g_contention_steps.load(), lost = g_overlap_lost.load();
        fprintf(stderr, "[CONTENTION-DBG] steps=%llu contention_steps=%llu (%.2f%%) overlap_bytes_lost=%llu (avg %.1f/contention-step)\n",
                (unsigned long long)steps, (unsigned long long)csteps,
                steps ? 100.0 * (double)csteps / (double)steps : 0.0,
                (unsigned long long)lost, csteps ? (double)lost / (double)csteps : 0.0);
    }

    // Merge all threads' contigs into one contig list + global cid remap,
    // exactly as done for Method A (see kmer_anchor_pg.cpp for the identical
    // pattern and its correctness reasoning). Remap place[].cid from
    // thread-local to global FIRST, using each thread's still-intact
    // contig_members, then move everything into the merged lists — doing
    // the remap after the move would iterate already-moved-from vectors.
    std::vector<uint32_t> cid_remap_base((size_t)growth_threads, 0);
    {
        uint32_t running_cid = 0;
        for (int t = 0; t < growth_threads; ++t) {
            cid_remap_base[(size_t)t] = running_cid;
            running_cid += (uint32_t)tstate[(size_t)t].contigs.size();
        }
    }
    for (int t = 0; t < growth_threads; ++t) {
        uint32_t base = cid_remap_base[(size_t)t];
        for (auto& members : tstate[(size_t)t].contig_members)
            for (uint32_t rid : members) place[rid].cid += base;
    }

    std::vector<std::string> contigs;
    std::vector<std::vector<uint32_t>> contig_members;
    contigs.reserve(total_contigs_dbg);
    contig_members.reserve(total_contigs_dbg);
    for (int t = 0; t < growth_threads; ++t) {
        for (auto& c : tstate[(size_t)t].contigs) contigs.push_back(std::move(c));
        for (auto& m : tstate[(size_t)t].contig_members) contig_members.push_back(std::move(m));
    }

    // ── Concatenate contigs -> pg, compute absolute per-contig base offsets ─────
    std::vector<uint32_t> contig_base(contigs.size());
    r.pg.reserve([&]{ size_t t = 0; for (auto& c : contigs) t += c.size(); return t; }());
    for (size_t c = 0; c < contigs.size(); ++c) {
        contig_base[c] = (uint32_t)r.pg.size();
        r.pg += contigs[c];
    }

    // ── Expose placements for in-process reference-free calling (arcs --call) ───
    // Same contract as build_multicontig_pg's call_out (see chain_encoder.h):
    // every original read gets (contig id, contig-local pos, RC flag). Grown
    // reads use place[i] directly (already global post-merge cid/off/view —
    // view IS the RC flag, same convention record_mapped already uses).
    // Fallback-mapped and never-placed reads are filled in further below,
    // once their outcomes are known.
    if (call_out) {
        call_out->contigs = contigs;
        call_out->read_cid.assign(n, 0);
        call_out->read_pos.assign(n, 0);
        call_out->read_rc.assign(n, 0);
        for (size_t i = 0; i < n; ++i) {
            if (!place[i].used) continue;
            call_out->read_cid[i] = place[i].cid;
            call_out->read_pos[i] = place[i].off;
            call_out->read_rc[i]  = place[i].view;
        }
    }

    // ── Emit placements for every read that was resolved during growth ─────────
    // Sorted by absolute pg position (not original read index): pg_pos then
    // delta-codes to almost nothing downstream, same reason the default
    // assembler emits in gpos-sorted order. chain_order carries the mapping
    // back to original index, so this is purely an emission-order choice —
    // lossless either way, this ordering just compresses far better.
    std::vector<uint32_t> resolved;
    resolved.reserve(n);
    for (size_t i = 0; i < n; ++i) if (place[i].used) resolved.push_back((uint32_t)i);
    std::sort(resolved.begin(), resolved.end(), [&](uint32_t a, uint32_t b) {
        uint32_t pa = contig_base[place[a].cid] + place[a].off;
        uint32_t pb = contig_base[place[b].cid] + place[b].off;
        return pa < pb;
    });
    if (getenv("ARCS_VODBG_EXT_DEBUG")) {
        long bad = 0, tiny_mm = 0, huge_mm = 0;
        for (uint32_t i : resolved) {
            uint32_t abs_pos = contig_base[place[i].cid] + place[i].off;
            const std::string& target = view_seq(i, place[i].view);
            if (abs_pos + target.size() > r.pg.size()) { ++bad; continue; }
            int mm = 0;
            for (size_t j = 0; j < target.size(); ++j)
                if (r.pg[abs_pos + j] != target[j] && is_acgt_strict(target[j])) ++mm;
            if (mm > (int)(target.size() / 2)) ++huge_mm; else if (mm > 0) ++tiny_mm;
        }
        fprintf(stderr, "[EXT-DBG] resolved=%zu out_of_bounds=%ld tiny_mm(1..half)=%ld huge_mm(>half)=%ld\n",
                resolved.size(), bad, tiny_mm, huge_mm);
    }
    for (uint32_t i : resolved) {
        uint32_t abs_pos = contig_base[place[i].cid] + place[i].off;
        const std::string& target = view_seq(i, place[i].view);
        record_mapped(r, i, abs_pos, (int)place[i].view, target, (int)target.size());
    }

    // ── Fallback placement pass for reads growth never touched (too short for
    // K0, or would only have extended a contig backward) — search the FINAL,
    // fixed pg directly. Rare in practice; always lossless regardless (falls
    // through to record_append, a plain new reference segment, if no good
    // placement is found).
    //
    // Parallelized identically to Method A's proven index-build fix (see
    // src/flat_kmer_index.h): confirmed via direct reading of PgRC2's real
    // source (ReadsMatchers.cpp's CopMEMReadsApproxMatcher::executeMatching,
    // matching mode 'c', the NORMAL-level default for this exact step —
    // mapping unresolved/LQ reads onto the finished pg) that PgRC2 itself
    // solves this SAME sub-problem the SAME way: build one shared hash index
    // over the pg once (its own genCummMultithreaded two-phase parallel
    // counting-sort — the same shape as our own), then run a plain
    // parallel-for over every unresolved read, each thread independently
    // querying that one shared, already-built, read-only index. This is not
    // a namesake reuse of our own pattern — it's confirmed to be the actual
    // competitor's real architecture for this exact step. record_mapped/
    // record_append themselves stay serial (they mutate shared vectors) —
    // only the search-and-verify work per read is parallelized; results are
    // computed per-thread into a private array, then emitted in original
    // unresolved order afterward, so output is identical to the old serial
    // version's ordering, not just its content. ─────────────────────────────
    std::vector<uint32_t> unresolved;
    for (size_t i = 0; i < n; ++i) if (!place[i].used) unresolved.push_back((uint32_t)i);

    if (!unresolved.empty()) {
        int Lpg = (int)r.pg.size();
        int fallback_threads = (int)std::thread::hardware_concurrency();
        if (fallback_threads < 1) fallback_threads = 1;
        if (const char* s = getenv("ARCS_VODBG_FALLBACK_THREADS")) { int v = atoi(s); if (v >= 1) fallback_threads = v; }

        FlatKmerIndex pgidx;
        if (Lpg >= K0) {
            size_t n_positions = (size_t)(Lpg - K0 + 1);
            pgidx = build_flat_kmer_index_parallel(n_positions, BUCKET_CAP, fallback_threads,
                [&](size_t lo, size_t hi, std::vector<std::pair<uint64_t,uint64_t>>& out) {
                    for (size_t p = lo; p < hi; ++p) {
                        uint64_t km;
                        if (!pack_kmer(r.pg, (int)p, K0, km)) continue;
                        out.push_back({km, (uint64_t)p});
                    }
                });
        }

        struct FallbackResult { int32_t bestpos; uint8_t bestview; int bestmm; };
        std::vector<FallbackResult> results(unresolved.size());

        auto search_one = [&](uint32_t rid) -> FallbackResult {
            int Lr0 = (int)reads[rid].seq.size();
            int maxmm_for_read = (MAXMM >= 0) ? MAXMM : std::max(1, Lr0 / 3);
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
                    pgidx.for_each(km, [&](uint64_t val) {
                        int32_t p = (int32_t)val;
                        int32_t cs = p - off;
                        if (cs < 0 || cs + Lr > Lpg) return;
                        int mm = 0;
                        for (int j = 0; j < Lr && mm < bestmm; ++j)
                            if (r.pg[(size_t)(cs + j)] != R[(size_t)j]) ++mm;
                        if (mm < bestmm) { bestmm = mm; bestpos = cs; bestview = view; }
                    });
                }
            }
            return {bestpos, bestview, bestmm};
        };

        {
            std::vector<std::thread> ths;
            ths.reserve((size_t)fallback_threads);
            for (int t = 0; t < fallback_threads; ++t) {
                size_t lo = unresolved.size() * (size_t)t / (size_t)fallback_threads;
                size_t hi = unresolved.size() * (size_t)(t + 1) / (size_t)fallback_threads;
                ths.emplace_back([&, lo, hi]{
                    for (size_t idx = lo; idx < hi; ++idx) results[idx] = search_one(unresolved[idx]);
                });
            }
            for (auto& th : ths) th.join();
        }

        for (size_t idx = 0; idx < unresolved.size(); ++idx) {
            uint32_t rid = unresolved[idx];
            const FallbackResult& res = results[idx];
            int Lr0 = (int)reads[rid].seq.size();
            int maxmm_for_read = (MAXMM >= 0) ? MAXMM : std::max(1, Lr0 / 3);
            if (res.bestpos >= 0 && res.bestmm <= maxmm_for_read) {
                const std::string& target = view_seq(rid, res.bestview);
                record_mapped(r, rid, (uint32_t)res.bestpos, (int)res.bestview, target, (int)target.size());
                if (call_out) {
                    // Invert absolute pg position back to (contig id, local offset).
                    // upper_bound gives the first base > abs_pos; the owning contig is
                    // the one just before that. NOTE: the fallback match itself only
                    // ever checked against the GLOBAL r.pg boundary (cs + Lr > Lpg),
                    // never against individual contig boundaries — so a match can
                    // legitimately straddle two contigs in r.pg (harmless for
                    // record_mapped, which only needs an absolute position). CallData's
                    // contract requires ONE contig id + a local offset that stays
                    // entirely inside it, so a straddling match can't be represented
                    // this way — found by a dedicated correctness test
                    // (tests/test_vodbg_calldata.cpp), not assumed. Give it its own
                    // singleton contig instead, exactly like a never-placed read: still
                    // a fully valid, self-consistent CallData entry, just not sharing
                    // the growth contig it happened to land in.
                    auto it = std::upper_bound(contig_base.begin(), contig_base.end(), (uint32_t)res.bestpos);
                    size_t cid = (size_t)(it - contig_base.begin()) - 1;
                    uint32_t local_off = (uint32_t)res.bestpos - contig_base[cid];
                    if ((size_t)local_off + target.size() <= contigs[cid].size()) {
                        call_out->read_cid[rid] = (uint32_t)cid;
                        call_out->read_pos[rid] = local_off;
                        call_out->read_rc[rid]  = res.bestview;
                    } else {
                        uint32_t new_cid = (uint32_t)call_out->contigs.size();
                        call_out->contigs.push_back(target);
                        call_out->read_cid[rid] = new_cid;
                        call_out->read_pos[rid] = 0;
                        call_out->read_rc[rid]  = res.bestview;
                    }
                }
            } else {
                record_append(r, rid, reads[rid].seq, (int)reads[rid].seq.size());
                if (call_out) {
                    // Never placed anywhere — becomes its own singleton contig, same
                    // convention build_multicontig_pg uses for a "new contig start" read
                    // (pos=0, rc=0; record_append always stores the read forward, never RC).
                    uint32_t cid = (uint32_t)call_out->contigs.size();
                    call_out->contigs.push_back(reads[rid].seq);
                    call_out->read_cid[rid] = cid;
                    call_out->read_pos[rid] = 0;
                    call_out->read_rc[rid]  = 0;
                }
            }
        }
    }
    if (call_out) call_out->valid = true;

    return r;
}
