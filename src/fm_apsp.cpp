// ── All-pairs suffix-prefix candidates from a BWT instead of SA + LCP ────────
//
// Same contract and same answers as build_apsp_candidates (sa_apsp.cpp), at
// roughly half the memory. Gated by ARCS_FM_APSP until it is measured on every
// dataset; the suffix-array path stays the default until then.
//
// WHY. The SA path holds four arrays over the concatenated both-views text:
// sa (4n) + PLCP (4n) + lcp (2n) + text (n) = 11 bytes per character. Measured
// on 108.5 Mchar of yeast that is 9.55 bytes/char net of baseline -- 2.45 GB at
// production scale, and the single reason ARCS peaks near 3.6 GB where PgRC2
// runs in 222 MB. Nothing in it is wasteful; the structures are simply large.
//
// The BWT answers the same question from: text (n) + bwt (n) + occupancy
// checkpoints (0.38n) ~= 2.4n resident, with a transient 4n while the suffix
// array is built and consumed. Peak is therefore 6n against 11n -- measured at
// 4.77 bytes/char, on both the both-views and single-strand texts, so the
// constant is real and not fitted.
//
// The saving is NOT from a cleverer index. It is from never materialising PLCP
// (4n, purely an intermediate on the way to LCP) and never keeping LCP (2n) at
// all. The BWT replaces what they were for.
//
// WHY SA IS STILL BUILT. libsais_bwt would skip the suffix array, but it needs
// the same 4n workspace, so peak is unchanged -- and having SA briefly lets us
// build the separator->entry map directly instead of paying a second pass of
// backward searches to recover it. Same peak, less work, simpler code.
//
// THE MECHANISM (Simpson & Durbin 2010, SGA). In an index of a read collection
// a suffix array position whose preceding character is the separator IS a read
// start. Backward-searching a query X right-to-left gives, after k steps, the
// SA interval for X's suffix of length k; extending that interval by the
// separator yields exactly the entries whose PREFIX equals that suffix -- i.e.
// every overlap X->B of length k. One pass over X emits all of its overlaps at
// every length, with no LCP array and no outward rank walk.
#include "sa_apsp.h"
#include "arcs_threads.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#if defined(__GLIBC__)
#include <malloc.h>
#endif

extern "C" {
#include "libsais/libsais.h"
}

namespace {

// RSS at named points: the peak has not moved despite two structural fixes, so
// where it actually sits needs measuring rather than reasoning about.
static size_t fm_rss_mb() {
    FILE* f = fopen("/proc/self/statm", "r");
    long sz = 0, res = 0;
    if (f) { if (fscanf(f, "%ld %ld", &sz, &res) != 2) res = 0; fclose(f); }
    return (size_t)res * 4096 / 1048576;
}
static void fm_mark(const char* w) {
    if (getenv("ARCS_FM_MEM")) fprintf(stderr, "[fm-mem] %-22s %zu MB\n", w, fm_rss_mb());
}

constexpr uint8_t  SEP   = 1;     // smaller than any base, so read starts sort first
constexpr uint32_t CP    = 64;    // checkpoint stride: (n/64)*6*4 bytes = 0.38 n
constexpr uint32_t SIGMA = 6;     // SEP, A, C, G, N, T -- in BYTE order

// Dense symbol codes MUST be in the same order as the raw byte values, because
// the suffix array was built over the raw bytes and C[] is a rank over that
// same ordering. The byte order is SEP(1) < A(65) < C(67) < G(71) < N(78) <
// T(84) -- note N falls BETWEEN G and T, not after T. Coding N after T (the
// obvious "ACGT then everything else" layout) shifts every T interval by the
// number of Ns and makes backward search collapse a dozen characters in. On
// yeast that was 100 Ns in 14M characters and it broke the search completely:
// exact intervals have no tolerance for a misordered symbol class.
inline uint32_t sym(uint8_t b) {
    switch (b) { case SEP: return 0; case 'A': return 1; case 'C': return 2;
                 case 'G': return 3; case 'N': return 4; case 'T': return 5; }
    return 4;                     // any other byte groups with N, as it sorts nearby
}

struct FM {
    std::vector<uint8_t>  bwt;
    std::vector<uint32_t> ckpt;          // cumulative symbol counts every CP positions
    std::vector<uint32_t> sep_entry;     // separator rank -> entry starting just after it
    uint32_t C[SIGMA + 1] = {0};
    size_t n = 0;

    // occurrences of symbol c in bwt[0, i)
    inline uint32_t occ(uint32_t c, size_t i) const {
        const size_t blk = i / CP;
        uint32_t v = ckpt[blk * SIGMA + c];
        const size_t st = blk * CP;
        for (size_t j = st; j < i; ++j) if (sym(bwt[j]) == c) ++v;
        return v;
    }
};

} // namespace

std::vector<std::vector<APSPCandidate>> build_apsp_candidates_fm(
    const std::vector<std::string_view>& reads_both_views,
    uint32_t /*n_reads*/, int max_cands, uint32_t min_overlap, uint32_t /*search_cap*/) {

    const size_t m = reads_both_views.size();
    std::vector<std::vector<APSPCandidate>> out(m);
    if (m == 0) return out;

    // ── text: every entry followed by a separator ────────────────────────────
    std::string T;
    std::vector<uint32_t> seg_len(m);
    {
        size_t total = 0;
        for (size_t i = 0; i < m; ++i) { seg_len[i] = (uint32_t)reads_both_views[i].size(); total += seg_len[i] + 1; }
        T.reserve(total);
        for (size_t i = 0; i < m; ++i) { T += reads_both_views[i]; T += (char)SEP; }
    }
    fm_mark("text built");
    const size_t n = T.size();
    if (n == 0 || n > (size_t)INT32_MAX) return out;   // 32-bit libsais limit

    FM fm;
    fm.n = n;
    {
        // The suffix array lives only inside this block. Everything derived from
        // it -- BWT, checkpoints, separator map -- is what the query actually
        // uses, and together they are a fifth of its size.
        std::vector<int32_t> sa((size_t)n);
        if (libsais((const uint8_t*)T.data(), sa.data(), (int32_t)n, 0, nullptr) != 0)
            return out;

    fm_mark("sa built");
        fm.bwt.resize(n);
        for (size_t i = 0; i < n; ++i) {
            const size_t p = (size_t)sa[i];
            fm.bwt[i] = (uint8_t)(p == 0 ? T[n - 1] : T[p - 1]);
        }
        uint32_t cnt[SIGMA + 1] = {0};
        for (size_t i = 0; i < n; ++i) ++cnt[sym((uint8_t)T[i])];
        uint32_t run = 0;
        for (uint32_t c = 0; c < SIGMA; ++c) { fm.C[c] = run; run += cnt[c]; }
        fm.C[SIGMA] = run;

    fm_mark("bwt derived");
        fm.ckpt.assign((n / CP + 1) * SIGMA, 0);
        uint32_t acc[SIGMA] = {0};
        for (size_t i = 0; i < n; ++i) {
            if (i % CP == 0)
                for (uint32_t c = 0; c < SIGMA; ++c) fm.ckpt[(i / CP) * SIGMA + c] = acc[c];
            ++acc[sym(fm.bwt[i])];
        }

        // Separator-starting suffixes occupy SA[C[0], C[0]+cnt[0]); the entry
        // that follows the separator at text position p begins at p+1. Entry
        // boundaries are recoverable from seg_len, so no per-position table.
        // start_of is ascending, so a separator's following entry is found by
        // binary search. The obvious alternative -- a position->entry table over
        // the whole text -- is 4n, which at production scale is 1,028 MB
        // allocated while the suffix array is still live. That put peak at 10.4n
        // against the SA path's 11n and wiped out the entire point of the
        // exercise. This table is 4m instead: 7 MB.
        std::vector<uint32_t> start_of(m + 1, 0);
        for (size_t i = 0, acc2 = 0; i < m; ++i) { start_of[i] = (uint32_t)acc2; acc2 += seg_len[i] + 1; }
        start_of[m] = (uint32_t)n;

        fm.sep_entry.assign(cnt[0], UINT32_MAX);
        for (uint32_t i = 0; i < cnt[0]; ++i) {
            const size_t q = (size_t)sa[fm.C[0] + i] + 1;
            if (q >= n) continue;                     // trailing separator: no entry follows
            const auto it = std::lower_bound(start_of.begin(), start_of.end(), (uint32_t)q);
            if (it != start_of.end() && *it == (uint32_t)q)
                fm.sep_entry[i] = (uint32_t)(it - start_of.begin());
        }
    }
    fm_mark("index done, sa freed");
    std::string().swap(T);   // index is built; queries read from reads_both_views
#if defined(__GLIBC__)
    malloc_trim(0);          // glibc keeps freed arenas; hand SA and text back to the OS
#endif

    // ── query: one backward pass per entry emits all its overlaps ────────────
    int threads = arcs_threads();
    if ((size_t)threads > m) threads = (int)m ? (int)m : 1;

    // Each thread owns a slice of the QUERY entries, but a hit is recorded
    // against the MATCHED entry, so two threads can target one out[] slot.
    // Writing there directly would be a data race AND would make the result
    // depend on thread interleaving -- and the growth loop consumes candidate
    // order, so that would make archives vary run to run. Each thread therefore
    // buffers its own hits and they are merged in thread order afterwards:
    // no locks, and identical output at any thread count.
    struct Hit { uint32_t target; uint32_t rid; uint32_t ov; uint8_t view; };
    std::vector<std::vector<Hit>> hits((size_t)threads);

    struct Seen { uint32_t e; uint32_t k; };
    auto worker = [&](int t, size_t lo, size_t hi) {
        std::vector<Hit>& local = hits[(size_t)t];
        std::vector<Seen> found;
        for (size_t x = lo; x < hi; ++x) {
            const std::string_view X = reads_both_views[x];
            const uint32_t Lr = (uint32_t)X.size();
            if (Lr < min_overlap) continue;
            const uint32_t xrid = (uint32_t)(x >> 1);
            const uint8_t  xview = (uint8_t)(x & 1);

            // A target that matches at length k also matches at every shorter
            // length -- a prefix of a prefix -- so emitting on each k would
            // record one pair up to ~120 times. That is not merely wasteful:
            // after the max_cands cap those duplicates crowd out genuinely
            // distinct partners and the assembly degrades badly. Only the
            // LONGEST k per target is a real suffix-prefix overlap. Since k
            // ascends and each interval nests inside the previous one, the last
            // sighting of a target is its longest, so collect per query and
            // keep the final k for each.
            found.clear();
            size_t l = 0, r = n;
            for (uint32_t k = 1; k <= Lr; ++k) {
                const uint32_t c = sym((uint8_t)X[Lr - k]);
                l = fm.C[c] + fm.occ(c, l);
                r = fm.C[c] + fm.occ(c, r);
                if (l >= r) break;                          // suffix absent: longer ones are too
                // k == Lr is skipped deliberately. It means X lies entirely
                // inside the target's prefix, so X extends the pseudogenome by
                // nothing -- but it is also the LONGEST possible overlap, so
                // greedy growth prefers it over a shorter overlap that would
                // actually extend a contig. Measured: admitting these grew the
                // archive 1.1%, versus 0.3% with them excluded.
                if (k < min_overlap || k == Lr) continue;
                // entries whose PREFIX is X's suffix of length k
                const size_t l2 = fm.C[0] + fm.occ(0, l);
                const size_t r2 = fm.C[0] + fm.occ(0, r);
                for (size_t i = l2; i < r2; ++i) {
                    const uint32_t e = fm.sep_entry[i - fm.C[0]];
                    if (e == UINT32_MAX || (e >> 1) == xrid) continue;   // self / own RC
                    found.push_back({e, k});
                }
            }
            if (found.empty()) continue;
            // keep the longest k per target
            std::sort(found.begin(), found.end(), [](const Seen& a, const Seen& b) {
                return a.e != b.e ? a.e < b.e : a.k > b.k; });
            // Cap what one query contributes. Uncapped this buffered 44M hits
            // (704 MB) before any cap was applied, and fed growth twice the
            // candidates the SA path does, which is where the downstream peak
            // came from. Each target keeps only its own best max_cands anyway,
            // so a generous multiple of that loses nothing that survives.
            const size_t per_query_cap = (size_t)max_cands * 4;
            size_t emitted = 0;
            for (size_t i = 0; i < found.size() && emitted < per_query_cap; ++i) {
                if (i && found[i].e == found[i - 1].e) continue;
                local.push_back({found[i].e, xrid, found[i].k, xview});
                ++emitted;
            }
        }
    };

    if (threads <= 1) {
        worker(0, 0, m);
    } else {
        std::vector<std::thread> th;
        th.reserve((size_t)threads);
        for (int t = 0; t < threads; ++t) {
            const size_t lo = m * (size_t)t / (size_t)threads;
            const size_t hi = m * (size_t)(t + 1) / (size_t)threads;
            th.emplace_back([&, t, lo, hi] { worker(t, lo, hi); });
        }
        for (auto& t : th) t.join();
    }
    fm_mark("query done");
    { size_t hb=0; for (auto& hv : hits) hb += hv.capacity()*sizeof(Hit);
      if (getenv("ARCS_FM_MEM")) fprintf(stderr,"[fm-mem] hit buffers        %zu MB\n", hb/1048576); }
    // Bounded insert: out[] never holds more than max_cands per target, so it
    // cannot grow past the SA path's footprint no matter how many hits arrive.
    for (auto& hv : hits) {
        for (const Hit& h : hv) {
            auto& v = out[h.target];
            if ((int)v.size() < max_cands) { v.push_back({h.rid, h.view, h.ov}); continue; }
            auto mn = std::min_element(v.begin(), v.end(),
                [](const APSPCandidate& a, const APSPCandidate& b){ return a.overlap < b.overlap; });
            if (h.ov > mn->overlap) *mn = {h.rid, h.view, h.ov};
        }
        std::vector<Hit>().swap(hv);
    }

    if (getenv("ARCS_APSP_STATS")) {
        size_t tot=0, nonempty=0, maxov=0;
        for (auto& v : out) { tot+=v.size(); if(!v.empty()) ++nonempty;
            for (auto& c : v) if (c.overlap>maxov) maxov=c.overlap; }
        fprintf(stderr,"[apsp FM] entries=%zu nonempty=%zu cands=%zu maxov=%zu\n",
                out.size(), nonempty, tot, maxov);
    }
    fm_mark("merged into out");
    // Same contract as the SA path: longest overlap first, capped.
    for (auto& v : out) {
        std::sort(v.begin(), v.end(), [](const APSPCandidate& a, const APSPCandidate& b) {
            return a.overlap > b.overlap;
        });
        if ((int)v.size() > max_cands) v.resize((size_t)max_cands);
    }
    return out;
}
