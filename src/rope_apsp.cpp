// ── APSP candidates from an incrementally-built rope BWT ─────────────────────
//
// Same contract as build_apsp_candidates (sa_apsp.cpp), at a twenty-fourth of
// the index memory. Gated by ARCS_ROPE_APSP until measured on every dataset.
//
// WHY. The suffix-array path holds sa(4n) + PLCP(4n) + lcp(2n) + text(n) = 11
// bytes per character of the both-views text: 2,730 MB at production scale, and
// the single reason ARCS peaks at 3.6 GB. Our own BWT path (fm_apsp.cpp) got
// that to 4.77 bytes/char by never building PLCP, but it still pays 4n of
// suffix-array workspace in one shot.
//
// The rope (ropebwt2's B+-tree, vendored from ropebwt3 under MIT) builds the
// BWT by INSERTING sequences one batch at a time, so peak tracks the compressed
// index rather than any multiple of n. Measured on yeast_sub, 1,702,550 entries
// over 257 Mchar: the index is 113 MB and builds in 11.6 s, against 2,730 MB
// and 12.8 s. Same answers, 24x less memory, no slower to build.
//
// THE COST is query speed: rank2a on a run-length-compressed rope is ~342 ns
// against ~90 ns for a flat checkpoint array, and the walk needs two per
// position. That is the trade this path makes, and why it is opt-in.
//
// TWO THINGS THAT ARE EASY TO GET WRONG, both found by verifying candidates
// against the actual strings rather than trusting counts:
//
//   1. mr_insert_multi wants each sequence REVERSED and NUL-terminated. The
//      resulting index is over the forward sequences.
//   2. The sentinel block of F is ordered by the sequence FOLLOWING the
//      sentinel -- i.e. lexicographically -- NOT by insertion order. Assuming
//      insertion order made every single candidate wrong (good=4, bad=247,022)
//      while still producing plausible-looking counts. ropebwt3 keeps an r2i[]
//      table for this; sorting the entries recovers the same mapping without
//      vendoring their sampled-suffix-array layer.
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
#include <mutex>
#if defined(__GLIBC__)
#include <malloc.h>
#endif

extern "C" {
#include "ropebwt/mrope.h"
}

namespace {

// rope alphabet: 0 = sentinel, 1..5 = A C G T N (their ordering, and the index
// is built consistently with it, so it defines the collation)
inline uint8_t nt6(char c) {
    switch (c) { case 'A': return 1; case 'C': return 2;
                 case 'G': return 3; case 'T': return 4; }
    return 5;
}

static size_t rope_rss_mb() {
    FILE* f = fopen("/proc/self/statm", "r"); long sz = 0, res = 0;
    if (f) { if (fscanf(f, "%ld %ld", &sz, &res) != 2) res = 0; fclose(f); }
    return (size_t)res * 4096 / 1048576;
}
static void rmark(const char* w) {
    if (getenv("ARCS_ROPE_MEM")) fprintf(stderr, "[rope-mem] %-20s %zu MB\n", w, rope_rss_mb());
}

} // namespace

std::vector<std::vector<APSPCandidate>> build_apsp_candidates_rope(
    const std::vector<std::string_view>& reads_both_views,
    uint32_t /*n_reads*/, int max_cands, uint32_t min_overlap, uint32_t /*search_cap*/) {

    const size_t m = reads_both_views.size();
    std::vector<std::vector<APSPCandidate>> out(m);
    if (m == 0) return out;

    // ── build: insert reversed sequences in batches ──────────────────────────
    // Batch size bounds the transient input buffer; the rope itself grows with
    // the compressed index, which is the whole point.
    size_t batch_bytes = 32u << 20;
    if (const char* e = getenv("ARCS_ROPE_BATCH_MB")) {
        long long v = atoll(e); if (v > 0) batch_bytes = (size_t)v << 20;
    }
    rmark("entry");
    mrope_t* mr = mr_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN, MR_SO_IO);
    {
        std::vector<uint8_t> buf;
        buf.reserve(batch_bytes + 512);
        for (size_t i = 0; i < m; ++i) {
            const std::string_view S = reads_both_views[i];
            for (size_t k = S.size(); k-- > 0; ) buf.push_back(nt6(S[k]));
            buf.push_back(0);
            if (buf.size() >= batch_bytes || i + 1 == m) {
                mr_insert_multi(mr, (int64_t)buf.size(), buf.data(), 1);
                buf.clear();
            }
        }
    }

    rmark("rope built");
    // C[]: cumulative symbol counts, read straight off the full-range rank
    int64_t C[7] = {0};
    {
        int64_t ok[6], ol[6];
        mr_rank2a(mr, 0, mr_get_tot(mr), ok, ol);
        int64_t run = 0;
        for (int c = 0; c < 6; ++c) { C[c] = run; run += ol[c]; }
        C[6] = run;
    }

    rmark("C[] done");
    // sentinel rank -> entry (see note 2 in the header). The order must be the
    // ROPE's collation, not the bytes': nt6 puts N last (A<C<G<T<N) whereas the
    // raw bytes put it between G and T ('N'=78 < 'T'=84). Sorting by bytes
    // misplaces every entry whose comparison reaches an N, which on real data
    // cost 12% of archive while still verifying clean on an N-free sample --
    // the standalone spike filtered N-containing reads, so it could not see it.
    auto nt6_less = [&](uint32_t a, uint32_t b) {
        const std::string_view A = reads_both_views[a], B = reads_both_views[b];
        const size_t n = std::min(A.size(), B.size());
        for (size_t i = 0; i < n; ++i) {
            const uint8_t ca = nt6(A[i]), cb = nt6(B[i]);
            if (ca != cb) return ca < cb;
        }
        return A.size() < B.size();
    };
    std::vector<uint32_t> sent2ent(m);
    for (uint32_t i = 0; i < m; ++i) sent2ent[i] = i;
    std::stable_sort(sent2ent.begin(), sent2ent.end(), nt6_less);

#if defined(__GLIBC__)
    malloc_trim(0);
#endif

    rmark("sent2ent done");
    // ── query: one backward pass per entry, sentinel step at each length ─────
    int threads = arcs_threads();
    if ((size_t)threads > m) threads = (int)m ? (int)m : 1;

    struct Hit { uint32_t target; uint32_t rid; uint32_t ov; uint8_t view; };
    struct Seen { uint32_t e; uint32_t k; };

    // Accumulating every hit before capping cost 827 MB -- 31M hits, and the
    // vector growth overshoots on top. Threads flush into out[] periodically
    // instead, under a lock that is contended only once per few hundred
    // thousand hits.
    //
    // Flushing makes insertion order depend on thread timing, so the bounded
    // insert must be order-INDEPENDENT or archives would vary run to run. It is:
    // the kept set is the max_cands largest under a TOTAL order (overlap, then
    // rid, then view), which is the same set whatever sequence they arrive in.
    auto worse = [](const APSPCandidate& a, const APSPCandidate& b) {
        if (a.overlap != b.overlap) return a.overlap < b.overlap;
        if (a.rid     != b.rid)     return a.rid     > b.rid;
        return a.view > b.view;
    };
    std::mutex out_mu;
    auto flush = [&](std::vector<Hit>& local) {
        std::lock_guard<std::mutex> lk(out_mu);
        for (const Hit& h : local) {
            auto& v = out[h.target];
            const APSPCandidate cand{h.rid, h.view, h.ov};
            if ((int)v.size() < max_cands) { v.push_back(cand); continue; }
            auto mn = std::min_element(v.begin(), v.end(), worse);
            if (worse(*mn, cand)) *mn = cand;
        }
        local.clear();
    };
    size_t flush_at = 1u << 18;                       // ~256k hits: 4 MB per thread
    if (const char* e = getenv("ARCS_ROPE_FLUSH")) { long long v = atoll(e); if (v > 0) flush_at = (size_t)v; }

    auto worker = [&](size_t lo, size_t hi) {
        std::vector<Hit> local;
        local.reserve(flush_at + 1024);
        std::vector<Seen> found;
        for (size_t x = lo; x < hi; ++x) {
            const std::string_view X = reads_both_views[x];
            const uint32_t Lr = (uint32_t)X.size();
            if (Lr < min_overlap) continue;
            const uint32_t xrid = (uint32_t)(x >> 1);
            const uint8_t  xview = (uint8_t)(x & 1);

            found.clear();
            int64_t k = 0, l = C[6];
            for (uint32_t i = 0; i < Lr; ++i) {
                const uint8_t c = nt6(X[Lr - 1 - i]);
                int64_t ok[6], ol[6];
                mr_rank2a(mr, k, l, ok, ol);
                k = C[c] + ok[c];
                l = C[c] + ol[c];
                if (l <= k) break;
                const uint32_t klen = i + 1;
                if (klen < min_overlap) continue;
                // a single occurrence is X's own suffix: nothing can start with
                // it, and no longer length can either
                if (l - k == 1 && klen < Lr) break;
                int64_t sk[6], sl[6];
                mr_rank2a(mr, k, l, sk, sl);
                for (int64_t r = C[0] + sk[0]; r < C[0] + sl[0]; ++r) {
                    if (r < 0 || (size_t)r >= m) continue;
                    const uint32_t e = sent2ent[(size_t)r];
                    if ((e >> 1) == xrid) continue;             // self / own RC
                    found.push_back({e, klen});
                }
            }
            if (found.empty()) continue;
            // longest length per target -- a target matching at length k also
            // matches at every shorter one
            std::sort(found.begin(), found.end(), [](const Seen& a, const Seen& b) {
                return a.e != b.e ? a.e < b.e : a.k > b.k; });
            for (size_t i = 0; i < found.size(); ++i) {
                if (i && found[i].e == found[i - 1].e) continue;
                local.push_back({found[i].e, xrid, found[i].k, xview});
            }
            if (local.size() >= flush_at) flush(local);
        }
        if (!local.empty()) flush(local);
    };

    if (threads <= 1) {
        worker(0, m);
    } else {
        std::vector<std::thread> th;
        th.reserve((size_t)threads);
        for (int t = 0; t < threads; ++t) {
            const size_t lo = m * (size_t)t / (size_t)threads;
            const size_t hi = m * (size_t)(t + 1) / (size_t)threads;
            th.emplace_back([&, lo, hi] { worker(lo, hi); });
        }
        for (auto& t : th) t.join();
    }
    rmark("query done");
    mr_destroy(mr);

    rmark("rope freed");
    rmark("query+merge done");
    for (auto& v : out)
        std::sort(v.begin(), v.end(), [](const APSPCandidate& a, const APSPCandidate& b) {
            return a.overlap > b.overlap; });
    return out;
}
