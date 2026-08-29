// ── APSP candidates from a descending merge sweep: no index at all ──────────
//
// Same contract as build_apsp_candidates (sa_apsp.cpp). Gated by
// ARCS_SWEEP_APSP until it is measured on every dataset.
//
// WHY. Every other backend here answers "which reads overlap" by building a
// structure over the read TEXT -- a suffix array (11 bytes/char), a BWT with
// checkpoints (4.77), a rope (0.4). They differ by 20x among themselves and all
// have the wrong shape: they are random-access indexes serving a
// sequential-access algorithm. Greedy growth consumes the longest overlap first
// and never revisits; it never needed "the longest overlap for every pair at
// every length, on demand".
//
// PgRC2 has no suffix array anywhere. GreedySwipingPackedOverlapPseudoGenome-
// Generator sweeps the overlap length downward and merge-joins two sorted index
// lists at each length: reads ordered by prefix against reads ordered by
// suffix-at-offset. Their working set is two uint32 arrays plus packed reads --
// about 8 bytes per READ rather than 11 per character. Measured as a standalone
// prototype (benchmark/reimpl/38_merge_sweep.cpp) that shape holds: 219 MB peak
// against our 3,602 MB, which is PgRC2's real 222 MB almost exactly.
//
// This implementation needs even less than theirs, because it never copies the
// sequences: reads_both_views already points at live storage (Read::seq and the
// RC arena), so the only allocation is two uint32 index arrays -- 14 MB for
// 1.7M entries -- plus the output table every backend has to produce.
//
// SORTED ORDER IS INCREMENTAL, NOT RE-SORTED. Entries sorted by suffix-at-off
// are grouped by the symbol at `off` into contiguous blocks, and within a block
// they are already in suffix-at-(off+1) order. So the next offset is a k-way
// merge of those blocks, O(m) per step, instead of an O(m log m) re-sort. That
// is the difference between PgRC2's 1.5 s and the 110 s a naive re-sorting
// prototype took.
//
// DESCENDING ORDER IS WHAT BOUNDS MEMORY. Sweeping longest-overlap-first means
// the first max_cands candidates found for a target ARE its best ones, so out[]
// is capped from the first iteration and needs no hit buffers -- the rope path
// needed 827 MB of those before it was made to flush. A target whose list is
// full is dropped from the prefix side entirely, so that list shrinks as the
// sweep proceeds, which is PgRC2's readsLeft effect.
#include "sa_apsp.h"
#include "arcs_threads.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <thread>
#include <chrono>

namespace {

// Symbol codes MUST follow BYTE order, because the prefix list is sorted with
// memcmp and the block decomposition has to agree with that ordering. The byte
// order is A(65) < C(67) < G(71) < N(78) < T(84) -- note N sits between G and T,
// not after T. Getting this wrong is silent: it produces plausible counts and
// wrong pairs, which is exactly how it went wrong in fm_apsp.cpp.
constexpr int NSYM = 6;
inline int symcode(unsigned char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2;
                 case 'N': return 3; case 'T': return 4; }
    return 5;
}

// Compare two entries' suffixes starting at `o`, over `len` bytes.
//
// The first byte decides three quarters of DNA comparisons, so checking it
// inline before falling into memcmp removes most of the call overhead -- and
// this runs ~1.3 billion times over a full sweep, so the call overhead was the
// dominant cost rather than the comparing.
inline int suffix_cmp(const char* px, const char* py, uint32_t len) {
    if (*px != *py) return (int)(unsigned char)*px - (int)(unsigned char)*py;
    return len > 1 ? memcmp(px + 1, py + 1, len - 1) : 0;
}

} // namespace

std::vector<std::vector<APSPCandidate>> build_apsp_candidates_sweep(
    const std::vector<std::string_view>& views,
    uint32_t n_reads, int max_cands, uint32_t min_overlap, uint32_t search_cap) {

    const size_t m = views.size();
    std::vector<std::vector<APSPCandidate>> out(m);
    if (m == 0) return out;

    // The sweep is defined over a single overlap length per step, which only
    // makes sense when every entry has the same length -- PgRC2 uses
    // PackedConstantLengthReadsSet for the same reason. Variable-length input
    // (rare outside Illumina) falls back to the suffix array rather than
    // silently assembling differently.
    const uint32_t L = (uint32_t)views[0].size();
    for (size_t i = 1; i < m; ++i)
        if (views[i].size() != L)
            return build_apsp_candidates(views, n_reads, max_cands, min_overlap, search_cap);
    if (L < min_overlap) return out;

    const bool VERIFY = getenv("ARCS_SWEEP_VERIFY") != nullptr;
    size_t bad = 0, good = 0;
    const int nthreads = VERIFY ? 1 : arcs_threads();   // counters are unsynchronised
    double t_join = 0, t_merge = 0;
    const bool TIMING = getenv("ARCS_SWEEP_TRACE") != nullptr;
    auto now = [] { return std::chrono::steady_clock::now(); };

    // ── prefix side: entries ordered by their sequence ───────────────────────
    std::vector<uint32_t> byPrefix(m);
    std::iota(byPrefix.begin(), byPrefix.end(), 0u);
    auto cmp_full = [&](uint32_t a, uint32_t b) {
        const int c = suffix_cmp(views[a].data(), views[b].data(), L);
        return c != 0 ? c < 0 : a < b;
    };
    std::sort(byPrefix.begin(), byPrefix.end(), cmp_full);

    // ── suffix side ──────────────────────────────────────────────────────────
    //
    // THE SWEEP RUNS FROM THE SHORTEST OVERLAP TO THE LONGEST, which is the
    // opposite of PgRC2 and needs explaining.
    //
    // Sorting entries by suffix-at-off is, given them already sorted by
    // suffix-at-(off+1), just a STABLE COUNTING SORT on the single symbol at
    // position off -- O(m), no comparisons at all, since suffix-at-off is
    // symbol[off] followed by suffix-at-(off+1). That only works while `off`
    // DECREASES, i.e. while the overlap length increases.
    //
    // Sweeping the other way (longest overlap first, as PgRC2 does) means each
    // step is a k-way merge of the symbol blocks instead. Measured, that merge
    // was 41.7 s of a 42.9 s sweep and is inherently serial, so it, not the
    // join, was the whole cost.
    //
    // What the descending direction bought was work-bounding: a target whose
    // list was full could be dropped. Measured, that only shrank the prefix side
    // from 1,703,834 to 1,096,941 across 112 of 126 offsets, so it was worth
    // about 20% -- far less than the merge cost. Going upward instead needs an
    // order-INDEPENDENT bounded insert (keep the max_cands largest under a total
    // order) so a candidate arriving late can displace an earlier, shorter one.
    // That also makes the result identical at any thread count.
    const uint32_t off_max = L - min_overlap;
    std::vector<uint32_t> order(m);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        const int c = suffix_cmp(views[a].data() + off_max, views[b].data() + off_max, min_overlap);
        return c != 0 ? c < 0 : a < b;
    });
    std::vector<uint32_t> next_order(m);
    std::vector<uint8_t>  symbuf(m);       // symbol column, extracted once per step

    const size_t pfx_n = m;
    std::vector<uint32_t>& pfx = byPrefix;

    // keep the max_cands largest under a TOTAL order, so arrival order cannot
    // change the kept set
    auto worse = [](const APSPCandidate& a, const APSPCandidate& b) {
        if (a.overlap != b.overlap) return a.overlap < b.overlap;
        if (a.rid     != b.rid)     return a.rid     > b.rid;
        return a.view > b.view;
    };

    for (uint32_t step = 0; step <= off_max; ++step) {
        const uint32_t off = off_max - step;
        const uint32_t ov  = L - off;

        auto j0 = now();
        // Merge-join: suffix(A, off, ov) == prefix(B, 0, ov).
        //
        // Threads partition the PREFIX side, so each writes only to targets in
        // its own range -- no conflicting writes, no hit buffers, and no
        // dependence on interleaving. The matching suffix range is found by
        // binary search, both lists being sorted on the same key.
        //
        // Both sides carry duplicate keys, so the inner scan rescans the equal
        // range from a saved start; a plain two-pointer walk would hide a B from
        // every later A sharing its key (PgRC2's curPreIt serves the same role).
        auto join_range = [&](size_t plo, size_t phi) {
            if (plo >= phi) return;
            auto first_suffix_ge = [&](size_t p) -> size_t {
                if (p >= pfx_n) return m;
                const char* key = views[pfx[p]].data();
                size_t lo = 0, hi = m;
                while (lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    if (suffix_cmp(views[order[mid]].data() + off, key, ov) < 0) lo = mid + 1;
                    else hi = mid;
                }
                return lo;
            };
            size_t ia = first_suffix_ge(plo);
            const size_t ia_end = first_suffix_ge(phi);
            size_t ib = plo;
            while (ia < ia_end && ib < phi) {
                const uint32_t A = order[ia];
                const char* pa = views[A].data() + off;
                while (ib < phi && suffix_cmp(pa, views[pfx[ib]].data(), ov) > 0) ++ib;
                if (ib >= phi) break;
                for (size_t j = ib; j < phi; ++j) {
                    const uint32_t B = pfx[j];
                    if (suffix_cmp(pa, views[B].data(), ov) != 0) break;
                    if ((A >> 1) == (B >> 1)) continue;      // self / own reverse complement
                    if (VERIFY) {
                        if (memcmp(views[A].data() + (L - ov), views[B].data(), ov) == 0) ++good; else ++bad;
                    }
                    const APSPCandidate cand{(uint32_t)(A >> 1), (uint8_t)(A & 1), ov};
                    auto& v = out[B];
                    if ((int)v.size() < max_cands) { v.push_back(cand); continue; }
                    auto mn = std::min_element(v.begin(), v.end(), worse);
                    if (worse(*mn, cand)) *mn = cand;
                }
                ++ia;
            }
        };

        if (nthreads <= 1 || pfx_n < 4096) {
            join_range(0, pfx_n);
        } else {
            std::vector<std::thread> th;
            th.reserve((size_t)nthreads);
            for (int t = 0; t < nthreads; ++t) {
                const size_t plo = pfx_n * (size_t)t / (size_t)nthreads;
                const size_t phi = pfx_n * (size_t)(t + 1) / (size_t)nthreads;
                th.emplace_back([&, plo, phi] { join_range(plo, phi); });
            }
            for (auto& t : th) t.join();
        }
        if (TIMING) t_join += std::chrono::duration<double>(now() - j0).count();

        if (off == 0) break;

        // ── advance to off-1: stable counting sort on the symbol at off-1 ────
        //
        // Two details matter for speed. The symbol is extracted once into a
        // contiguous array rather than read again during the scatter: reading
        // views[order[i]][off-1] is a random hop into scattered sequence
        // storage, and doing it twice per element cost about half the sort.
        // And the sort is parallel -- per-thread counts, prefix-summed by
        // (symbol, thread) so each thread scatters into a disjoint range. That
        // ordering keeps it stable, which the correctness of the whole sweep
        // depends on: the entries within a symbol block must stay in
        // suffix-at-off order.
        auto m0 = now();
        {
            const uint32_t o = off - 1;
            const int T = (nthreads > 1 && m >= 65536) ? nthreads : 1;
            std::vector<uint32_t> cnt((size_t)T * NSYM, 0);
            auto chunk_lo = [&](int t) { return m * (size_t)t / (size_t)T; };

            auto count_chunk = [&](int t) {
                uint32_t* c = &cnt[(size_t)t * NSYM];
                for (size_t i = chunk_lo(t); i < chunk_lo(t + 1); ++i) {
                    const uint8_t sc = (uint8_t)symcode((unsigned char)views[order[i]][o]);
                    symbuf[i] = sc;
                    ++c[sc];
                }
            };
            if (T == 1) count_chunk(0);
            else {
                std::vector<std::thread> th; th.reserve((size_t)T);
                for (int t = 0; t < T; ++t) th.emplace_back(count_chunk, t);
                for (auto& x : th) x.join();
            }
            // prefix sum by (symbol, thread) -- symbol major keeps blocks
            // contiguous, thread minor keeps the sort stable
            std::vector<uint32_t> base((size_t)T * NSYM, 0);
            uint32_t run = 0;
            for (int c2 = 0; c2 < NSYM; ++c2)
                for (int t = 0; t < T; ++t) {
                    base[(size_t)t * NSYM + c2] = run;
                    run += cnt[(size_t)t * NSYM + c2];
                }
            auto scatter_chunk = [&](int t) {
                uint32_t* b = &base[(size_t)t * NSYM];
                for (size_t i = chunk_lo(t); i < chunk_lo(t + 1); ++i)
                    next_order[b[symbuf[i]]++] = order[i];
            };
            if (T == 1) scatter_chunk(0);
            else {
                std::vector<std::thread> th; th.reserve((size_t)T);
                for (int t = 0; t < T; ++t) th.emplace_back(scatter_chunk, t);
                for (auto& x : th) x.join();
            }
            order.swap(next_order);
        }
        if (TIMING) t_merge += std::chrono::duration<double>(now() - m0).count();
    }

    // the contract wants longest first
    for (auto& v : out)
        std::sort(v.begin(), v.end(), [](const APSPCandidate& a, const APSPCandidate& b) {
            return a.overlap > b.overlap;
        });

    if (TIMING) fprintf(stderr, "[sweep] join %.1f s  countsort %.1f s  threads=%d\n", t_join, t_merge, nthreads);
    if (VERIFY) fprintf(stderr, "[sweep-verify] good=%zu bad=%zu\n", good, bad);
    return out;
}
