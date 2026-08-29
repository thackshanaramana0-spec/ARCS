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

    // ── prefix side: entries ordered by their sequence ───────────────────────
    std::vector<uint32_t> byPrefix(m);
    std::iota(byPrefix.begin(), byPrefix.end(), 0u);
    auto cmp_full = [&](uint32_t a, uint32_t b) {
        const int c = suffix_cmp(views[a].data(), views[b].data(), L);
        return c != 0 ? c < 0 : a < b;
    };
    std::sort(byPrefix.begin(), byPrefix.end(), cmp_full);

    // ── suffix side: same order at offset 0, advanced one symbol per step ────
    std::vector<uint32_t> order(byPrefix);
    std::vector<uint32_t> next_order(m);

    size_t pfx_n = m;                       // live prefix entries (shrinks)
    std::vector<uint32_t> pfx(byPrefix);    // compacted working copy

    for (uint32_t off = 0; off + min_overlap <= L; ++off) {
        const uint32_t ov = L - off;

        // merge-join: suffix(A, off, ov) == prefix(B, 0, ov)
        //
        // Both sides carry duplicate keys, so this cannot be a plain two-pointer
        // walk: advancing past a B that is already full would hide it from every
        // later A sharing the key. The equal range is rescanned from a saved
        // start, which is what PgRC2's curPreIt is for.
        size_t ia = 0, ib = 0;
        while (ia < m && ib < pfx_n) {
            const uint32_t A = order[ia];
            const char* pa = views[A].data() + off;
            while (ib < pfx_n && suffix_cmp(pa, views[pfx[ib]].data(), ov) > 0) ++ib;
            if (ib >= pfx_n) break;
            for (size_t j = ib; j < pfx_n; ++j) {
                const uint32_t B = pfx[j];
                if (suffix_cmp(pa, views[B].data(), ov) != 0) break;
                if ((A >> 1) == (B >> 1)) continue;          // self / own reverse complement
                auto& v = out[B];
                if ((int)v.size() >= max_cands) continue;
                if (VERIFY) {
                    // the claim: A's suffix of length ov IS B's prefix
                    if (memcmp(views[A].data() + (L - ov), views[B].data(), ov) == 0) ++good; else ++bad;
                }
                // descending sweep: this is already the best remaining overlap
                // for B, so the list stays sorted with no final sort
                v.push_back({(uint32_t)(A >> 1), (uint8_t)(A & 1), ov});
            }
            ++ia;
        }

        // drop targets whose lists are full: they cannot be improved by any
        // shorter overlap, and removing them shrinks the join for every later
        // step (PgRC2's readsLeft)
        size_t w = 0;
        for (size_t i = 0; i < pfx_n; ++i)
            if ((int)out[pfx[i]].size() < max_cands) pfx[w++] = pfx[i];
        pfx_n = w;
        if (getenv("ARCS_SWEEP_TRACE") && (off % 16 == 0 || pfx_n == 0))
            fprintf(stderr, "[sweep] off=%u ov=%u prefix_live=%zu\n", off, ov, pfx_n);
        if (pfx_n == 0) break;

        if (off + 1 + min_overlap > L) break;

        // ── advance the suffix order from `off` to `off+1` ───────────────────
        // order is sorted by suffix-at-off, so it is grouped by symbol at `off`
        // into contiguous blocks, and each block is internally in
        // suffix-at-(off+1) order. A k-way merge of the blocks therefore yields
        // suffix-at-(off+1) order in O(m) rather than O(m log m).
        uint32_t bstart[NSYM + 1];
        {
            size_t i = 0;
            for (int s = 0; s < NSYM; ++s) {
                bstart[s] = (uint32_t)i;
                while (i < m && symcode((unsigned char)views[order[i]][off]) == s) ++i;
            }
            bstart[NSYM] = (uint32_t)i;
            // any entries left over would mean order was not grouped as assumed
            if (i != m) {
                std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
                    const int c = suffix_cmp(views[a].data() + off + 1, views[b].data() + off + 1, L - off - 1);
                    return c != 0 ? c < 0 : a < b;
                });
                continue;
            }
        }
        uint32_t head[NSYM], end[NSYM];
        for (int s = 0; s < NSYM; ++s) { head[s] = bstart[s]; end[s] = bstart[s + 1]; }
        const uint32_t noff = off + 1;
        const uint32_t nlen = L - noff;
        for (size_t k = 0; k < m; ++k) {
            int best = -1;
            for (int s = 0; s < NSYM; ++s) {
                if (head[s] >= end[s]) continue;
                if (best < 0) { best = s; continue; }
                const uint32_t x = order[head[s]], y = order[head[best]];
                const int c = suffix_cmp(views[x].data() + noff, views[y].data() + noff, nlen);
                if (c < 0 || (c == 0 && x < y)) best = s;
            }
            next_order[k] = order[head[best]++];
        }
        order.swap(next_order);
    }

    if (VERIFY) fprintf(stderr, "[sweep-verify] good=%zu bad=%zu\n", good, bad);
    return out;
}
