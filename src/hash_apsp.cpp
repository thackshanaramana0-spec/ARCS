// ── Hash-based all-pairs suffix-prefix candidates ("fast_upgrade") ───────────
//
// Drop-in alternative to build_apsp_candidates (sa_apsp.cpp) with identical
// semantics and output, trading exactness for memory. Selected by
// ARCS_FAST_UPGRADE=1; the suffix-array path remains the default.
//
// WHY THIS EXISTS. The SA path allocates ~13 bytes per character of the
// concatenated read text, and that text holds BOTH strand views of every read
// -- 257 Mchar for 128 Mbase of yeast reads, so ~3.3 GB before the assembly
// proper begins. That is the whole reason ARCS peaks near 4.3 GB where PgRC2
// runs in 232 MB. Nothing about the suffix array is wasteful; the structure
// itself is simply large.
//
// A seed index answers the same question for a fraction of the memory. Every
// entry contributes ONE key (the 32 bases at its start), so the index holds 2n
// postings -- ~1.7M for the same dataset -- instead of a 257M-entry SA plus
// LCP plus rank arrays. Peak drops by roughly the full 3.3 GB.
//
// WHAT IT COSTS. The suffix array finds the MAXIMAL overlap for every pair,
// with no lower bound. A 32-base seed cannot see an overlap shorter than 32,
// and it finds an overlap only when the first 32 bases of the prefix side are
// error-free. So this returns a subset of the SA's candidates: fewer, and
// biased toward long clean overlaps. On assemblies where short overlaps carry
// real signal the pseudogenome will be longer. That is the trade, and it is why
// this is opt-in rather than the default.
//
// Seeds are 32 bases so that a seed packs into exactly one uint64 at 2 bits per
// base. Seed equality is then integer equality -- no hash collisions to verify
// away -- and sliding the seed one base is a shift-or rather than a rescan.
// Walking a read's suffix offsets outward from the longest means the first
// verified hit at each offset is that pair's longest overlap.
#include "arcs_threads.h"
#include "sa_apsp.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <thread>
#include <unordered_map>
#include <vector>
#include <string_view>

namespace {

// Seed width. 32 fills a uint64 exactly, but it also sets the SHORTEST overlap
// this path can see, and that floor is why it returns fewer candidates than the
// suffix array (0.31% archive on yeast). A narrower seed sees shorter overlaps
// at the cost of specificity -- more candidates to verify. ARCS_HASH_SEED lets
// that be measured rather than assumed; prototype work on a standalone
// assembler put the optimum near 16.
static uint32_t seed_bases() {
    static const uint32_t v = [] {
        uint32_t d = 32;
        if (const char* e = getenv("ARCS_HASH_SEED")) { int x = atoi(e); if (x >= 8 && x <= 32) d = (uint32_t)x; }
        return d;
    }();
    return v;
}

inline int base2(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1;
                 case 'G': return 2; case 'T': return 3; }
    return -1;                                       // N or anything else: unseedable
}

// Pack the 32 bases at p. Fails (false) on any non-ACGT byte, which is correct
// rather than merely convenient: a seed spanning an N would compare equal to
// some real 32-mer under any 2-bit encoding.
inline bool pack_seed(const char* p, uint32_t w, uint64_t& out) {
    uint64_t k = 0;
    for (uint32_t i = 0; i < w; ++i) {
        const int v = base2(p[i]);
        if (v < 0) return false;
        k = (k << 2) | (uint64_t)v;
    }
    out = k;
    return true;
}

struct Found { uint32_t e; uint32_t rid; uint32_t overlap; uint8_t view; };

} // namespace

std::vector<std::vector<APSPCandidate>> build_apsp_candidates_hash(
    const std::vector<std::string_view>& reads_both_views,
    uint32_t n_reads, int max_cands, uint32_t min_overlap, uint32_t /*search_cap*/) {

    const size_t m = reads_both_views.size();        // 2 * n_reads
    std::vector<std::vector<APSPCandidate>> out(m);
    if (m == 0) return out;

    // A seed shorter than the seed width can never be found, so raise the floor
    // rather than silently returning nothing for overlaps below it.
    const uint32_t SEED_BASES = seed_bases();
    const uint64_t SEED_MASK = (SEED_BASES >= 32) ? ~0ULL : ((1ULL << (2 * SEED_BASES)) - 1);
    const uint32_t min_ov = std::max(min_overlap, SEED_BASES);

    // ── index: first SEED_BASES of every entry -> entry ids ──────────────────
    std::unordered_map<uint64_t, std::vector<uint32_t>> prefix_idx;
    prefix_idx.reserve(m * 2);
    for (uint32_t e = 0; e < (uint32_t)m; ++e) {
        const std::string_view E = reads_both_views[e];
        if (E.size() < SEED_BASES) continue;
        uint64_t k;
        if (pack_seed(E.data(), SEED_BASES, k)) prefix_idx[k].push_back(e);
    }

    // ── for each entry x, find entries whose PREFIX matches a SUFFIX of x ────
    // out[] is indexed by the PREFIX side (see sa_apsp.h), so a hit at offset
    // `off` in x is recorded against the matched entry e, carrying x's read and
    // view as the candidate. Entries are laid out as 2*rid + view.
    int threads = arcs_threads();
    if (threads < 1) threads = 1;
    if (const char* s = getenv("ARCS_VODBG_APSP_THREADS")) { int v = atoi(s); if (v >= 1) threads = v; }
    if ((size_t)threads > m) threads = std::max(1, (int)m);

    // Each thread collects its own hits; they are merged in thread order below
    // so the result does not depend on thread scheduling.
    std::vector<std::vector<Found>> hits((size_t)threads);

    auto worker = [&](int t, size_t lo, size_t hi) {
        std::vector<Found>& local = hits[(size_t)t];
        for (size_t x = lo; x < hi; ++x) {
            const std::string_view X = reads_both_views[x];
            if (X.size() < min_ov) continue;
            const uint32_t xrid  = (uint32_t)(x >> 1);
            const uint8_t  xview = (uint8_t)(x & 1);
            const uint32_t maxoff = (uint32_t)X.size() - min_ov;

            int kept = 0;
            uint64_t k = 0;
            bool rolling = false;
            for (uint32_t off = 1; off <= maxoff; ++off) {
                if (off + SEED_BASES > X.size()) break;
                if (!rolling) {
                    if (!pack_seed(X.data() + off, SEED_BASES, k)) continue;   // N inside the window
                    rolling = true;
                } else {
                    const int v = base2(X[off + SEED_BASES - 1]);
                    if (v < 0) { rolling = false; continue; }
                    k = ((k << 2) | (uint64_t)v) & SEED_MASK;
                }
                auto it = prefix_idx.find(k);
                if (it == prefix_idx.end()) continue;

                const uint32_t ov = (uint32_t)X.size() - off;      // x's suffix length
                for (uint32_t e : it->second) {
                    if ((e >> 1) == xrid) continue;                // no self / own-RC overlap
                    const std::string_view E = reads_both_views[e];
                    if (E.size() < ov) continue;
                    if (memcmp(X.data() + off, E.data(), ov) != 0) continue;
                    local.push_back({e, xrid, ov, xview});
                    if (++kept >= max_cands) break;
                }
                if (kept >= max_cands) break;
            }
        }
    };

    {
        std::vector<std::thread> th;
        th.reserve((size_t)threads);
        for (int t = 0; t < threads; ++t) {
            const size_t lo = m * (size_t)t / (size_t)threads;
            const size_t hi = m * (size_t)(t + 1) / (size_t)threads;
            th.emplace_back([&, t, lo, hi] { worker(t, lo, hi); });
        }
        for (auto& h : th) h.join();
    }

    for (const auto& bucket : hits)
        for (const Found& f : bucket)
            out[f.e].push_back({f.rid, f.view, f.overlap});

    // Same contract as the SA path: longest overlap first, one entry per
    // distinct read, capped at max_cands.
    for (auto& v : out) {
        std::sort(v.begin(), v.end(), [](const APSPCandidate& a, const APSPCandidate& b) {
            if (a.overlap != b.overlap) return a.overlap > b.overlap;
            if (a.rid     != b.rid)     return a.rid     < b.rid;      // stable tie-break
            return a.view < b.view;
        });
        std::vector<APSPCandidate> kept;
        kept.reserve(std::min((size_t)max_cands, v.size()));
        for (const APSPCandidate& c : v) {
            bool dup = false;
            for (const APSPCandidate& k : kept) if (k.rid == c.rid) { dup = true; break; }
            if (dup) continue;
            kept.push_back(c);
            if ((int)kept.size() >= max_cands) break;
        }
        v.swap(kept);
    }

    if (getenv("ARCS_VODBG_TIMING")) {
        size_t total = 0;
        for (const auto& v : out) total += v.size();
        fprintf(stderr, "[VB-TIMING]   hash_apsp: entries=%zu index_keys=%zu candidates=%zu\n",
                m, prefix_idx.size(), total);
    }
    (void)n_reads;
    return out;
}
