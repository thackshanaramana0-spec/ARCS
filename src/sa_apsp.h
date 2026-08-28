#pragma once
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

// ── Suffix array + LCP array (Kasai) ─────────────────────────────────────────
// Standard O(n log n) prefix-doubling construction (Karp–Miller–Rosenberg
// ranking) using radix/counting sort per round — not std::sort-based
// comparison sort, which would cost an extra log n factor and isn't
// tractable at real dataset scale (hundreds of millions of characters).
struct SuffixArray {
    // Storage stays uint32_t deliberately (not widened to uint64_t) even
    // though construction internally uses libsais64 for texts beyond
    // INT32_MAX chars: real position/length values here never exceed
    // ~4.29B (UINT32_MAX) for any realistic dataset this project processes
    // (confirmed by direct memory-cost analysis: widening to uint64_t would
    // roughly double SA+LCP's already-dominant share of peak RAM for no
    // realistic benefit). build_libsais enforces this ceiling explicitly
    // (throws rather than silently truncating) if ever actually exceeded.
    std::vector<uint32_t> sa;   // sa[i] = start position of the i-th suffix in sorted order
    // uint16, not uint32. Every LCP here is a match between two read suffixes
    // in a text whose reads are separated by a byte no base can equal, so the
    // values live in read-length territory (~150) while a uint32 reserves four
    // bytes for each. At 257 Mchar that halved array is ~490 MB off the
    // construction peak. Values are saturated at 65535 on the way in, which is
    // lossless for every use here: overlaps are capped at the read length, and
    // chain-pg already refuses reads longer than 65535 bp.
    std::vector<uint16_t> lcp; // lcp[i] = LCP(suffix at sa[i-1], suffix at sa[i]); lcp[0] = 0

    void build(const std::string& text);

    // OpenMP-parallel construction via libsais (vendored, third_party/libsais):
    // O(n) induced-sorting (SA-IS), not O(n log n) prefix-doubling — faster
    // even single-threaded, plus real parallel SA *and* LCP construction
    // (PLCP-based, not Kasai's — Kasai's has an inherent sequential
    // dependency; PLCP does not). Confirmed via direct empirical testing
    // (test_libsais_apsp.cpp) that although libsais's raw (sa,lcp) arrays use
    // a different tie-breaking convention than build() above at
    // separator-adjacent positions, the resulting APSP overlap candidates
    // computed from them are IDENTICAL to brute-force-verified correct
    // results — the tie-break difference doesn't affect correctness, only
    // build() and build_libsais() aren't byte-identical to each other.
    // threads=0 lets libsais pick its own OpenMP default.
    //
    // Transparently dispatches to libsais64 (int64_t-indexed) when text.size()
    // exceeds INT32_MAX — found via a real crash on real GIAB HG002 data
    // (~3.8B chars): the 32-bit libsais API's own int32_t indices overflow
    // there, and — more importantly — this SAME function's own prior version
    // computed `n` via a premature signed-int cast, which silently wrapped
    // negative before ever reaching a size check, producing the exact
    // "cannot create std::vector larger than max_size()" crash observed.
    // Fixed at the root: size is now computed via size_t/int64_t throughout,
    // and the 32-vs-64-bit libsais choice is made from that safe value, never
    // from an already-overflowed signed int.
    void build_libsais(const std::string& text, int threads);
};

// ── All-Pairs Suffix-Prefix (APSP) candidate table ───────────────────────────
// For every (read, view) "prefix side", finds the best-overlapping (read,
// view) "suffix side" partners — i.e. entries whose sequence, taken as a
// whole, ends exactly where the shared match ends (a genuine suffix-prefix
// overlap, not just any shared substring) — using the suffix array's LCP
// structure directly, via Gusfield's classic technique: walk outward from a
// read-start's rank in SA order, tracking the running-minimum LCP, and check
// at each step whether that running LCP reaches exactly the end of the
// neighboring suffix's owning read. This finds the TRUE maximal overlap (no
// fixed anchor-length blind spot, unlike k-mer seeding), which is what
// actually gives PgRC2-class assembly its edge over anchor-based methods.
struct APSPCandidate { uint32_t rid; uint8_t view; uint32_t overlap; };

// reads_both_views[i] for i in [0, 2*n): forward and RC sequence of read i/2
// (interleaved: index 2*rid + view). Returns, for each such entry, up to
// `max_cands` best distinct-read overlap partners (sorted by overlap desc).
std::vector<std::vector<APSPCandidate>> build_apsp_candidates(
    const std::vector<std::string_view>& reads_both_views,
    uint32_t n_reads, int max_cands, uint32_t min_overlap, uint32_t search_cap);

// Same contract, seed-hash implementation instead of a suffix array: far less
// memory (2n index postings vs a 257 Mchar SA + LCP + rank on yeast), but it
// cannot see overlaps below its 32-base seed width and needs the prefix side's
// first 32 bases error-free, so it returns a SUBSET of the above. Opt-in via
// ARCS_FAST_UPGRADE=1. See hash_apsp.cpp.
std::vector<std::vector<APSPCandidate>> build_apsp_candidates_hash(
    const std::vector<std::string_view>& reads_both_views,
    uint32_t n_reads, int max_cands, uint32_t min_overlap, uint32_t search_cap);
