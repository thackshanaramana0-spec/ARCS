#pragma once
// ARCS_REPEAT_ELIM (opt-in): self-referential exact-repeat elimination over an
// already-assembled pseudogenome string.
//
// Why: ARCS's placement pass builds the pg via a single greedy forward walk —
// once a read fails to place well (or the walk has moved on), its sequence is
// appended fresh even if an exact (or near-exact) copy of that stretch already
// exists elsewhere, non-adjacently, in the pg (transposons, multi-copy genes,
// repetitive elements). The placement k-mer index never re-checks "does this
// content already exist anywhere in the finished pg" as a distinct question
// from "does this read overlap here for placement purposes" — this pass asks
// exactly that question, once, after assembly, purely as a compression step.
//
// How (SA/LPF-exact, replacing an earlier copMEM2-*style* sparse-hash
// approach): build ONE suffix array + LCP array over
// T = pg + SEP + reverse_complement(pg) + SEP (the same forward+RC
// concatenation trick sa_apsp.cpp already uses for APSP overlap discovery,
// via the same vendored libsais). For every pg position i, an outward walk
// in SA-rank order from rank_of[i] (Gusfield's technique — running-min LCP
// tracked per direction, same pattern as build_apsp_candidates) finds the
// longest available match to any strictly-earlier position, in EITHER
// orientation: a neighbor in the forward half of T is a plain repeat
// candidate (src < i required — non-overlapping-with-dest is enforced by
// capping the usable length at i - src); a neighbor in the RC half of T at
// offset q corresponds to source range end s_end = n - q (fixed regardless
// of match length — see repeat_elim.cpp for the derivation), so it's usable
// exactly when s_end <= i. This is a real completeness guarantee (any match
// of ANY length >= min_match_len at ANY position is found, not just ones
// that happen to land on a sampled/hashed anchor) — replacing the earlier
// approach's two independent, real sources of missed matches: sparse
// minimizer sampling (a match could exist between two positions that never
// shared a sampled minimizer) and fixed hash-bucket capacity (BUCKET_CAP,
// "keep earliest, skip insert" on overflow silently dropped real
// candidates). Per-position work is embarrassingly parallel (each query only
// reads the shared, already-built SA/LCP/rank_of tables) — computed for all
// positions up front, multithreaded, before the necessarily-serial greedy
// left-to-right parse consumes the precomputed results and jumps ahead by
// each accepted match's length, exactly as before.
//
// Also checks, at every position, whether the current window is the
// reverse-complement of an earlier stretch (palindromic/RC-duplicated
// repeats — inverted transposons, RC-oriented multi-copy elements — which a
// forward-only scan can never find since RC content never equals its own
// forward form byte-for-byte). Both orientations are considered at every
// position; whichever extends longer wins.
//
// A single outward SA-rank walk is bounded (search_cap) as a safety net for
// pathological low-complexity regions (long homopolymer/STR runs where LCP
// stays high across many neighbors) — same rationale, same technique, as
// sa_apsp.cpp's own search_cap. This is not a completeness trade-off in
// practice: the walk already stops the moment the running-min LCP drops
// below min_match_len, which happens almost immediately for real DNA
// content; search_cap only matters for degenerate repeat-dense stretches.
//
// This is a pure, self-contained transform: repeat_elim_decode reproduces the
// exact original pg string bit-for-bit, so nothing downstream (read placement
// offsets, `arcs export`/`coverage`/`query`, variant calling) is aware this
// pass ever ran — it only ever sees the fully-restored pg.

#include <cstdint>
#include <string>
#include <vector>

struct RepeatMatch {
    uint32_t run_len;       // literal bytes immediately preceding this match
    uint32_t src;           // start offset in the *original* pg (always < this match's dest)
    uint32_t len;           // match length (>= min_match_len)
    uint8_t  is_rc = 0;     // 1 = this stretch is reverse_complement(pg[src:src+len)), not a plain copy
};

// Attempts repeat elimination. Returns false (leaving outputs untouched) only
// if pg.size() exceeds what fits in the uint32_t offset fields used on the wire.
bool repeat_elim_encode(const std::string& pg, uint32_t min_match_len,
                         std::string& literal_out,
                         std::vector<RepeatMatch>& matches_out);

// Same, but with an explicit cap on the number of distinct minimizer entries
// the sparse index targets (normally derived internally from pg.size()).
// Exposed only so tests can force the sparse (W>1) code path on small inputs
// without needing a multi-hundred-MB string; production code should use the
// 4-argument overload above.
bool repeat_elim_encode_cap(const std::string& pg, uint32_t min_match_len,
                            uint64_t cap_entries,
                            std::string& literal_out,
                            std::vector<RepeatMatch>& matches_out);

// Inverse of repeat_elim_encode. `matches` must be in the same order
// repeat_elim_encode produced them in (increasing original-text order).
std::string repeat_elim_decode(const std::string& literal,
                                const std::vector<RepeatMatch>& matches,
                                uint32_t orig_len);
