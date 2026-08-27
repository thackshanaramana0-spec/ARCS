#pragma once
// ARCS-DNA: pseudogenome-native DNA compressor.
// Uses pre-seeded FCM (orders 6 + 12), RC-aware context mixing,
// adaptive log-loss weights, and rANS coding.
// No subprocess. No temp files. Embedded directly in ARCS.
//
// Genuinely different from GeCo3:
//   - Pre-seeding from reads (GeCo3 is blind to reads)
//   - rANS coder (GeCo3 uses range/arithmetic)
//   - Adaptive log-loss mixing (GeCo3 uses gamma)
//   - Sparse order-12 table (GeCo3 uses fixed arrays)
//   - Chain-start structural priors
#include <string>
#include <vector>
#include <cstdint>

// Compress a raw DNA pseudogenome string.
// reads:        original reads (used to pre-seed context model)
// chain_order:  chain_order[i] = original read index at chain position i
// pg_pos:       pg_pos[i] = start position in pg of chain position i
// seed:         chars from the preceding block's pg tail. Warms both the
//               context shift-registers AND the FCM frequency tables, so
//               split-block coding starts in the correct state (lossless
//               parallel-block support).
//
// NEGATIVE RESULT — do not re-attempt block-parallel coding for speed.
// The coder is single-threaded and dominates both sides (~8.7s encode, ~8.1s
// decode on a 16.3 MB yeast pg; 72% of decompress). Splitting it into blocks
// parallelises almost perfectly but costs ratio, and seeding does not buy the
// ratio back. Measured on that real pg (bytes, whole-pg = 3,432,287):
//
//     N= 2 cold                     3,573,136   +140,849
//     N= 4 cold                     3,686,154   +253,867
//     N= 8 cold                     3,795,467   +363,180
//     N=12 cold                     3,823,581   +391,294
//     N= 4, 256 KB shared seed      3,663,445   +231,158
//     N= 4,   4 MB shared seed      3,628,994   +196,707
//
// The seed tested is the strongest one a DECODER can actually have: block 0 is
// decoded first and every other block seeds from its tail, so phase 2 still
// runs fully parallel. Even 4 MB of it recovers only 22% of the 4-way loss --
// the high-order FCM tables warm too slowly relative to block size -- while
// adding 1.1s of seeding work back.
//
// The trade is refused on size, not on effort: 4-way costs 196,707 B against
// margins over PgRC2 of 402,417 B (yeast), 421,806 B (E. coli) and 150,884 B
// (P. falciparum) -- half the first, most of the last. Every dataset has to
// keep winning, so a split that erases one margin is not available at any
// speed. Harness: /tmp/split/{split,seed2}.cpp pattern (links dna_coder.cpp.o).
// Returns compressed bytes (self-contained; includes pg_len header).
std::vector<uint8_t> dna_encode(
    const std::string&              pg,
    const std::vector<std::string>& reads,
    const std::vector<uint32_t>&    chain_order,
    const std::vector<uint32_t>&    pg_pos,
    const std::string&              seed = "");

// Decompress bytes produced by dna_encode.
// pg_len is embedded in the first 8 bytes of data.
// seed: same seed used at encode time (seeds context shift-registers).
std::string dna_decode(const std::vector<uint8_t>& data,
                       const std::string&           seed = "");

// ── LZMA-ASCII pseudogenome backend (fast-decode mode) ───────────────────────
// Alternative to adaptive FCM: plain LZMA-9 on raw ASCII pg text.
// Used by ARCS_PG_FAST_DECODE=1. Trades ~0.23% ratio (GIAB) for 2.8× faster
// decompress (1.0s vs 2.8s) — LZMA streaming replay vs FCM adaptive re-training.
// Same 8-byte pg_len header as dna_encode for symmetric format handling.
std::vector<uint8_t> vle_encode_pg(const std::string& pg);
std::string          vle_decode_pg(const std::vector<uint8_t>& data);

// ── 2-bit + LZMA pg codec (format 0x08, fast-decode) ─────────────────────────
// Z4 encoding A=0,C=1,T=2,G=3 (complement = +2 mod 4). Packs 4 bases/byte
// MSB-first, then LZMA-9. Decode: ~10ms vs 26s FCM for 12MB pg (DS2).
// Ratio: +2.2% vs FCM on DS2 (M2 measurement, trial24). Non-ACTG pg falls
// back to vle_encode_pg (LZMA-ASCII, format 0x07-style payload).
std::vector<uint8_t> pg_encode_2bit(const std::string& pg);
std::string          pg_decode_2bit(const std::vector<uint8_t>& data);

// ── Cross-stream quality: pseudogenome "surprise" signal (idea B) ─────────────
// Computes, for each base of the pseudogenome, a quantized local-predictability
// ("surprise") bucket in [0, nbuckets): a low-order causal finite-context model
// scans the pg once and, before observing each base, records how confidently it
// predicted that base from the preceding `order` bases. Bucket 0 = highly
// predictable (e.g. homopolymer / low-complexity), higher = less predictable.
//
// It is a pure deterministic function of `pg` only, so the decoder reproduces it
// bit-identically from the decoded pseudogenome — no bytes are stored. Used as an
// extra quality-coder context (sequencer confidence tends to vary with local
// sequence composition). See rans_model.h::N_QUAL_SURPRISE_BINS.
std::vector<uint8_t> compute_pg_surprise(const std::string& pg,
                                         int order    = 4,
                                         int nbuckets = 4);
