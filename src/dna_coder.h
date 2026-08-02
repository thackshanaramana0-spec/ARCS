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
// seed:         up to max-FCM-order chars from the preceding block's pg tail;
//               seeds the context shift-registers so split-block coding starts
//               in the correct context state (lossless parallel-block support).
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
