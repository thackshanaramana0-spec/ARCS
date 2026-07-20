#pragma once
// ── Adaptive context-model quality coder (fqzcomp-style) ──────────────────────
// A fully ADAPTIVE quality-score coder: the frequency model updates online and is
// NEVER transmitted, so rich context costs zero archive bytes (the key lesson from
// the SOTA — see trial20_arcs_dna_crossstream/QUALITY_SOTA.md). Each quality value
// is coded by binary decomposition with adaptive 12-bit bit-models selected by a
// fqzcomp-class context: {q_{i-1}, max(q_{i-2},q_{i-3}), (q_{i-2}==q_{i-3}),
// running-delta, position, mismatch-flag}. Encoder and decoder update identically,
// so it is exactly lossless and portable (integer-only carryless arithmetic coder).
//
// rq[idx][j]        = raw Phred value (0..42) of read idx at position j.
// order[k]          = read index visited at traversal step k (chain/DFS order).
// dev_sets[idx][j]  = mismatch-vs-pg flag for read idx position j (orig-read frame).
// L                 = read length.
#include <vector>
#include <string>
#include <cstdint>

// seqs (optional, the GAME-CHANGER lever): per-read base sequence in the SAME
// frame/index as rq (orig-read frame). When supplied, quality is conditioned on
// the local sequence context (bases around each position) — the physical cause of
// quality difficulty (homopolymers/motifs), which reference-free per-read tools
// cannot exploit. Measured −10.65% on held-out data (see QUALITY_SOTA.md). The
// decoder passes the reconstructed sequence (available before quality), so it is
// lossless and stores nothing. When null, behaviour is the quality-history model.
std::vector<uint8_t> qual_cm_encode(
    const std::vector<std::vector<uint8_t>>& rq,
    const std::vector<uint32_t>&             order,
    const std::vector<std::vector<bool>>&    dev_sets,
    int                                      L,
    const std::vector<std::string>*          seqs = nullptr);

// Decodes a stream produced by qual_cm_encode into rq_out (indexed by read index,
// same frame as `order`). rq_out is sized [max(order)+1][L]. Returns false on a
// malformed stream. seqs must match what the encoder used (same frame as rq_out).
bool qual_cm_decode(
    const std::vector<uint8_t>&              data,
    const std::vector<uint32_t>&             order,
    const std::vector<std::vector<bool>>&    dev_sets,
    int                                      L,
    std::vector<std::vector<uint8_t>>&       rq_out,
    const std::vector<std::string>*          seqs  = nullptr,
    const std::vector<int>*                  rlens = nullptr);   // per-read lengths
