#pragma once
#include "common.h"
#include "../third_party/rans_byte.h"
#include <vector>
#include <cstdint>

// ── Zigzag + varint encoding ──────────────────────────────────────────────────
// Converts signed integers to unsigned for varint encoding.
// Zigzag: n → 2n if n>=0, -2n-1 if n<0 (interleaves pos/neg values)
inline uint32_t zigzag_encode(int32_t n) {
    return (uint32_t)((n << 1) ^ (n >> 31));
}
inline int32_t zigzag_decode(uint32_t n) {
    return (int32_t)((n >> 1) ^ -(int32_t)(n & 1));
}

// Variable-length integer encoding (unsigned, little-endian groups of 7 bits)
inline void write_varint(std::vector<uint8_t>& out, uint64_t v) {
    while (v >= 0x80) {
        out.push_back((uint8_t)(v | 0x80));
        v >>= 7;
    }
    out.push_back((uint8_t)v);
}

inline uint64_t read_varint(const uint8_t*& p, const uint8_t* end) {
    uint64_t v = 0; int shift = 0;
    while (p < end) {
        uint8_t b = *p++;
        v |= (uint64_t)(b & 0x7F) << shift;
        if (!(b & 0x80)) break;
        shift += 7;
    }
    return v;
}

// ── SE position encoder ───────────────────────────────────────────────────────
// Encodes sorted mapped positions as zigzag-varint deltas.
// Input: list of (read_index, ref_position) pairs sorted by position.
struct PositionEntry {
    uint32_t read_idx;
    pos_t    ref_pos;
};

std::vector<uint8_t> encode_se_positions(const std::vector<PositionEntry>& entries);
std::vector<PositionEntry> decode_se_positions(const uint8_t* data, size_t len,
                                                size_t n_entries);

// ── PE position encoder ────────────────────────────────────────────────────────
// Encodes paired-end positions using:
//   - R1: sorted, delta-encoded
//   - R2: delta from R1 + insert_mean (typically ±30bp range → ~6 bits vs 25 bits)
//
// Insert size statistics estimated from first 1000 mapped pairs.
struct PEPositionStats {
    int32_t  insert_mean = 350; // mean R2_pos - R1_pos
    int32_t  insert_std  = 30;  // std of above
    uint32_t n_samples   = 0;
};

PEPositionStats estimate_insert_stats(
    const std::vector<pos_t>& r1_positions,
    const std::vector<pos_t>& r2_positions,
    size_t max_samples = 1000
);

struct PEPositionEntry {
    uint32_t pair_idx;
    pos_t    r1_pos;
    pos_t    r2_pos;
    bool     r1_rc;
    bool     r2_rc;
};

std::vector<uint8_t> encode_pe_positions(const std::vector<PEPositionEntry>& entries,
                                          const PEPositionStats& stats);
std::vector<PEPositionEntry> decode_pe_positions(const uint8_t* data, size_t len,
                                                   size_t n_entries,
                                                   const PEPositionStats& stats);
