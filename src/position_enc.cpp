#include "position_enc.h"
#include <algorithm>
#include <cmath>
#include <numeric>

// ── SE positions ──────────────────────────────────────────────────────────────
std::vector<uint8_t> encode_se_positions(const std::vector<PositionEntry>& entries) {
    std::vector<uint8_t> out;
    out.reserve(entries.size() * 3); // estimate

    write_varint(out, entries.size());

    pos_t prev = 0;
    for (const auto& e : entries) {
        // Store read_idx as varint
        write_varint(out, e.read_idx);
        // Store position delta as zigzag varint
        int32_t delta = (int32_t)e.ref_pos - (int32_t)prev;
        write_varint(out, zigzag_encode(delta));
        prev = e.ref_pos;
    }
    return out;
}

std::vector<PositionEntry> decode_se_positions(const uint8_t* data, size_t len,
                                                size_t n_entries) {
    const uint8_t* p   = data;
    const uint8_t* end = data + len;

    size_t n = (size_t)read_varint(p, end);
    (void)n_entries; // n_entries is a hint; actual count from stream

    std::vector<PositionEntry> result;
    result.reserve(n);

    pos_t prev = 0;
    for (size_t i = 0; i < n && p < end; ++i) {
        uint32_t idx = (uint32_t)read_varint(p, end);
        int32_t  delta = zigzag_decode((uint32_t)read_varint(p, end));
        pos_t    pos   = (pos_t)((int32_t)prev + delta);
        result.push_back({idx, pos});
        prev = pos;
    }
    return result;
}

// ── PE insert size estimation ──────────────────────────────────────────────────
PEPositionStats estimate_insert_stats(const std::vector<pos_t>& r1_pos,
                                       const std::vector<pos_t>& r2_pos,
                                       size_t max_samples) {
    PEPositionStats stats;
    size_t n = std::min({r1_pos.size(), r2_pos.size(), max_samples});
    if (n == 0) return stats;

    // Compute insert sizes: R2_pos - R1_pos
    // Only use pairs where R2_pos > R1_pos (expected orientation)
    std::vector<int32_t> inserts;
    inserts.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        int32_t ins = (int32_t)r2_pos[i] - (int32_t)r1_pos[i];
        if (ins > 0 && ins < 2000) // sanity bounds
            inserts.push_back(ins);
    }

    if (inserts.empty()) return stats;

    double mean = 0.0;
    for (int32_t x : inserts) mean += x;
    mean /= inserts.size();

    double var = 0.0;
    for (int32_t x : inserts) var += (x - mean) * (x - mean);
    var /= inserts.size();

    stats.insert_mean = (int32_t)std::round(mean);
    stats.insert_std  = (int32_t)std::ceil(std::sqrt(var));
    stats.n_samples   = (uint32_t)inserts.size();
    return stats;
}

// ── PE positions ──────────────────────────────────────────────────────────────
std::vector<uint8_t> encode_pe_positions(const std::vector<PEPositionEntry>& entries,
                                          const PEPositionStats& stats) {
    std::vector<uint8_t> out;
    out.reserve(entries.size() * 5);

    // Write stats header
    auto write32s = [&](int32_t v) {
        write_varint(out, zigzag_encode(v));
    };
    write32s(stats.insert_mean);
    write32s(stats.insert_std);
    write_varint(out, entries.size());

    pos_t prev_r1 = 0;
    for (const auto& e : entries) {
        // R1: delta-encoded sorted position + rc flag
        int32_t r1_delta = (int32_t)e.r1_pos - (int32_t)prev_r1;
        write_varint(out, zigzag_encode(r1_delta));
        out.push_back(e.r1_rc ? 1 : 0);

        // R2: delta from (R1 + insert_mean) + rc flag
        // This exploits insert size constraint: R2_delta is typically ±3σ ≈ ±90bp
        // Encoding ±90 range costs ~8 bits vs 25 bits for absolute position
        int32_t r2_expected = (int32_t)e.r1_pos + stats.insert_mean;
        int32_t r2_delta    = (int32_t)e.r2_pos - r2_expected;
        write_varint(out, zigzag_encode(r2_delta));
        out.push_back(e.r2_rc ? 1 : 0);

        // Pair index
        write_varint(out, e.pair_idx);
        prev_r1 = e.r1_pos;
    }
    return out;
}

std::vector<PEPositionEntry> decode_pe_positions(const uint8_t* data, size_t len,
                                                   size_t /*n_entries*/,
                                                   const PEPositionStats& /*hint*/) {
    const uint8_t* p   = data;
    const uint8_t* end = data + len;

    PEPositionStats stats;
    stats.insert_mean = zigzag_decode((uint32_t)read_varint(p, end));
    stats.insert_std  = zigzag_decode((uint32_t)read_varint(p, end));
    size_t n          = (size_t)read_varint(p, end);

    std::vector<PEPositionEntry> result;
    result.reserve(n);

    pos_t prev_r1 = 0;
    for (size_t i = 0; i < n && p < end; ++i) {
        int32_t r1_delta = zigzag_decode((uint32_t)read_varint(p, end));
        bool    r1_rc    = (*p++ != 0);
        pos_t   r1_pos   = (pos_t)((int32_t)prev_r1 + r1_delta);

        int32_t r2_delta   = zigzag_decode((uint32_t)read_varint(p, end));
        bool    r2_rc      = (*p++ != 0);
        pos_t   r2_pos     = (pos_t)((int32_t)r1_pos + stats.insert_mean + r2_delta);

        uint32_t pair_idx  = (uint32_t)read_varint(p, end);
        result.push_back({pair_idx, r1_pos, r2_pos, r1_rc, r2_rc});
        prev_r1 = r1_pos;
    }
    return result;
}
