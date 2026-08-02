#include "akc.h"
#include <zlib.h>
#include <numeric>
#include <cstring>

// ── zlib-9 compression ────────────────────────────────────────────────────────
size_t AKCRouter::zlib9_compress(const std::vector<uint8_t>& data) {
    if (data.empty()) return 0;

    uLongf bound = compressBound((uLong)data.size());
    std::vector<uint8_t> buf(bound);
    uLongf out_len = bound;

    int ret = compress2(buf.data(), &out_len,
                        data.data(), (uLong)data.size(),
                        Z_BEST_COMPRESSION); // level 9
    if (ret != Z_OK) return data.size(); // fallback: no compression
    return (size_t)out_len;
}

// ── AKCRouter::classify ───────────────────────────────────────────────────────
// Returns diagnostic regime. ROUTING is binary: score < thr_amplicon (0.15) → BSC path.
DataRegime AKCRouter::classify(float score, float thr_amplicon, float thr_wgs) {
    // thr_amplicon = 0.15 (routing boundary: < 0.15 -> BSC, >= 0.15 -> PgRC2)
    // Internal sub-regimes within BSC cluster (for logging only):
    if (score < 0.05f)         return DataRegime::AMPLICON;   // ultra-high redundancy
    if (score < thr_amplicon)  return DataRegime::TARGETED;   // targeted panels
    if (score < thr_wgs)       return DataRegime::WGS;        // shotgun WGS
    return DataRegime::METAGENOMIC;                            // multi-organism
}

// ── AKCRouter::compute ────────────────────────────────────────────────────────
AKCResult AKCRouter::compute(const std::vector<Read>& pilot_reads) const {
    if (pilot_reads.empty())
        return {0.5f, DataRegime::WGS, false};

    // Concatenate sequences only (AKC uses sequence, not quality)
    std::vector<uint8_t> concat;
    concat.reserve(pilot_reads.size() * pilot_reads[0].seq.size());
    for (const auto& r : pilot_reads)
        for (char c : r.seq) concat.push_back((uint8_t)c);

    size_t raw_size  = concat.size();
    size_t comp_size = zlib9_compress(concat);

    float score = (raw_size > 0) ? (float)comp_size / raw_size : 1.0f;
    DataRegime regime = classify(score);
    bool use_bsc = (regime == DataRegime::AMPLICON ||
                    regime == DataRegime::TARGETED);

    return {score, regime, use_bsc};
}

// ── AKCRouter::compute_from_file ─────────────────────────────────────────────
AKCResult AKCRouter::compute_from_file(const std::string& fastq_path,
                                        size_t n_pilot) const {
    auto pilot = sample_stratified(fastq_path, n_pilot, 3);
    return compute(pilot);
}
