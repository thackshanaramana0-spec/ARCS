#pragma once
#include "common.h"
#include "fastq_io.h"
#include <zlib.h>
#include <string>

// ── AKC Entropy-Class Router — Adaptive Backend Selection ────────────────────
//
// ALGORITHM:
//   Sample --pilot reads (default 5,000) stratified from the FASTQ.
//   Compress their concatenated sequences with zlib-9.
//   AKC score = compressed_size / raw_size.
//
// ROUTING (binary decision on score):
//   score < 0.15 (amplicon/targeted) → BSC/LZMA-9 amplicon path (use_bsc = true)
//   score >= 0.15 (WGS/metagenomic)  → LSH k-NN MST path (use_bsc = false)
//
// DIAGNOSTIC REGIMES (for `arcs akc` output and logging — not used for routing):
//   AMPLICON    score < 0.05  : ultra-high redundancy (viral amplicon, 16S)
//   TARGETED    0.05–0.15     : targeted panels, exome
//   WGS         0.15–0.35     : shotgun WGS
//   METAGENOMIC score >= 0.35 : multi-organism
//
// THEORY (MDL principle, Grunwald 2007): pilot compressed size is a proxy for
//   the description length of the full dataset under each codec.
//   Selecting min(size_bsc, size_lzma) on the pilot is equivalent to LOO-CV
//   model selection and is consistent (converges to the true optimal codec).
//   Same principle: Zstandard --adapt (RFC 8878), ClickHouse adaptive codec.

class AKCRouter {
public:
    AKCRouter() = default;

    // Compute AKC from pre-loaded pilot reads.
    AKCResult compute(const std::vector<Read>& pilot_reads) const;

    // Compute AKC from FASTQ file (reads stratified pilot internally).
    AKCResult compute_from_file(const std::string& fastq_path,
                                size_t n_pilot = 1000) const;

    // Classify score to regime
    static DataRegime classify(float score,
                               float thr_amplicon = 0.15f,
                               float thr_wgs      = 0.35f);

private:
    // Compress a byte buffer with zlib level 9. Returns compressed size.
    static size_t zlib9_compress(const std::vector<uint8_t>& data);
};
