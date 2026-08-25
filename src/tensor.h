#pragma once
#include <string>

// ── GU-256 Tensor Exporter ────────────────────────────────────────────────────
// Reads an ARCS archive directly (no FASTQ needed) and emits ML-ready feature
// tensors from the placement records already stored for lossless decompression.
//
// Outputs:
//   {prefix}_sample.json  - sample-level feature vector (always)
//   {prefix}_reads.tsv    - per-read feature matrix     (if --reads flag)
//   {prefix}_coverage.bin - uint32 or binned coverage   (if --coverage flag)
//
// All features are reference-free: no reference genome required at any stage.

struct TensorOptions {
    bool emit_reads    = true;   // write per-read TSV
    bool emit_coverage = false;  // write coverage array (large for WGS)
    size_t coverage_bins = 50000; // max bins for binned coverage (auto-bins if pg > 50MB)
};

int run_tensor_export(const std::string& archive_path,
                      const std::string& output_prefix,
                      const TensorOptions& opts = TensorOptions{});
