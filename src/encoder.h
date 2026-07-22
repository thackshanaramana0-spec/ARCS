#pragma once
#include "common.h"
#include "fastq_io.h"
#include "akc.h"
#include "debruijn.h"
#include "minimizer.h"
#include "mapper.h"
#include "hamming_mst.h"
#include "mst_encoder.h"
#include "rans_model.h"
#include "position_enc.h"
#include "container.h"
#include <string>
#include <functional>

// ── Progress callback ─────────────────────────────────────────────────────────
struct EncodeProgress {
    size_t total_reads    = 0;
    size_t mapped_reads   = 0;
    size_t unmapped_reads = 0;
    size_t genome_len     = 0;
    int    k_used         = 0;
    float  akc_score      = 0.0f;
    DataRegime regime     = DataRegime::WGS;
    double elapsed_sec    = 0.0;
    std::string pg_method; // "debruijn" or "greedy_scs" or "amplicon"
};

using ProgressCallback = std::function<void(const EncodeProgress&)>;

// ── ARCSEncoder ────────────────────────────────────────────────────────────────
// Full lossless FASTQ compression pipeline.
// Single-end mode: compress(se_path, output_path)
// Paired-end mode: compress_pe(r1_path, r2_path, output_path)

class ARCSEncoder {
public:
    explicit ARCSEncoder(const ARCSParams& params = ARCSParams{});

    // Set progress callback (optional)
    void set_progress_callback(ProgressCallback cb) { progress_cb_ = cb; }

    // Fused compress-and-call: when set (chain-pg mode), the assembler's read→
    // consensus placements are captured here so a single assembly pass serves BOTH
    // the archive and reference-free variant calling — no second assembly, no
    // decompress round-trip. The reads are also captured (calling needs them).
    // Filled during compress(); consumed by the caller afterwards.
    struct CallData* call_capture_ = nullptr;
    std::vector<Read>* call_reads_ = nullptr;

    // Optional hook launched as a 3rd async task alongside quality+names encoding.
    // Set before compress() to overlap calling with the encoding phase (zero overhead).
    std::function<void()> post_assembly_call_hook_;

    // Compress single-end FASTQ
    EncodeProgress compress(const std::string& input_path,
                            const std::string& output_path);

    // Compress paired-end FASTQ
    EncodeProgress compress_pe(const std::string& r1_path,
                               const std::string& r2_path,
                               const std::string& output_path);

private:
    ARCSParams params_;
    ProgressCallback progress_cb_;
    bool crlf_ = false;      // input used CRLF line endings → reproduce on decode
    bool final_nl_ = true;   // input ended with a trailing newline

    // ── Chunked mode: compress input in independent MST chunks ───────────────
    // Each chunk of chunk_size reads is compressed as a self-contained inner
    // ARCS archive.  All inner archives are packed into one CHUNK_DATA blob
    // inside an outer (wrapper) archive.  Peak memory = O(chunk_size * L).
    void compress_chunked(const std::string& input_path,
                          const std::string& output_path,
                          EncodeProgress& prog);

    // ── WGS path — SCS pseudogenome (current) ────────────────────────────────
    void encode_wgs_se(const std::vector<Read>& reads,
                       ARCSWriter& writer,
                       EncodeProgress& prog);

    // ── WGS path — MST tree encoding (new, beats SCS by 22-50%) ─────────────
    void encode_wgs_mst(const std::vector<Read>& reads,
                        ARCSWriter& writer,
                        EncodeProgress& prog);

    void encode_wgs_chain(const std::vector<Read>& reads,
                          ARCSWriter& writer,
                          EncodeProgress& prog);

    void encode_wgs_chain_pg(const std::vector<Read>& reads,
                              ARCSWriter& writer,
                              EncodeProgress& prog);

    void encode_wgs_pe(const std::vector<Read>& r1_reads,
                       const std::vector<Read>& r2_reads,
                       ARCSWriter& writer,
                       EncodeProgress& prog);

    // ── Amplicon path (BSC) ───────────────────────────────────────────────────
    void encode_amplicon(const std::vector<Read>& reads,
                         ARCSWriter& writer,
                         EncodeProgress& prog);

    // ── Shared helpers ────────────────────────────────────────────────────────
    // Encode mismatches with context rANS. Returns bytes.
    std::vector<uint8_t> encode_mismatches(
        const std::vector<Read>& reads,
        const std::vector<MapResult>& mappings,
        const std::string& genome,
        ContextModel& mismatch_model
    ) const;

    // Encode quality with context rANS. Returns bytes.
    std::vector<uint8_t> encode_quality(
        const std::vector<Read>& reads,
        const std::vector<MapResult>& mappings,
        ContextModel& quality_model
    ) const;

    // Encode read names. Returns compressed bytes.
    // chain_order: if provided (chain-pg path), enables paired-end name dedup (format 0x06).
    //   chain_order[k] = original read index at SCS position k.
    // orig_reads: original reads before SCS reorder (needed for base-name extraction).
    std::vector<uint8_t> encode_names(
        const std::vector<Read>& reads,
        const std::vector<uint32_t>* chain_order = nullptr,
        const std::vector<Read>*     orig_reads  = nullptr) const;

    // Encode strand flags (1 bit per read). Returns bytes.
    std::vector<uint8_t> encode_strands(const std::vector<MapResult>& mappings) const;

    // Train and finalize context models on pilot reads.
    void train_models(const std::vector<Read>& reads,
                      const std::vector<MapResult>& mappings,
                      const std::string& genome,
                      ContextModel& mismatch_model,
                      ContextModel& quality_model) const;

    // Load all reads from file (streaming into memory).
    // For large files (>500MB), caller should use streaming encoder (TODO).
    std::vector<Read> load_reads(const std::string& path) const;

    // BSC compression via subprocess call — original FASTA format (legacy WGS path).
    std::vector<uint8_t> bsc_compress_sequences(
        const std::vector<Read>& reads) const;

    // BSC compression of raw byte payload (V7 amplicon: no FASTA overhead).
    // First byte of returned blob = codec flag: 0x01=LZMA fallback, 0x02=BSC.
    std::vector<uint8_t> bsc_compress_raw(const std::vector<uint8_t>& raw) const;

    // BSC decompression of a blob written by bsc_compress_raw.
    static std::vector<uint8_t> bsc_decompress_raw(const std::vector<uint8_t>& blob);
};
