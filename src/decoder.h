#pragma once
#include "common.h"
#include "container.h"
#include "fastq_io.h"
#include "rans_model.h"
#include "position_enc.h"
#include <string>

// ── ARCSDecoder ────────────────────────────────────────────────────────────────
// Decodes an ARCS container back to FASTQ.
// Guarantees: decoded FASTQ is byte-exact with original (lossless).

class ARCSDecoder {
public:
    ARCSDecoder() = default;

    // Decompress ARCS archive to FASTQ.
    void decompress(const std::string& input_path,
                    const std::string& output_path);

    // Info: print container metadata without decompressing.
    void print_info(const std::string& input_path);

private:
    // Decode WGS path (old SCS format)
    void decode_wgs(const ARCSReader& rdr, FASTQWriter& out);

    // Decode WGS MST path (V7 format — MST_TREE + MST_DELTAS + delta-rANS quality)
    void decode_wgs_mst(const ARCSReader& rdr, FASTQWriter& out);

    // Decode WGS chain path (trial11 — CHAIN_VSEQ + MST_DELTAS + CHAIN_STARTS)
    void decode_wgs_chain(const ARCSReader& rdr, FASTQWriter& out);

    // Decode WGS chain-pg path (CHAIN_PG_SEQ + CHAIN_PG_POS + CHAIN_PG_AUX)
    void decode_wgs_chain_pg(const ARCSReader& rdr, FASTQWriter& out);

    // Decode amplicon path (BSC)
    void decode_amplicon(const ARCSReader& rdr, FASTQWriter& out);

    // Apply mismatches to reconstruct original sequence from pseudogenome substring.
    // ref: substring of pseudogenome at read's mapping position
    // mm_offsets, mm_bases: mismatch list (positions + original bases)
    // rc: whether read was reverse-complemented before mapping
    std::string reconstruct_seq(const std::string& ref_substr,
                                 const std::vector<uint32_t>& mm_offsets,
                                 const std::vector<uint8_t>&  mm_bases,
                                 bool rc) const;

    // Decode quality stream for one read using context model.
    std::string decode_quality(
        const uint8_t*& ptr,
        const uint8_t* end,
        RansDecoder& dec,
        const ContextModel& q_model,
        const std::vector<bool>& is_dev,
        int read_len) const;
};
