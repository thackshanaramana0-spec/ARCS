#pragma once
#include "common.h"
#include "mst_encoder.h"
#include <vector>
#include <string>
#include <cstdint>

// ── Greedy chain encoder (trial11) ────────────────────────────────────────────
//
// Fixes trial10's two fatal bugs:
//   Bug 1: min-k-mer sort → ~1.5% delta rate on DS-1. Fix: k-NN graph + greedy chains.
//   Bug 2: CHAIN_READ_PERM cost 200+ KB. Fix: output in chain order, no permutation.
//
// Algorithm:
//   1. Build k-NN overlap graph via MSTSequenceEncoder::get_knn_edges() (reuses all
//      LSH/minimizer infrastructure from the MST encoder — same quality graph).
//   2. Greedy chain extraction: process edges ascending by distance; use Union-Find
//      to prevent cycles. Each read gets at most 1 predecessor and 1 successor.
//   3. Walk chains (chain-starts first, follow successors); delta-encode each read
//      against its predecessor using align_reads() (proper shift + substitutions).
//   4. Store CHAIN_STARTS bitset only (12.5 KB) — NO CHAIN_READ_PERM.
//
// Overhead now matches SPRING: chain-start bitset = 12.5 KB fixed cost.
// Delta stream quality >= SPRING because k-NN graph finds globally better overlaps.

struct ChainEncodeResult {
    // Verbatim sequences + N-positions (separate LZMA blob for better compression)
    std::vector<uint8_t>  vseq_bytes;     // col_vseq + col_vN (chain-start reads only)
    // Delta-field columns: RC flags, shifts, n_subs, sub_pos, sub_bases, extras, N-pos
    std::vector<uint8_t>  delta_bytes;    // all delta fields (NO verbatim sequences)
    std::vector<uint8_t>  chain_starts;   // packed bitset: bit i=1 if chain pos i is a chain-start
    // No sorted_order field — reads are stored and output in chain order, no permutation.
    std::vector<uint8_t>  rc_flags;       // packed bitset: bit i=1 if read at chain pos i used RC
    std::vector<uint8_t>  quality_bytes;  // rANS quality stream
    std::vector<uint8_t>  quality_model;  // serialized ContextModel
    size_t n_reads        = 0;
    size_t n_chain_starts = 0;
    size_t n_deltas       = 0;
    int    read_len       = 0;

    // ── chain_order: original read index at each chain position (chain-pg mode) ──
    std::vector<uint32_t> chain_order;   // chain_order[i] = original read idx at chain pos i

    // ── Pseudogenome output (populated when ChainSequenceEncoder build_pg=true) ──
    bool          has_pg = false;
    std::string   pg;                    // pseudogenome (ACGT only, N → 'A')
    // Per-chain-position flat streams (length = n_reads, chain traversal order):
    std::vector<uint32_t> pg_pos;        // pg start position for read at chain pos i
    std::vector<uint8_t>  pg_rc;         // RC flag per chain pos (0 or 1)
    std::vector<uint16_t> pg_readlen;    // per-read length (supports variable-length reads)
    // Sequence mismatch data (pg frame): for each read, n_mm mismatches vs pg
    std::vector<uint32_t> pg_mm_counts;  // mm count per chain-order read
    std::vector<uint16_t> pg_mm_pos_flat;// mismatch positions in enc_seq frame (0..L-1)
    std::vector<uint8_t>  pg_mm_base_flat;// mismatch bases (2-bit code: A=0,C=1,G=2,T=3)
    // N-position data (enc_seq frame): non-ACGT bases masked to 'A' in pg. The
    // ORIGINAL byte (N, lowercase acgt/n, IUPAC R/Y/…, any non-ACGT) is stored in
    // pg_N_char_flat so decode restores it EXACTLY (true byte-losslessness).
    std::vector<uint32_t> pg_N_counts;   // N count per chain-order read
    std::vector<uint16_t> pg_N_pos_flat; // N positions in enc_seq frame (0..L-1)
    std::vector<uint8_t>  pg_N_char_flat;// literal original byte at each N position
    // Quality deviation data (orig-read frame): needed for quality context in decoder
    // For shift>=0 deltas: qmm = mm converted to orig frame.
    // For shift<0 deltas: qmm = chain alignment sub_positions converted to orig frame.
    // For chain_starts: qmm is empty.
    std::vector<uint32_t> pg_qmm_counts;  // quality-mm count per chain-order read
    std::vector<uint16_t> pg_qmm_pos_flat;// quality-mm positions in orig-read frame (0..L-1)
};

class ChainSequenceEncoder {
public:
    explicit ChainSequenceEncoder(const MSTConfig& cfg = MSTConfig{}, bool build_pg = false)
        : cfg_(cfg), build_pg_(build_pg) {}

    ChainEncodeResult encode(const std::vector<Read>& reads);

private:
    MSTConfig cfg_;
    bool      build_pg_ = false;
};

// ── Dedup pseudogenome assembler (chain-pg alternative) ──────────────────────
// Builds a compact pseudogenome by greedy reference assembly with global
// approximate remapping: each read is mapped (position + orientation +
// mismatches) onto an already-built region of the pg when a low-mismatch match
// exists (a re-observation of the same locus), otherwise appended as new
// reference. This collapses the ~30x coverage redundancy the greedy chain
// leaves behind, shrinking the pg toward genome size.
//
// Produces the same ChainEncodeResult pg_* fields the chain-pg encoder/decoder
// consume, so the existing serialization, quality path, and decoder work
// unchanged. Losslessness holds by construction: applying the recorded
// mismatches + N restores to pg[pos..pos+L] yields RC^rc(orig_read), and RC is
// involutive so the decoder recovers orig_read exactly.
ChainEncodeResult build_dedup_pg(const std::vector<Read>& reads);

// ── Multi-contig greedy assembler (chain-pg alternative) ─────────────────────
// Grows many independent contigs at once: each read extends whichever contig's
// end it overlaps (via a global k-mer index), maps internally when fully
// contained, or starts a new contig. Because contigs grow independently there
// is no O(n^2) shifting, and — unlike the single-pass frontier assembler — a
// read may extend ANY contig's end, not just the one global frontier. This
// captures interior-overlap reads that the single-pass version fragments into
// appends, shrinking the pg toward genome size.
//
// Own design (greedy multi-contig growth), distinct from the overlap-layout /
// approximate-mapping pipeline in the reference; produces the same
// ChainEncodeResult pg_* fields, so serialization/quality/decoder are unchanged.
// Lossless by the same RC-involutive argument as build_dedup_pg.
// Reference-free variant calling needs the read→consensus placement the assembler
// already computes (contig id + position + strand per read). When a CallData* is
// passed, build_multicontig_pg fills it with the final (post-merge) contigs and
// per-original-read placement — the same data ARCS_PILEUP_SAM emits — so the
// caller can build its pileup with no external aligner. Off (nullptr) by default.
struct CallData {
    std::vector<std::string> contigs;   // consensus contigs (post-merge)
    std::vector<uint32_t>    read_cid;  // [orig read idx] → contig id
    std::vector<uint32_t>    read_pos;  // [orig read idx] → start position in contig (contig frame)
    std::vector<uint8_t>     read_rc;   // [orig read idx] → reverse-complement flag
    bool valid = false;
};
ChainEncodeResult build_multicontig_pg(const std::vector<Read>& reads,
                                       CallData* call_out = nullptr);

// ── Chain decoder ─────────────────────────────────────────────────────────────
class ChainSequenceDecoder {
public:
    // Decode sequences from chain encoding.
    // sorted_order[i] = original read index at chain position i.
    // In trial11 this is always identity {0,1,...,n-1} (chain order = output order).
    std::vector<std::string> decode_sequences(
        const std::vector<uint8_t>&  vseq_bytes,        // verbatim seqs + N-pos
        const std::vector<uint8_t>&  delta_bytes,        // delta fields only
        const std::vector<uint8_t>&  chain_starts_packed,
        const std::vector<uint32_t>& sorted_order,
        const std::vector<uint8_t>&  rc_flags_packed,
        size_t n_reads,
        int    read_len,
        std::vector<uint32_t>*           parents_out  = nullptr,
        std::vector<std::vector<bool>>*  dev_sets_out = nullptr
    ) const;
};
