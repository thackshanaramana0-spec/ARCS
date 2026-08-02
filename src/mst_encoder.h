#pragma once
#include "common.h"
#include "hamming_mst.h"
#include "../third_party/rans_byte.h"
#include "rans_model.h"
#include <vector>
#include <string>
#include <cstdint>

// ── MST-based sequence encoder ─────────────────────────────────────────────────
//
// Replaces the linear SCS pseudogenome approach with a TREE-based encoding.
// Theory: SCS forces linear chains. MST finds globally optimal connections
//         (each read connects to its closest neighbor, not just its sequential
//         neighbour). Proven 22-50% better than SCS/SPRING by Hamming-Shifting
//         graph paper (PLOS Comp Bio 2021).
//
// Our advance over Mstcom (which was impractical at 112GB/11552s):
//   - Minimizer-based k-NN: O(n) not O(n²) pairwise
//   - Approximate MST on k-NN subgraph: fast, near-optimal
//   - Delta encoding with rANS: compact storage
//
// Encoding per read in DFS order:
//   Root: verbatim sequence (L bases × 2 bits)
//   Non-root: (shift, substitution_list, non-overlap_prefix, non-overlap_suffix)
//     shift     : how many bases child is shifted relative to parent (zigzag varint)
//     n_subs    : number of substitutions in overlap (varint)
//     sub_pos[] : substitution positions in overlap (delta-varint)
//     sub_base[]: substitution bases (rANS coded with P(base|pos_bin))
//     prefix    : non-overlapping prefix of child (if shift > 0) (rANS coded)
//     suffix    : non-overlapping suffix of child (if shift < 0) (rANS coded)

struct ReadAlignResult {
    int     shift;           // child is shifted by `shift` bases relative to parent
                             // shift > 0: child extends to the right of parent
                             // shift < 0: child extends to the left of parent
                             // shift = 0: same position, just substitutions
    int     overlap_len;     // length of aligned overlap region
    int     n_subs;          // Hamming distance in overlap
    std::vector<uint32_t> sub_positions; // within-overlap offsets (0-based from overlap start)
    std::vector<uint8_t>  sub_bases;     // correct bases at those positions (in child)
    bool    rc_child;        // child was reverse-complemented before aligning
};

// Align child to parent. Tries all shifts in [-max_shift, +max_shift] and
// both strand orientations of child. Returns the best alignment.
ReadAlignResult align_reads(const std::string& parent,
                             const std::string& child,
                             int max_shift = 50,
                             int max_mm    = 10);

// ── MST encoder result ─────────────────────────────────────────────────────────
struct MSTEncodeResult {
    std::vector<uint8_t> tree_bytes;        // parent indices (compressed)
    std::vector<uint8_t> delta_bytes;       // shift + sub lists (metadata only, no bases)
    std::vector<uint8_t> seq_text_bytes;    // ACGT chars: root seqs + sub bases + nonoverlap segs
    std::vector<uint8_t> quality_bytes;     // quality rANS stream
    std::vector<uint8_t> quality_model;     // serialized ContextModel for quality
    std::vector<uint8_t> rc_flags;          // RC flags (bitpacked)

    size_t n_reads           = 0;
    int    read_len          = 0;
    size_t n_root_reads      = 0;
    size_t n_delta_reads     = 0;
    double avg_subs_per_read = 0;
    double avg_overlap       = 0;

    // DFS traversal order (read index → DFS position)
    // Needed by decoder to reorder quality back to original order.
    std::vector<uint32_t> dfs_order;
};

// ── MST encoder config ─────────────────────────────────────────────────────────
struct MSTConfig {
    int  k_neighbors = 20;   // k-NN graph: k nearest neighbors per read
    int  n_lsh       = 8;    // LSH hash functions for k-NN
    int  max_shift   = 50;   // maximum shift (bases) to try when aligning
    int  max_mm      = -1;   // max substitutions accepted (-1 = read_len/10)
    int  w           = 10;   // minimizer window for k-NN seeding
    int  mink        = 15;   // minimizer k for k-NN seeding
    bool try_rc      = true; // also try reverse complement of child
    int  pilot       = 5000; // pilot reads for quality model training
};

// ── Main encoder class ─────────────────────────────────────────────────────────
class MSTSequenceEncoder {
public:
    explicit MSTSequenceEncoder(const MSTConfig& cfg = MSTConfig{});

    // Encode reads using MST tree structure.
    MSTEncodeResult encode(const std::vector<Read>& reads);

    // Encode quality using rANS with MST deviation context.
    // dev_sets[i] = is_dev[j] = true if position j of read dfs_order[i]
    //               is a substitution from its MST parent.
    // Returns (quality_bytes, model_bytes).
    // encode_quality_rans: parent-prediction residuals + LZMA-9
    // parents[i] = tree parent of read i (UINT32_MAX for roots)
    // aln_shifts[i] = shift of read i relative to parent (0 for roots)
    // surprise (optional, idea B): surprise[read_idx][pos] = pg local-predictability
    // bucket for that base (see dna_coder.h::compute_pg_surprise). When non-null it
    // is mixed into the quality context (mode 0x05/0x06); when null, behaviour and
    // output bytes are identical to before. Decoder recomputes the same buckets.
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
    encode_quality_rans(const std::vector<Read>&              reads,
                        const std::vector<uint32_t>&          dfs_order,
                        const std::vector<uint32_t>&          parents,
                        const std::vector<int>&               aln_shifts,
                        const std::vector<std::vector<bool>>& dev_sets,
                        int                                   read_len,
                        const std::vector<std::vector<uint8_t>>* surprise = nullptr) const;

    // Public accessor used by ChainSequenceEncoder (chain_encoder.cpp).
    std::vector<KNNEdge> get_knn_edges(const std::vector<std::string>& seqs) const {
        return build_knn(seqs);
    }

private:
    MSTConfig cfg_;

    // Build k-NN graph using minimizer-bucket approximate nearest neighbors.
    // For each pair in the same minimizer bucket, compute exact Hamming+shift distance.
    // Keep top-k neighbors per read.
    std::vector<KNNEdge> build_knn(const std::vector<std::string>& seqs) const;

    // Encode the tree structure (parent array) compactly.
    // Format: [n_reads:varint] [parent_0:varint] ... [parent_{n-1}:varint]
    // Root nodes have parent = UINT32_MAX.
    std::vector<uint8_t> encode_tree(const std::vector<uint32_t>& parents,
                                      size_t n_reads) const;

    // Encode all read deltas in DFS order.
    // Also outputs:
    //   dev_sets_out[read_idx][pos] = true if position is a substitution from MST parent
    //   aln_shifts_out[read_idx]    = alignment shift of read relative to parent
    std::vector<uint8_t> encode_deltas(const std::vector<Read>&         reads,
                                        const std::vector<uint32_t>&     parents,
                                        const std::vector<uint32_t>&     dfs_order,
                                        std::vector<bool>&               rc_flags_out,
                                        std::vector<std::vector<bool>>&  dev_sets_out,
                                        std::vector<int>&                aln_shifts_out,
                                        MSTEncodeResult&                 stats) const;
};

// ── MST decoder ────────────────────────────────────────────────────────────────
class MSTSequenceDecoder {
public:
    // Decode reads from MST encoding.
    // Returns reads in ORIGINAL ORDER.
    // Optionally outputs dfs_order (DFS traversal, used for quality decoding).
    // Optionally outputs dfs_order, parents, and dev_sets.
    // dev_sets[i][j] = true iff position j of read i is a substitution from its MST parent.
    // Reconstructed directly from the stored sub_positions — no extra archive data needed.
    std::vector<std::string> decode_sequences(
        const std::vector<uint8_t>& tree_bytes,
        const std::vector<uint8_t>& delta_bytes,
        const std::vector<uint8_t>& rc_flags,
        size_t                      n_reads,
        int                         read_len,
        std::vector<uint32_t>*           dfs_order_out = nullptr,
        std::vector<uint32_t>*           parents_out   = nullptr,
        std::vector<std::vector<bool>>*  dev_sets_out  = nullptr,
        const std::vector<uint8_t>*      seq_text      = nullptr
    ) const;

    // Decode quality stream back to per-read quality strings in ORIGINAL ORDER.
    // dev_sets: per-read deviation flags from decode_sequences() (used for is_dev context).
    // mode byte is first byte of quality_model_bytes:
    //   0x00 = LZMA fallback (DFS-order concatenated raw Phred)
    //   0x01 = context-rANS (q_prev_exact × is_dev × pos_bin — 860 contexts)
    //   0x03 = parallel blocked context-rANS
    //   0x05 = mode 0x01 + surprise context (idea B); requires surprise non-null
    //   0x06 = mode 0x03 + surprise context (idea B); requires surprise non-null
    // surprise (optional): surprise[read_idx][pos] recomputed by the caller from
    // the decoded pseudogenome — must match what the encoder used.
    std::vector<std::string> decode_quality(
        const std::vector<uint8_t>&        quality_bytes,
        const std::vector<uint8_t>&        quality_model_bytes,
        const std::vector<uint32_t>&       dfs_order,
        const std::vector<uint32_t>&       parents,
        const std::vector<std::vector<bool>>& dev_sets,
        size_t                             n_reads,
        int                                read_len,
        const std::vector<std::vector<uint8_t>>* surprise = nullptr
    ) const;
};
