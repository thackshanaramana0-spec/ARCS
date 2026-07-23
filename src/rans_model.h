#pragma once
#include "common.h"
#include "../third_party/rans_byte.h"
#include <vector>
#include <unordered_map>
#include <functional>

// ── Context key ───────────────────────────────────────────────────────────────
// Encodes the context for a context-dependent probability model as a uint32_t.
// Layout: bits[15:8] = context_id, bits[7:0] = not used (for alignment).
using ContextKey = uint32_t;

// ── Context-dependent frequency model ────────────────────────────────────────
// Maps ContextKey → FreqTable (one per context).
// Used for:
//   - Mismatch codec: P(base | pos_bin[10], gc_bin[4])
//   - Quality codec:  P(q | q_prev[43], is_dev[2], pos_bin[20])

class ContextModel {
public:
    ContextModel(int alphabet_size, int n_pos_bins, int extra_dims = 1)
        : alpha_(alphabet_size), n_pos_(n_pos_bins), extra_(extra_dims) {}

    // Add observation: symbol sym in context key ctx.
    void observe(ContextKey ctx, uint8_t sym);

    // Finalize: build FreqTables from accumulated counts.
    void finalize();

    // Encode one symbol. Returns encoded bytes (appended to out).
    void encode_sym(std::vector<uint8_t>& out,
                    RansEncoder& enc,
                    ContextKey ctx,
                    uint8_t sym) const;

    // Decode one symbol given current decoder state.
    uint8_t decode_sym(RansDecoder& dec,
                       const uint8_t*& ptr,
                       const uint8_t* end,
                       ContextKey ctx) const;

    // Encode/decode a full sequence using the context model.
    // ctx_fn(i, sym_so_far) returns the ContextKey for position i.
    std::vector<uint8_t> encode_sequence(
        const std::vector<uint8_t>& symbols,
        const std::function<ContextKey(size_t, const std::vector<uint8_t>&)>& ctx_fn
    ) const;

    std::vector<uint8_t> decode_sequence(
        const uint8_t* data, size_t data_len, size_t n_symbols,
        const std::function<ContextKey(size_t, const std::vector<uint8_t>&)>& ctx_fn
    ) const;

    // Full freq-table serialization (current format, ~90 bytes/context raw).
    std::vector<uint8_t> serialize() const;
    void deserialize(const uint8_t* data, size_t len);

    // Compact TopK serialization: stores only K significant symbols per context.
    // For a 43-symbol quality model, K≈7 vs 43 full entries → ~8x smaller.
    // Format: [n_ctx:4][alpha:2][per ctx: key:4, K:1, K×{sym:1,freq:2}]
    //   Remaining (alpha-K) symbols share the residual probability uniformly.
    // sig_threshold: minimum freq to be considered "significant" (default 100/16384 ≈ 0.6%).
    std::vector<uint8_t> serialize_topk(uint16_t sig_threshold = 100) const;
    void deserialize_topk(const uint8_t* data, size_t len);

    int alphabet_size() const { return alpha_; }

private:
    int alpha_;
    int n_pos_;
    int extra_;

    // Counts per context
    std::unordered_map<ContextKey, std::vector<uint32_t>> counts_;
    // Finalized tables per context
    std::unordered_map<ContextKey, FreqTable> tables_;

    bool finalized_ = false;

    // Get or create default table for unknown context
    const FreqTable& get_table(ContextKey ctx) const;
    FreqTable default_table_; // uniform distribution fallback
};

// ── Mismatch context key builder ──────────────────────────────────────────────
// Context for mismatch codec: (pos_bin × gc_bin)
// pos_bin: position in read / (read_len / N_POS_BINS)
// gc_bin:  GC content of local window / 25% bins

static constexpr int N_POS_BINS_MM = 10;
static constexpr int N_GC_BINS     = 4;

inline ContextKey mismatch_ctx(int pos, int read_len, float gc_frac) {
    int pos_bin = std::min(N_POS_BINS_MM - 1,
                           (int)(pos * N_POS_BINS_MM / read_len));
    int gc_bin  = std::min(N_GC_BINS - 1, (int)(gc_frac * N_GC_BINS));
    return (ContextKey)(pos_bin * N_GC_BINS + gc_bin);
}

// ── Quality context key builder ───────────────────────────────────────────────
// Context for quality codec: (q_prev × is_dev × pos_bin)
// q_prev:  previous quality value (0-42, binned to 5 levels)
// is_dev:  is current position a mismatch (0 or 1)
// pos_bin: position bin (0-19)

static constexpr int N_POS_BINS_Q = 20;
static constexpr int N_QPREV_BINS = 5;  // bin quality 0-42 into 5 levels

// ── Cross-stream quality: sequence-model "surprise" context (idea B) ──────────
// Number of buckets for the pseudogenome local-predictability (surprise) signal
// used as an extra quality context dimension in the chain-pg path. The DNA
// model's confidence at each pg base is quantized into this many levels and
// mixed into the quality context. Both encoder and decoder compute it
// identically from the (identical) decoded pseudogenome, so it is lossless and
// costs zero stored bytes. See dna_coder.h::compute_pg_surprise.
static constexpr int N_QUAL_SURPRISE_BINS = 4;

inline int qprev_bin(uint8_t q) {
    // Illumina Phred quality levels:
    // 0-9 → 0 (very low), 10-19 → 1, 20-29 → 2, 30-36 → 3, 37-42 → 4
    if (q < 10) return 0;
    if (q < 20) return 1;
    if (q < 30) return 2;
    if (q < 37) return 3;
    return 4;
}

inline ContextKey quality_ctx(uint8_t q_prev, bool is_dev, int pos, int read_len) {
    int qb  = qprev_bin(q_prev);
    int pb  = std::min(N_POS_BINS_Q - 1, (int)(pos * N_POS_BINS_Q / read_len));
    int dev = is_dev ? 1 : 0;
    return (ContextKey)(qb * N_POS_BINS_Q * 2 + dev * N_POS_BINS_Q + pb);
}
