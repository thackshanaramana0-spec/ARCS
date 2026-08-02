#include "rans_model.h"
#include <cstring>
#include <stdexcept>
#include <functional>

// ── ContextModel ──────────────────────────────────────────────────────────────
void ContextModel::observe(ContextKey ctx, uint8_t sym) {
    ARCS_CHECK(!finalized_, "Cannot observe after finalize()");
    ARCS_CHECK(sym < (uint8_t)alpha_, "Symbol out of alphabet range");
    auto& cnt = counts_[ctx];
    if (cnt.empty()) cnt.assign(alpha_, 0);
    ++cnt[sym];
}

void ContextModel::finalize() {
    ARCS_CHECK(!finalized_, "Already finalized");

    // Build default table (uniform distribution)
    std::vector<uint32_t> uniform(alpha_, 1);
    default_table_.build(uniform.data(), alpha_);

    // Build per-context tables
    for (auto& [ctx, cnts] : counts_) {
        tables_[ctx].build(cnts.data(), alpha_);
    }
    finalized_ = true;
}

const FreqTable& ContextModel::get_table(ContextKey ctx) const {
    auto it = tables_.find(ctx);
    return (it != tables_.end()) ? it->second : default_table_;
}

void ContextModel::encode_sym(std::vector<uint8_t>& out,
                               RansEncoder& enc,
                               ContextKey ctx,
                               uint8_t sym) const {
    ARCS_CHECK(finalized_, "Must finalize() before encoding");
    const FreqTable& ft = get_table(ctx);
    ARCS_CHECK(sym < (uint8_t)ft.syms.size(), "Symbol out of range");
    enc.encode(out, ft.syms[sym]);
}

uint8_t ContextModel::decode_sym(RansDecoder& dec,
                                  const uint8_t*& ptr,
                                  const uint8_t* end,
                                  ContextKey ctx) const {
    ARCS_CHECK(finalized_, "Must finalize() before decoding");
    const FreqTable& ft = get_table(ctx);
    uint32_t cf  = dec.get_cumfreq();
    uint8_t  sym = ft.decode_sym(cf);
    dec.advance(ptr, end, ft.syms[sym]);
    return sym;
}

std::vector<uint8_t> ContextModel::encode_sequence(
    const std::vector<uint8_t>& symbols,
    const std::function<ContextKey(size_t, const std::vector<uint8_t>&)>& ctx_fn
) const {
    ARCS_CHECK(finalized_, "Must finalize() before encoding");

    // rANS encodes in reverse order
    RansEncoder enc;
    std::vector<uint8_t> out;
    out.reserve(symbols.size() + 16);

    // We need to encode forward but rANS outputs in reverse.
    // Solution: pre-compute context for each position, then encode in reverse.
    std::vector<ContextKey> contexts(symbols.size());
    std::vector<uint8_t> decoded_so_far;
    decoded_so_far.reserve(symbols.size());
    for (size_t i = 0; i < symbols.size(); ++i) {
        contexts[i] = ctx_fn(i, decoded_so_far);
        decoded_so_far.push_back(symbols[i]);
    }

    // Encode in reverse
    for (int i = (int)symbols.size() - 1; i >= 0; --i)
        enc.encode(out, get_table(contexts[i]).syms[symbols[i]]);
    enc.flush(out);
    std::reverse(out.begin(), out.end());
    return out;
}

std::vector<uint8_t> ContextModel::decode_sequence(
    const uint8_t* data, size_t data_len, size_t n_symbols,
    const std::function<ContextKey(size_t, const std::vector<uint8_t>&)>& ctx_fn
) const {
    ARCS_CHECK(finalized_, "Must finalize() before decoding");
    if (data_len < 4) return {};

    RansDecoder dec;
    const uint8_t* ptr = data;
    const uint8_t* end = data + data_len;
    dec.init(ptr);
    ptr += 4;

    std::vector<uint8_t> out;
    out.reserve(n_symbols);

    for (size_t i = 0; i < n_symbols; ++i) {
        ContextKey ctx = ctx_fn(i, out);
        uint8_t sym    = decode_sym(dec, ptr, end, ctx);
        out.push_back(sym);
    }
    return out;
}

// ── Serialization ─────────────────────────────────────────────────────────────
// Format: [n_contexts:4][alpha:2][for each context: ctx_key:4, freqs: alpha*2 bytes]
std::vector<uint8_t> ContextModel::serialize() const {
    ARCS_CHECK(finalized_, "Must finalize() before serialize()");

    std::vector<uint8_t> out;
    uint32_t n = (uint32_t)tables_.size();

    auto write32 = [&](uint32_t v) {
        out.push_back(v >> 24); out.push_back(v >> 16);
        out.push_back(v >>  8); out.push_back(v);
    };
    auto write16 = [&](uint16_t v) {
        out.push_back(v >> 8); out.push_back(v);
    };

    write32(n);
    write16((uint16_t)alpha_);

    // Sort by context key so adjacent contexts (same q_prev, similar distributions)
    // are written together → LZMA finds cross-context patterns → ~5-8x better compression.
    std::vector<std::pair<ContextKey, const FreqTable*>> sorted;
    sorted.reserve(tables_.size());
    for (const auto& [ctx, ft] : tables_)
        sorted.push_back({ctx, &ft});
    std::sort(sorted.begin(), sorted.end(), [](auto& a, auto& b){ return a.first < b.first; });

    for (const auto& [ctx, ftp] : sorted) {
        write32(ctx);
        for (int s = 0; s < alpha_; ++s)
            write16((uint16_t)ftp->syms[s].freq);
    }
    return out;
}

std::vector<uint8_t> ContextModel::serialize_topk(uint16_t sig_threshold) const {
    ARCS_CHECK(finalized_, "Must finalize() before serialize_topk()");

    std::vector<uint8_t> out;
    uint32_t n = (uint32_t)tables_.size();

    auto w32 = [&](uint32_t v) {
        out.push_back(v>>24); out.push_back(v>>16);
        out.push_back(v>>8);  out.push_back(v);
    };
    auto w16 = [&](uint16_t v) {
        out.push_back(v>>8); out.push_back(v);
    };

    w32(n);
    w16((uint16_t)alpha_);

    for (const auto& [ctx, ft] : tables_) {
        // Collect significant symbols (freq > threshold)
        std::vector<std::pair<uint16_t,uint8_t>> sig; // (freq, sym)
        for (int s = 0; s < alpha_; ++s) {
            uint16_t f = (uint16_t)ft.syms[s].freq;
            if (f > sig_threshold) sig.push_back({f, (uint8_t)s});
        }
        // Sort descending by freq for better LZMA compression of the list
        std::sort(sig.begin(), sig.end(), [](auto& a, auto& b){ return a.first > b.first; });

        w32(ctx);
        out.push_back((uint8_t)sig.size());
        for (auto& [f, s] : sig) {
            out.push_back(s);
            w16(f);
        }
    }
    return out;
}

void ContextModel::deserialize_topk(const uint8_t* data, size_t len) {
    const uint8_t* p   = data;
    const uint8_t* end = data + len;

    auto r32 = [&]() -> uint32_t {
        ARCS_CHECK(p + 4 <= end, "Truncated topk model");
        uint32_t v = ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|
                     ((uint32_t)p[2]<<8)|(uint32_t)p[3];
        p += 4; return v;
    };
    auto r16 = [&]() -> uint16_t {
        ARCS_CHECK(p + 2 <= end, "Truncated topk model");
        uint16_t v = ((uint16_t)p[0]<<8)|(uint16_t)p[1];
        p += 2; return v;
    };

    uint32_t n_ctx = r32();
    alpha_         = (int)r16();

    tables_.clear();
    counts_.clear();

    for (uint32_t c = 0; c < n_ctx; ++c) {
        ContextKey ctx = r32();
        uint8_t K = *p++;

        // Read K significant (sym, freq) pairs
        std::vector<uint32_t> freqs(alpha_, 0);
        uint32_t sum_sig = 0;
        for (int k = 0; k < K; ++k) {
            uint8_t sym  = *p++;
            uint16_t f   = r16();
            if (sym < (uint8_t)alpha_) { freqs[sym] = f; sum_sig += f; }
        }

        // Distribute residual probability uniformly among non-significant symbols
        int n_floor = 0;
        for (int s = 0; s < alpha_; ++s) if (freqs[s] == 0) ++n_floor;
        if (n_floor > 0) {
            uint32_t residual = SCALE - sum_sig;
            uint32_t floor_f  = residual / n_floor;
            uint32_t leftover = residual % n_floor;
            int assigned = 0;
            for (int s = 0; s < alpha_; ++s) {
                if (freqs[s] == 0) {
                    freqs[s] = floor_f + (assigned < (int)leftover ? 1 : 0);
                    ++assigned;
                }
            }
        } else if (sum_sig != SCALE) {
            // All symbols significant but sum mismatch — adjust last significant symbol
            int last_sig = alpha_ - 1;
            while (last_sig > 0 && freqs[last_sig] == 0) --last_sig;
            freqs[last_sig] += SCALE - sum_sig;
        }

        tables_[ctx].build_exact(freqs.data(), alpha_);
    }

    std::vector<uint32_t> uniform(alpha_, 1);
    default_table_.build(uniform.data(), alpha_);
    finalized_ = true;
}

void ContextModel::deserialize(const uint8_t* data, size_t len) {
    const uint8_t* p   = data;
    const uint8_t* end = data + len;

    auto read32 = [&]() -> uint32_t {
        ARCS_CHECK(p + 4 <= end, "Truncated model data");
        uint32_t v = ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|
                     ((uint32_t)p[2]<<8)|(uint32_t)p[3];
        p += 4; return v;
    };
    auto read16 = [&]() -> uint16_t {
        ARCS_CHECK(p + 2 <= end, "Truncated model data");
        uint16_t v = ((uint16_t)p[0]<<8)|(uint16_t)p[1];
        p += 2; return v;
    };

    uint32_t n_ctx = read32();
    alpha_         = read16();

    tables_.clear();
    counts_.clear();

    for (uint32_t c = 0; c < n_ctx; ++c) {
        ContextKey ctx = read32();
        std::vector<uint32_t> freqs(alpha_);
        for (int s = 0; s < alpha_; ++s)
            freqs[s] = read16();
        // Use build_exact: stored freqs already sum to SCALE. Re-normalising here
        // (via build()) would change the tables and desync the rANS decoder.
        tables_[ctx].build_exact(freqs.data(), alpha_);
    }

    // Build default uniform table
    std::vector<uint32_t> uniform(alpha_, 1);
    default_table_.build(uniform.data(), alpha_);
    finalized_ = true;
}
