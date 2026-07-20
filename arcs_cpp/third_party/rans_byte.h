// rans_byte.h — rANS (range Asymmetric Numeral Systems) byte-precision codec
// Based on Fabian Giesen's public domain implementation.
// Extended for context-model integration and streaming encode/decode.
// Original: https://github.com/rygorous/ryg_rans
// License: Public domain

#pragma once
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <vector>
#include <algorithm>

// ── rANS state ────────────────────────────────────────────────────────────────
// State is a 64-bit integer. Lower 32 bits hold the ANS state.
// Upper 32 bits hold buffered output bytes.

static constexpr uint32_t RANS_BYTE_L = (1u << 23);  // lower bound on normalized state
static constexpr int      SCALE_BITS  = 14;           // freq table precision: 1/2^14
static constexpr uint32_t SCALE       = (1u << SCALE_BITS);

using RansState = uint64_t;

// ── Symbol statistics (cumulative frequency table) ───────────────────────────
struct RansSymbol {
    uint32_t start; // cumulative frequency (sum of freqs of all prior symbols)
    uint32_t freq;  // frequency of this symbol (in [1, SCALE])
};

// ── Encoder (outputs bytes in reverse order) ──────────────────────────────────
struct RansEncoder {
    RansState state = RANS_BYTE_L;

    void init() { state = RANS_BYTE_L; }

    // Encode one symbol. Output bytes pushed to 'out' (reversed on flush).
    inline void encode(std::vector<uint8_t>& out, const RansSymbol& sym) {
        assert(sym.freq > 0 && sym.start + sym.freq <= SCALE);
        // Renormalize: output bytes until state is in valid range
        uint32_t freq = sym.freq;
        uint64_t x = state;
        // Upper bound for this symbol: x < (RANS_BYTE_L / freq) * SCALE * 256
        uint64_t x_max = ((RANS_BYTE_L >> SCALE_BITS) << 8) * freq;
        while (x >= x_max) {
            out.push_back((uint8_t)(x & 0xFF));
            x >>= 8;
        }
        // Encode: state = (x / freq) * SCALE + (x % freq) + sym.start
        state = ((x / freq) << SCALE_BITS) + (x % freq) + sym.start;
    }

    // Flush remaining state (4 bytes, MSB first)
    inline void flush(std::vector<uint8_t>& out) {
        uint32_t s = (uint32_t)state;
        out.push_back((uint8_t)(s >> 24));
        out.push_back((uint8_t)(s >> 16));
        out.push_back((uint8_t)(s >>  8));
        out.push_back((uint8_t)(s >>  0));
    }
};

// ── Decoder (reads bytes in forward order) ────────────────────────────────────
struct RansDecoder {
    RansState state = 0;

    // Initialize from byte stream (reads 4 bytes from ptr).
    // After the encoder's flush(MSB-first push) + full-buffer reverse, the state
    // bytes appear little-endian at the front of the stream — read them as such.
    inline void init(const uint8_t* ptr) {
        state  = (uint32_t)ptr[0];
        state |= (uint32_t)ptr[1] << 8;
        state |= (uint32_t)ptr[2] << 16;
        state |= (uint32_t)ptr[3] << 24;
    }

    // Get cumulative frequency of current symbol (before decoding)
    inline uint32_t get_cumfreq() const {
        return (uint32_t)(state & (SCALE - 1));
    }

    // Advance decoder after decoding a symbol
    inline void advance(const uint8_t*& ptr, const uint8_t* end,
                        const RansSymbol& sym) {
        uint32_t freq  = sym.freq;
        uint32_t cum   = sym.start;
        uint64_t x     = state;
        // Decode: x = freq * (x >> SCALE_BITS) + (x & (SCALE-1)) - cum
        x = freq * (x >> SCALE_BITS) + (x & (SCALE-1)) - cum;
        // Renormalize: read bytes until state is back in [RANS_BYTE_L, RANS_BYTE_L*256)
        while (x < RANS_BYTE_L && ptr < end) {
            x = (x << 8) | *ptr++;
        }
        state = x;
    }
};

// ── Frequency table builder ───────────────────────────────────────────────────
// Converts raw counts to SCALE-summing frequencies with Laplace smoothing.
struct FreqTable {
    static constexpr int MAX_SYMBOLS = 256;

    std::vector<RansSymbol> syms; // indexed by symbol value
    std::vector<uint8_t>    alias; // SCALE-length table: alias[cum_freq] → symbol

    // Build directly from already-normalised frequencies that sum to exactly SCALE.
    // Does NOT re-normalise — used by deserialize() so decode tables exactly match
    // the encode tables (any re-normalisation here would desync the rANS decoder).
    void build_exact(const uint32_t* freqs, int n_syms) {
        assert(n_syms > 0 && n_syms <= MAX_SYMBOLS);
        syms.resize(n_syms);
        alias.resize(SCALE);
        uint32_t cum = 0;
        for (int i = 0; i < n_syms; ++i) {
            syms[i] = {cum, freqs[i]};
            for (uint32_t c = cum; c < cum + freqs[i]; ++c)
                alias[c] = (uint8_t)i;
            cum += freqs[i];
        }
        assert(cum == SCALE);
    }

    // Build from raw counts. n_syms = alphabet size.
    // Laplace smoothing: each symbol gets at least 1.
    void build(const uint32_t* counts, int n_syms) {
        assert(n_syms > 0 && n_syms <= MAX_SYMBOLS);
        syms.resize(n_syms);
        alias.resize(SCALE);

        // Laplace smoothing (+1 to each count)
        uint64_t total = 0;
        std::vector<uint64_t> smoothed(n_syms);
        for (int i = 0; i < n_syms; ++i) {
            smoothed[i] = (uint64_t)counts[i] + 1;
            total += smoothed[i];
        }

        // Compute floor-quantised frequencies; each gets at least 1.
        // Track actual sum — it may be > SCALE due to floor()+max(1) interaction.
        std::vector<uint32_t> freqs(n_syms);
        uint32_t total_freq = 0;
        for (int i = 0; i < n_syms; ++i) {
            freqs[i] = std::max(1u, (uint32_t)((smoothed[i] * SCALE) / total));
            total_freq += freqs[i];
        }

        // Adjust sum to exactly SCALE.
        // Case 1: total_freq < SCALE — distribute positive remainder.
        // Case 2: total_freq > SCALE — reduce largest freqs (must stay >= 1).
        if (total_freq < SCALE) {
            uint32_t rem = SCALE - total_freq;
            for (int i = 0; rem > 0; i = (i + 1) % n_syms) {
                ++freqs[i]; --rem;
            }
        } else if (total_freq > SCALE) {
            uint32_t excess = total_freq - SCALE;
            // Reduce from symbols with highest freq (sorted by freq desc)
            std::vector<int> order(n_syms);
            for (int i = 0; i < n_syms; ++i) order[i] = i;
            std::sort(order.begin(), order.end(),
                      [&](int a, int b){ return freqs[a] > freqs[b]; });
            for (int round = 0; excess > 0; round = (round + 1) % n_syms) {
                int idx = order[round % n_syms];
                if (freqs[idx] > 1) { --freqs[idx]; --excess; }
            }
        }

        // Build cumulative table and alias lookup
        uint32_t cum = 0;
        for (int i = 0; i < n_syms; ++i) {
            syms[i] = {cum, freqs[i]};
            assert(cum + freqs[i] <= SCALE); // guarantee no overflow
            for (uint32_t c = cum; c < cum + freqs[i]; ++c)
                alias[c] = (uint8_t)i;
            cum += freqs[i];
        }
        assert(cum == SCALE);
    }

    // Decode: given cumulative freq, return symbol
    uint8_t decode_sym(uint32_t cumfreq) const {
        assert(cumfreq < SCALE);
        return alias[cumfreq];
    }
};

// ── Encode/decode helpers ─────────────────────────────────────────────────────

// Encode a block of symbols using one frequency table.
// Output is byte-reversed; caller must reverse before writing to stream.
inline std::vector<uint8_t> rans_encode_block(
    const uint8_t* symbols, size_t n,
    const FreqTable& ft)
{
    RansEncoder enc;
    std::vector<uint8_t> out;
    out.reserve(n + 8);
    // Encode in reverse order (rANS is a stack)
    for (int i = (int)n - 1; i >= 0; --i)
        enc.encode(out, ft.syms[symbols[i]]);
    enc.flush(out);
    std::reverse(out.begin(), out.end());
    return out;
}

// Decode a block of n symbols from byte stream.
inline std::vector<uint8_t> rans_decode_block(
    const uint8_t* data, size_t data_len,
    const FreqTable& ft,
    size_t n_symbols)
{
    std::vector<uint8_t> out(n_symbols);
    RansDecoder dec;
    const uint8_t* ptr = data;
    const uint8_t* end = data + data_len;
    dec.init(ptr);
    ptr += 4;
    for (size_t i = 0; i < n_symbols; ++i) {
        uint32_t cf  = dec.get_cumfreq();
        uint8_t  sym = ft.decode_sym(cf);
        out[i]       = sym;
        dec.advance(ptr, end, ft.syms[sym]);
    }
    return out;
}
