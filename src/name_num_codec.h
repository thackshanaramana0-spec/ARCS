// ── Order-0 range coder for numeric name columns (our own, no LZMA) ───────────
// Read-name coordinates (Illumina tile/X/Y) are high-cardinality, near-random
// integers. LZMA is a dictionary+match coder — it cannot reach the symbol-level
// entropy of such columns (measured: ~487 KB on DS5 coords vs the ~365 KB order-0
// floor). We verified (held-out) there is NO exploitable sub-order-0 structure —
// X and Y are independent, no grid, no digit pattern — so an order-0 entropy coder
// on each column's value distribution IS optimal. This static-frequency range
// coder hits that floor and removes the external LZMA from the name path.
//
// Column payload:  [varint n][varint n_distinct][varint TOT]
//                  dict:  n_distinct × varint(value - prev)      (sorted distinct)
//                  freqs: n_distinct × varint(scaled_freq)       (sum == TOT)
//                  [range-coded stream of symbol IDs]
#pragma once
#include <vector>
#include <cstdint>
#include <algorithm>

namespace namenum {

constexpr uint32_t RC_TOP = 1u << 24;
constexpr uint32_t RC_BOT = 1u << 16;

struct REnc {
    uint32_t low = 0, range = 0xFFFFFFFFu;
    std::vector<uint8_t> out;
    inline void encode(uint32_t cum, uint32_t f, uint32_t tot) {
        range /= tot; low += cum * range; range *= f;
        while ((low ^ (low + range)) < RC_TOP ||
               (range < RC_BOT && ((range = (0u - low) & (RC_BOT - 1)), true))) {
            out.push_back((uint8_t)(low >> 24)); low <<= 8; range <<= 8;
        }
    }
    inline void flush() { for (int i = 0; i < 4; ++i) { out.push_back((uint8_t)(low >> 24)); low <<= 8; } }
};

struct RDec {
    uint32_t low = 0, range = 0xFFFFFFFFu, code = 0;
    const uint8_t* p; const uint8_t* e;
    RDec(const uint8_t* P, const uint8_t* E) : p(P), e(E) {
        for (int i = 0; i < 4; ++i) code = (code << 8) | (p < e ? *p++ : 0u);
    }
    inline uint32_t getfreq(uint32_t tot) { range /= tot; return (code - low) / range; }
    inline void update(uint32_t cum, uint32_t f) {
        low += cum * range; range *= f;
        while ((low ^ (low + range)) < RC_TOP ||
               (range < RC_BOT && ((range = (0u - low) & (RC_BOT - 1)), true))) {
            code = (code << 8) | (p < e ? *p++ : 0u); low <<= 8; range <<= 8;
        }
    }
};

inline void put_varint(std::vector<uint8_t>& b, uint64_t v) {
    for (;;) { uint8_t x = (uint8_t)(v & 0x7f); v >>= 7; if (v) b.push_back(x | 0x80); else { b.push_back(x); break; } }
}
inline uint64_t get_varint(const uint8_t*& p, const uint8_t* e) {
    uint64_t v = 0; int s = 0;
    while (p < e) { uint8_t b = *p++; v |= (uint64_t)(b & 0x7f) << s; if (!(b & 0x80)) break; s += 7; }
    return v;
}

// Encode a numeric column (non-negative int64 values). Returns the payload, or an
// empty vector if the column is degenerate (caller should keep-smaller anyway).
inline std::vector<uint8_t> encode_column(const std::vector<int64_t>& vals) {
    const size_t n = vals.size();
    std::vector<uint8_t> out;
    if (n == 0) return out;
    // Distinct sorted values → dense symbol IDs.
    std::vector<int64_t> dv(vals);
    std::sort(dv.begin(), dv.end());
    dv.erase(std::unique(dv.begin(), dv.end()), dv.end());
    const size_t D = dv.size();
    // The range coder needs total ≤ RC_BOT so range/tot ≥ 1 (range ≥ RC_BOT after
    // renorm). With one freq ≥ 1 per symbol that requires D ≤ RC_BOT. Columns with
    // more distinct values than that (e.g. a unique per-read index) are not suited to
    // a single static order-0 table — bail so the caller keeps the varint form.
    if (D >= (size_t)RC_BOT) return out;
    // counts per id
    std::vector<uint32_t> cnt(D, 0);
    auto id_of = [&](int64_t v) { return (uint32_t)(std::lower_bound(dv.begin(), dv.end(), v) - dv.begin()); };
    for (int64_t v : vals) cnt[id_of(v)]++;
    // Normalize counts to total TOT = RC_BOT (each freq ≥ 1; D ≤ RC_BOT guarantees room).
    const uint32_t TOT = RC_BOT;
    std::vector<uint32_t> f(D, 0);
    uint64_t acc = 0;
    for (size_t i = 0; i < D; ++i) {
        uint32_t s = (uint32_t)(((uint64_t)cnt[i] * TOT) / n);
        if (s == 0) s = 1;
        f[i] = s; acc += s;
    }
    // Fix rounding so sum(f) == TOT exactly (adjust the most-frequent symbol).
    if (acc != TOT) {
        size_t mx = 0; for (size_t i = 1; i < D; ++i) if (f[i] > f[mx]) mx = i;
        int64_t diff = (int64_t)TOT - (int64_t)acc;
        if ((int64_t)f[mx] + diff < 1) {           // fall back: rescale evenly
            // extremely skewed edge — recompute with larger TOT headroom
            f.assign(D, 1); uint64_t rem = TOT - D;
            for (size_t i = 0; i < D && rem > 0; ++i) { uint64_t add = (uint64_t)cnt[i] * rem / n; f[i] += (uint32_t)add; }
            uint64_t s2 = 0; for (uint32_t v : f) s2 += v;
            size_t m2 = 0; for (size_t i = 1; i < D; ++i) if (f[i] > f[m2]) m2 = i;
            f[m2] += (uint32_t)(TOT - s2);
        } else {
            f[mx] = (uint32_t)((int64_t)f[mx] + diff);
        }
    }
    // Cumulative
    std::vector<uint32_t> cum(D + 1, 0);
    for (size_t i = 0; i < D; ++i) cum[i + 1] = cum[i] + f[i];
    // Header
    put_varint(out, n);
    put_varint(out, D);
    put_varint(out, TOT);
    int64_t prev = 0;
    for (int64_t v : dv) { put_varint(out, (uint64_t)(v - prev)); prev = v; }   // dict (sorted → nonneg deltas)
    for (uint32_t v : f) put_varint(out, v);                                     // freqs
    // Range-code the symbol IDs
    REnc enc;
    for (int64_t v : vals) { uint32_t id = id_of(v); enc.encode(cum[id], f[id], TOT); }
    enc.flush();
    out.insert(out.end(), enc.out.begin(), enc.out.end());
    return out;
}

// Decode a numeric column payload written by encode_column. Returns the values.
inline std::vector<int64_t> decode_column(const uint8_t* p, const uint8_t* e) {
    std::vector<int64_t> vals;
    uint64_t n = get_varint(p, e);
    uint64_t D = get_varint(p, e);
    uint32_t TOT = (uint32_t)get_varint(p, e);
    if (n == 0 || D == 0) return vals;
    std::vector<int64_t> dict(D);
    int64_t prev = 0;
    for (uint64_t i = 0; i < D; ++i) { prev += (int64_t)get_varint(p, e); dict[i] = prev; }
    std::vector<uint32_t> f(D);
    for (uint64_t i = 0; i < D; ++i) f[i] = (uint32_t)get_varint(p, e);
    std::vector<uint32_t> cum(D + 1, 0);
    for (uint64_t i = 0; i < D; ++i) cum[i + 1] = cum[i] + f[i];
    // slot → id table for O(1) symbol lookup
    std::vector<uint32_t> slot(TOT);
    for (uint64_t i = 0; i < D; ++i)
        for (uint32_t s = cum[i]; s < cum[i + 1]; ++s) slot[s] = (uint32_t)i;
    RDec dec(p, e);
    vals.resize((size_t)n);
    for (uint64_t k = 0; k < n; ++k) {
        uint32_t ff = dec.getfreq(TOT);
        if (ff >= TOT) ff = TOT - 1;
        uint32_t id = slot[ff];
        dec.update(cum[id], f[id]);
        vals[(size_t)k] = dict[id];
    }
    return vals;
}

} // namespace namenum
