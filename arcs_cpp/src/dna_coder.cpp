// ARCS-DNA: pseudogenome-native DNA compressor.
//
// Method: multi-order finite-context models combined by an integer logistic
// mixer (PAQ/cmix-style context mixing), coded with a carryless binary
// arithmetic coder. Bases are split into 2 binary decisions and predicted
// bit-by-bit.
//
// Why logistic mixing is the right tool here (measured, not assumed):
//   An unseen high-order context predicts p=0.5 → stretch(0.5)=0 → it adds
//   nothing to the logit sum and is auto-ignored. A seen context predicts a
//   confident p → large stretch → it dominates. This gives per-position model
//   gating for free — precisely what a single slow global weight could not do.
//
// All arithmetic is integer → encode and decode are bit-exact and portable.
//
// Public API (dna_coder.h): dna_encode / dna_decode.
// reads / chain_order / pg_pos are accepted for a later pre-seeding stage;
// unused here so that decode (which lacks them) stays perfectly symmetric.

#include "dna_coder.h"
#include <array>
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <cassert>
#include <unordered_map>

// ── Base helpers ──────────────────────────────────────────────────────────────
static inline int b2i(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:             return 0;   // N → A
    }
}
static inline char i2b(int b) { return "ACGT"[b & 3]; }

// ── squash / stretch (PAQ canonical, 12-bit probabilities) ────────────────────
// squash: logit d ∈ [-2047,2047] → probability ∈ [0,4095].
// stretch: inverse, probability ∈ [0,4095] → logit ∈ [-2047,2047].
// Thread-local so parallel chunk encode/decode (chunked mode) don't race on the
// lazy table build. g_stretch is constant once built; each thread builds its own.
static thread_local int  g_stretch[4096];
static thread_local bool g_tables_ready = false;

static int squash(int d) {
    static const int t[33] = {
        1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,1546,2047,
        2549,2994,3348,3607,3785,3901,3975,4022,4050,4068,4079,4085,4089,4092,4093,4094
    };
    if (d >  2047) return 4095;
    if (d < -2047) return 0;
    int w = d & 127;
    d = (d >> 7) + 16;
    return (t[d] * (128 - w) + t[d + 1] * w + 64) >> 7;
}

static void init_tables() {
    if (g_tables_ready) return;
    // Build stretch by inverting squash (monotonic).
    int pi = 0;
    for (int x = -2047; x <= 2047; ++x) {
        int p = squash(x);
        for (int j = pi; j <= p; ++j) g_stretch[j] = x;
        pi = p + 1;
    }
    for (int j = pi; j < 4096; ++j) g_stretch[j] = 2047;
    g_tables_ready = true;
}
static inline int stretch(int p) { return g_stretch[p]; }

// ── Carryless binary arithmetic coder (PAQ-style) ─────────────────────────────
// Probabilities are 12-bit P(bit==1) in [1,4095].
struct BinEncoder {
    uint32_t x1 = 0, x2 = 0xFFFFFFFFu;
    std::vector<uint8_t> out;

    void encode(int bit, int p1) {          // p1 = P(bit==1), 12-bit
        if (p1 < 1) p1 = 1; else if (p1 > 4095) p1 = 4095;
        uint32_t xmid = x1 + (uint32_t)(((uint64_t)(x2 - x1) * (uint32_t)p1) >> 12);
        if (bit) x2 = xmid; else x1 = xmid + 1;
        // Renormalize: flush bytes that have converged.
        while (((x1 ^ x2) & 0xFF000000u) == 0) {
            out.push_back((uint8_t)(x2 >> 24));
            x1 <<= 8;
            x2 = (x2 << 8) | 0xFF;
        }
    }

    void flush() {
        // Emit 4 bytes of x1 to disambiguate the final interval.
        for (int i = 0; i < 4; ++i) {
            out.push_back((uint8_t)(x1 >> 24));
            x1 <<= 8;
        }
    }
};

struct BinDecoder {
    uint32_t x1 = 0, x2 = 0xFFFFFFFFu, x = 0;
    const uint8_t* ptr;
    const uint8_t* end;

    BinDecoder(const uint8_t* p, const uint8_t* e) : ptr(p), end(e) {
        for (int i = 0; i < 4; ++i)
            x = (x << 8) | (ptr < end ? *ptr++ : 0u);
    }

    int decode(int p1) {                    // p1 = P(bit==1), 12-bit
        if (p1 < 1) p1 = 1; else if (p1 > 4095) p1 = 4095;
        uint32_t xmid = x1 + (uint32_t)(((uint64_t)(x2 - x1) * (uint32_t)p1) >> 12);
        int bit;
        if (x <= xmid) { bit = 1; x2 = xmid; }
        else           { bit = 0; x1 = xmid + 1; }
        while (((x1 ^ x2) & 0xFF000000u) == 0) {
            x1 <<= 8;
            x2 = (x2 << 8) | 0xFF;
            x = (x << 8) | (ptr < end ? *ptr++ : 0u);
        }
        return bit;
    }
};

// ── Finite-context models ─────────────────────────────────────────────────────
// Each model keeps 4 counts per context. Low orders use a dense table indexed
// by the packed context; high orders use a fixed-size direct-mapped hash table
// (collisions share stats — bounded RAM, same idea as GeCo3's cache-hash).
// Counts are 8-bit: they rescale (halve) once any reaches COUNT_RESCALE (≤255) and
// so never exceed 255, so uint8 holds identical values to the old uint16 — the coder
// is BIT-IDENTICAL — while HALVING every model table's memory (the pg model tables
// are the dominant RAM user on large pseudogenomes). COUNT_RESCALE is clamped to 255.
using Row4 = std::array<uint8_t, 4>;
static thread_local uint16_t COUNT_RESCALE = 255;   // lower = faster count aging; ≤255

// Smoothing α expressed as ALPHA_NUM/16 (default 8 → α=0.5). Smaller α makes
// high-order models bet more confidently when they have evidence.
static thread_local int ALPHA_NUM = 4;

// 12-bit P(bit==1) for a node from a count row.
//   node 0 : P(high bit==1) = P(base ∈ {G,T})
//   node 1 : P(low bit==1 | high==0) = P(base==C | base ∈ {A,C})
//   node 2 : P(low bit==1 | high==1) = P(base==T | base ∈ {G,T})
static inline int predict_row(const Row4& c, int node, int hb) {
    int n1, n0;
    if (node == 0)      { n1 = c[2] + c[3]; n0 = c[0] + c[1]; }
    else if (hb == 0)   { n1 = c[1];        n0 = c[0]; }
    else                { n1 = c[3];        n0 = c[2]; }
    // p = (n1 + α) / (n0 + n1 + 2α), α = ALPHA_NUM/16, scaled to 12 bits.
    // (A magic-reciprocal table for this division was tried and REVERTED: den varies
    // per call, so a 180 KB lookup table thrashed the model-table cache — net slower.
    // Unlike build_cf, where the divisor is invariant per symbol. See memory.)
    int num = (n1 * 16 + ALPHA_NUM) * 4096;
    int den = (n0 + n1) * 16 + 2 * ALPHA_NUM;
    int p = num / den;
    if (p < 1) p = 1; else if (p > 4095) p = 4095;
    return p;
}
static inline int argmax_row(const Row4& c) {
    int best = 0;
    for (int j = 1; j < 4; ++j) if (c[j] > c[best]) best = j;
    return best;
}

struct Model {
    int      order;      // context length in bases
    int      bits;       // table index bits
    bool     dense;      // dense (indexed by context) vs hashed
    uint64_t mask;       // (1 << (2*order)) - 1
    std::vector<Row4> tbl;
    Row4*    cur = nullptr;   // row selected for the current position
    int      ir_ctx_shift = 0;   // rcomp >> this = RC context (see ir_update_fast)
    int      ir_sym_shift = 0;   // rcomp >> this & 3 = RC symbol

    static thread_local int hash_bits;           // hashed-model table size (2^bits); TL for parallel chunks
    void init(int ord) {
        order = ord;
        mask  = (ord >= 32) ? ~0ULL : ((1ULL << (2 * ord)) - 1);
        dense = (ord <= 8);
        bits  = dense ? (2 * ord) : hash_bits;
        tbl.assign((size_t)1 << bits, Row4{0, 0, 0, 0});
        // rcomp holds 32 complemented bases with the newest at bits 62-63. The top
        // k bases (shifted down by 2*(32-k)) form the RC context, and the base just
        // below the k-window (2*(31-k)) is the RC symbol. Valid for order ≤ 31.
        ir_ctx_shift = (ord <= 31) ? 2 * (32 - ord) : 0;
        ir_sym_shift = (ord <= 31) ? 2 * (31 - ord) : 0;
    }

    inline size_t index(uint64_t hist) const {
        uint64_t ctx = hist & mask;
        if (dense) return (size_t)ctx;
        uint64_t h = (ctx + 1) * 0x9E3779B97F4A7C15ULL;
        return (size_t)(h >> (64 - bits));
    }

    // Select the row for this position's context (called once per base). The
    // prefetch issues the (likely cache-missing) table load early so the 8+
    // models' misses overlap while the mixer consumes the first prediction.
    inline void select(uint64_t hist) { cur = &tbl[index(hist)]; __builtin_prefetch(cur, 1, 1); }

    inline int predict(int node, int hb) const { return predict_row(*cur, node, hb); }

    // Register the observed base into the selected row.
    inline void update(int b) {
        Row4& c = *cur;
        if (++c[b] >= COUNT_RESCALE) {
            c[0] >>= 1; c[1] >>= 1; c[2] >>= 1; c[3] >>= 1;
        }
    }

    // Inverted-repeat update: record the reverse complement of the observed
    // (k+1)-mer so that repeats on the opposite strand also train this table.
    // Observed (k+1)-mer x[i-k..i] (context = last k bases, symbol = b); its RC
    // splits into rc_context = comp(x[i]..x[i-k+1]) predicting rc_symbol =
    // comp(x[i-k]). hist holds x[i-1](LSB)…x[i-k] before b is folded in.
    inline void ir_update(uint64_t hist, int b) {
        int k = order;
        uint64_t rc = (uint64_t)(b ^ 3);
        for (int j = 1; j <= k - 1; ++j)
            rc = (rc << 2) | (uint64_t)((int)((hist >> (2 * (j - 1))) & 3) ^ 3);
        int rc_sym = (int)((hist >> (2 * (k - 1))) & 3) ^ 3;
        Row4& c = tbl[index(rc)];
        if (++c[rc_sym] >= COUNT_RESCALE) {
            c[0] >>= 1; c[1] >>= 1; c[2] >>= 1; c[3] >>= 1;
        }
    }

    // Incremental IR update — O(1) instead of O(order). `rcomp` is the rolling
    // reverse-complement of history: rcomp = (rcomp>>2) | (comp(x[i])<<62), newest
    // base at the top. Shifting down by ir_ctx_shift lands the top k complemented
    // bases in the same bit layout the per-base loop in ir_update builds, so the
    // table index and symbol are value-identical to ir_update for i ≥ order. (Warm-up
    // positions where x[i-k] does not exist differ harmlessly and stay consistent
    // between encoder and decoder.)
    inline void ir_update_fast(uint64_t rcomp) {
        Row4& c = tbl[index(rcomp >> ir_ctx_shift)];
        int rc_sym = (int)((rcomp >> ir_sym_shift) & 3);
        if (++c[rc_sym] >= COUNT_RESCALE) {
            c[0] >>= 1; c[1] >>= 1; c[2] >>= 1; c[3] >>= 1;
        }
    }
};
thread_local int Model::hash_bits = 22;

// ── Substitution-tolerant context model (STCM) ────────────────────────────────
// Shares a base FCM's count table but reads it through a "tolerant" context
// pointer that follows a previously-seen near-identical region across
// substitutions. On a mismatch it advances the context with the model's own
// prediction (staying aligned to the reference copy) and counts the edit; when
// edits exceed a threshold it resyncs to the true recent history. This captures
// the substitutional nature of genomic near-repeats (sequencing errors, SNPs).
struct STCM {
    Model* base = nullptr;
    int    threshold = 2;
    uint64_t tctx = 0;
    int    edits = 0;
    Row4*  cur = nullptr;

    void attach(Model* m, int thr) { base = m; threshold = thr; tctx = 0; edits = 0; }
    inline void select() { cur = &base->tbl[base->index(tctx)]; }
    inline int  predict(int node, int hb) const { return predict_row(*cur, node, hb); }

    inline void update(int b, uint64_t hist) {
        int best = argmax_row(*cur);
        if (b == best) { if (edits > 0) --edits; }
        else           { ++edits; }
        if (edits > threshold) {
            tctx = hist & base->mask;          // resync to reality
            edits = 0;
        } else {
            int adv = (b == best) ? b : best;  // follow the reference copy
            tctx = ((tctx << 2) | (uint64_t)adv) & base->mask;
        }
    }
};

// ── Logistic mixer (lpaq-style integer mixer, context-selected weights) ───────
// Weight set selected by an arbitrary context index. Inputs are stretched
// model predictions. Each bit picks one weight row of length n.
struct Mixer {
    int n;                          // number of models (inputs)
    int nsets;                      // number of selectable weight rows
    std::vector<int16_t> w;         // nsets × n weights
    int  in[64];                    // stretched inputs for current bit
    int  cnt = 0;
    int  sel = 0;                   // selected row
    int  pr = 2048;

    void init(int nmodels, int sets) {
        n = nmodels; nsets = sets;
        w.assign((size_t)nsets * n, 0);
    }
    inline void set_ctx(int s) { sel = s; cnt = 0; }
    inline void add(int stretched_p) { in[cnt++] = stretched_p; }

    inline int mix() {
        const int16_t* wr = &w[(size_t)sel * n];
        int64_t dot = 0;
        for (int i = 0; i < cnt; ++i) dot += (int64_t)in[i] * wr[i];
        int d = (int)(dot >> 8);
        if (d >  2047) d =  2047; else if (d < -2047) d = -2047;
        pr = squash(d);
        if (pr < 1) pr = 1; else if (pr > 4095) pr = 4095;
        return pr;
    }

    int lr = 7;
    inline void update(int bit) {
        int err = ((bit << 12) - pr) * lr;
        int16_t* wr = &w[(size_t)sel * n];
        for (int i = 0; i < cnt; ++i) {
            int nw = wr[i] + ((in[i] * err + 0x8000) >> 16);
            if (nw < -32768) nw = -32768; else if (nw > 32767) nw = 32767;
            wr[i] = (int16_t)nw;
        }
    }
};

// ── APM / SSE (adaptive probability map, lpaq2-style) ─────────────────────────
// Refines a probability using a small context by interpolating a learned curve.
struct APM {
    std::vector<uint16_t> t;
    int idx = 0;
    void init(int nctx) {
        t.resize((size_t)nctx * 33);
        for (int i = 0; i < nctx; ++i)
            for (int j = 0; j < 33; ++j)
                t[(size_t)i * 33 + j] = (uint16_t)(squash((j - 16) * 128) * 16);
    }
    // Refine pr (12-bit) under context cx.
    inline int pp(int pr, int cx) {
        int s = stretch(pr);                 // -2047..2047
        int w = s & 127;
        idx = ((s >> 7) + 16) + cx * 33;
        return (t[idx] * (128 - w) + t[idx + 1] * w) >> 11;
    }
    inline void update(int bit) {
        int g = (bit << 16) + (bit << 7) - bit - bit;
        t[idx]     = (uint16_t)(t[idx]     + ((g - t[idx])     >> 7));
        t[idx + 1] = (uint16_t)(t[idx + 1] + ((g - t[idx + 1]) >> 7));
    }
};

// ── Model configuration ───────────────────────────────────────────────────────
// Tunable via env vars (ARCSDNA_ORDERS, ARCSDNA_MIXBITS, ARCSDNA_APMW,
// ARCSDNA_APMBITS, ARCSDNA_LR) for experimentation. Encode and decode read the
// same config within a run. Defaults are the shipped configuration.
#include <cstdlib>
#include <sstream>

// Shipped defaults are the empirically-tuned winning configuration. Encode and
// decode both use these when no env override is set, so separate compress /
// decompress processes stay consistent. Env vars exist only for research.
struct DnaConfig {
    std::vector<int> orders = {2, 4, 8, 11, 14, 18, 22, 26};
    std::vector<std::pair<int,int>> stcm;  // (order, edit-threshold) tolerant models
    int mix_bits = 2;    // history bits (beyond node) selecting mixer weight row
    int apm_w    = 2;    // APM blend weight numerator (0 = APM off)
    int apm_bits = 8;    // history bits selecting APM context
    int lr       = 1;    // mixer learning rate
    int ir       = 1;    // inverted-repeat training on all models (0/1)
};

static thread_local DnaConfig g_cfg;
static thread_local bool      g_cfg_loaded = false;
static thread_local bool      g_hbits_from_env = false;

// Auto-size hashed-model tables to the pg length (unless overridden by env).
// Load factor stays low so high-order contexts rarely collide; both encode and
// decode know N, so they pick identical sizes.
static void autosize_hash_bits(size_t N) {
    if (g_hbits_from_env) return;
    int b = 18;
    while (((size_t)1 << b) < 2 * N && b < 24) ++b;  // ≥2× N entries, capped at 2^24
    Model::hash_bits = b;
}

static void load_cfg() {
    if (g_cfg_loaded) return;
    if (const char* s = getenv("ARCSDNA_ORDERS")) {
        std::vector<int> v; std::stringstream ss(s); std::string tok;
        while (std::getline(ss, tok, ',')) if (!tok.empty()) v.push_back(atoi(tok.c_str()));
        if (!v.empty()) g_cfg.orders = v;
    }
    if (const char* s = getenv("ARCSDNA_STCM")) {
        // format: "order:thresh,order:thresh" e.g. "14:2,18:3,22:4"
        std::vector<std::pair<int,int>> v; std::stringstream ss(s); std::string tok;
        while (std::getline(ss, tok, ',')) {
            size_t c = tok.find(':');
            if (c != std::string::npos)
                v.emplace_back(atoi(tok.substr(0,c).c_str()), atoi(tok.substr(c+1).c_str()));
        }
        g_cfg.stcm = v;
    }
    if (const char* s = getenv("ARCSDNA_MIXBITS")) g_cfg.mix_bits = atoi(s);
    if (const char* s = getenv("ARCSDNA_APMW"))    g_cfg.apm_w    = atoi(s);
    if (const char* s = getenv("ARCSDNA_APMBITS")) g_cfg.apm_bits = atoi(s);
    if (const char* s = getenv("ARCSDNA_LR"))      g_cfg.lr       = atoi(s);
    if (const char* s = getenv("ARCSDNA_IR"))      g_cfg.ir       = atoi(s);
    if (const char* s = getenv("ARCSDNA_RESCALE")) { int v = atoi(s); if (v < 2) v = 2; if (v > 255) v = 255; COUNT_RESCALE = (uint16_t)v; }
    if (const char* s = getenv("ARCSDNA_HBITS")) { Model::hash_bits = atoi(s); g_hbits_from_env = true; }
    if (const char* s = getenv("ARCSDNA_ALPHA"))   ALPHA_NUM = atoi(s);
    g_cfg_loaded = true;
}

static inline int mixer_sel(int node, uint64_t hist, int mix_bits) {
    int m = (1 << mix_bits) - 1;
    return node * (1 << mix_bits) + (int)(hist & m);
}
static inline int apm_sel(int node, uint64_t hist, int apm_bits) {
    int m = (1 << apm_bits) - 1;
    return node * (1 << apm_bits) + (int)(hist & m);
}

// ── Shared predictive core ────────────────────────────────────────────────────
// One implementation drives both encode and decode; the only difference is the
// per-bit coding step, injected via BitCoder. This guarantees the encoder and
// decoder walk identical model/mixer/APM state.
namespace {

struct EncBitCoder {
    BinEncoder* e;
    inline int code(int p, int known_bit) { e->encode(known_bit, p); return known_bit; }
};
struct DecBitCoder {
    BinDecoder* d;
    inline int code(int p, int /*ignored*/) { return d->decode(p); }
};

// Build models + STCMs from config. STCM shares the table of an FCM of the same
// order (created if absent). Pointers are stable because all models are created
// before any STCM attaches.
struct Engine {
    std::vector<Model> models;
    std::vector<STCM>  stcms;
    Mixer mixer;
    APM   apm;
    int   n_inputs, MB, AB, AW;

    void build() {
        for (int o : g_cfg.orders) { models.emplace_back(); models.back().init(o); }
        std::vector<int> base_idx;
        for (auto& sc : g_cfg.stcm) {
            int idx = -1;
            for (int i = 0; i < (int)models.size(); ++i)
                if (models[i].order == sc.first) { idx = i; break; }
            if (idx < 0) { models.emplace_back(); models.back().init(sc.first); idx = (int)models.size()-1; }
            base_idx.push_back(idx);
        }
        stcms.resize(g_cfg.stcm.size());
        for (size_t k = 0; k < g_cfg.stcm.size(); ++k)
            stcms[k].attach(&models[base_idx[k]], g_cfg.stcm[k].second);

        n_inputs = (int)(models.size() + stcms.size());
        MB = g_cfg.mix_bits; AB = g_cfg.apm_bits; AW = g_cfg.apm_w;
        mixer.init(n_inputs, 3 * (1 << MB)); mixer.lr = g_cfg.lr;
        apm.init(3 * (1 << AB));
    }
};

template<class BitCoder>
inline int code_node(Engine& e, BitCoder& coder, int node, int hb,
                      uint64_t hist, int known_bit) {
    e.mixer.set_ctx(mixer_sel(node, hist, e.MB));
    for (auto& m : e.models) e.mixer.add(stretch(m.predict(node, hb)));
    for (auto& s : e.stcms)  e.mixer.add(stretch(s.predict(node, hb)));
    int pm = e.mixer.mix();
    int pc = e.AW ? (pm * (4 - e.AW) + e.apm.pp(pm, apm_sel(node, hist, e.AB)) * e.AW) >> 2 : pm;
    int bit = coder.code(pc, known_bit);
    e.mixer.update(bit);
    if (e.AW) e.apm.update(bit);
    return bit;
}

template<bool ENCODE, class BitCoder>
void run_core(const std::string* pin, std::string* pout, size_t N, BitCoder& coder) {
    Engine e; e.build();
    // O(1) IR needs order ≤ 31 (rcomp holds 32 bases); fall back otherwise.
    bool use_fast_ir = true;
    for (auto& m : e.models) if (m.order > 31) use_fast_ir = false;
    uint64_t hist  = 0;
    uint64_t rcomp = 0;   // rolling reverse-complement of history for O(1) IR update
    for (size_t i = 0; i < N; ++i) {
        int b_known = ENCODE ? b2i((*pin)[i]) : 0;

        for (auto& m : e.models) m.select(hist);
        for (auto& s : e.stcms)  s.select();

        int hb = code_node(e, coder, 0, 0, hist, ENCODE ? ((b_known >> 1) & 1) : 0);
        int lb = code_node(e, coder, 1 + hb, hb, hist, ENCODE ? (b_known & 1) : 0);
        int b  = (hb << 1) | lb;

        if (!ENCODE) (*pout)[i] = i2b(b);

        rcomp = (rcomp >> 2) | ((uint64_t)(b ^ 3) << 62);
        for (auto& m : e.models) {
            m.update(b);
            if (g_cfg.ir) { if (use_fast_ir) m.ir_update_fast(rcomp); else m.ir_update(hist, b); }
        }
        for (auto& s : e.stcms)  s.update(b, hist);
        hist = (hist << 2) | (uint64_t)b;
    }
}

} // namespace

// ── Public API ────────────────────────────────────────────────────────────────
std::vector<uint8_t> dna_encode(
    const std::string&              pg,
    const std::vector<std::string>& /*reads*/,
    const std::vector<uint32_t>&    /*chain_order*/,
    const std::vector<uint32_t>&    /*pg_pos*/)
{
    init_tables();
    load_cfg();
    const size_t N = pg.size();

    std::vector<uint8_t> out(8, 0);
    for (int i = 0; i < 8; ++i) out[i] = (uint8_t)((uint64_t)N >> ((7 - i) * 8));
    if (N == 0) return out;

    autosize_hash_bits(N);
    BinEncoder enc;
    enc.out.reserve(N / 3 + 16);
    EncBitCoder coder{&enc};
    run_core<true>(&pg, nullptr, N, coder);
    enc.flush();

    out.insert(out.end(), enc.out.begin(), enc.out.end());
    return out;
}

std::string dna_decode(const std::vector<uint8_t>& data) {
    init_tables();
    load_cfg();
    if (data.size() < 8) return {};

    uint64_t N = 0;
    for (int i = 0; i < 8; ++i) N = (N << 8) | data[i];
    if (N == 0) return {};

    autosize_hash_bits(N);
    BinDecoder dec(data.data() + 8, data.data() + data.size());
    DecBitCoder coder{&dec};
    std::string pg(N, 'A');
    run_core<false>(nullptr, &pg, N, coder);
    return pg;
}

// ── compute_pg_surprise (idea B) ──────────────────────────────────────────────
// Deterministic low-order causal FCM confidence bucketer. Depends only on `pg`,
// so encoder and decoder derive identical buckets from the identical (decoded)
// pseudogenome — no bytes stored, fully lossless. See dna_coder.h for rationale.
std::vector<uint8_t> compute_pg_surprise(const std::string& pg, int order, int nbuckets) {
    const size_t N = pg.size();
    std::vector<uint8_t> out(N, 0);
    if (N == 0 || nbuckets < 1) return out;
    if (order < 1) order = 1;
    if (order > 16) order = 16;                 // 2 bits/base must fit in 32 bits
    const int    mid  = nbuckets / 2;
    const int    maxb = nbuckets - 1;
    const uint64_t mask = (order * 2 >= 64) ? ~0ULL : ((1ULL << (order * 2)) - 1);

    std::unordered_map<uint64_t, std::array<uint32_t, 4>> tbl;
    tbl.reserve(N / 4 + 16);

    uint64_t ctx = 0;
    int      have = 0;                          // valid bases currently in ctx

    for (size_t i = 0; i < N; ++i) {
        const int b = b2i(pg[i]);               // A/C/G/T → 0..3, else -1
        int bucket;
        if (b < 0 || have < order) {
            bucket = mid;                        // N base or warm-up: neutral
        } else {
            auto it = tbl.find(ctx);
            if (it == tbl.end()) {
                bucket = maxb;                   // never-seen context = max surprise
            } else {
                const auto& c = it->second;
                const uint32_t tot = c[0] + c[1] + c[2] + c[3];
                // Laplace-smoothed probability the model gave the base we saw.
                const double p = (double)(c[b] + 1) / (double)(tot + 4);
                // Higher confidence → lower surprise bucket. Thresholds tuned for
                // nbuckets==4; for other counts they degrade gracefully via clamp.
                if      (p >= 0.85) bucket = 0;
                else if (p >= 0.60) bucket = 1;
                else if (p >= 0.35) bucket = 2;
                else                bucket = 3;
                if (bucket > maxb) bucket = maxb;
            }
        }
        out[i] = (uint8_t)bucket;

        // Causal update: fold the observed base into the model AFTER emitting.
        if (b >= 0) {
            tbl[ctx][b]++;
            ctx  = ((ctx << 2) | (uint64_t)b) & mask;
            have = (have < order) ? have + 1 : order;
        } else {
            ctx = 0; have = 0;                   // reset context across N runs
        }
    }
    return out;
}
