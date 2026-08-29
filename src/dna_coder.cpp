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
#include "container.h"
#include <array>
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <chrono>

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
    //
    // This division has now been attacked TWICE and is not the bottleneck; do
    // not try a third time without new evidence.
    //   1st attempt: reciprocal table indexed by `den`, which spans a large
    //      sparse range (~180 KB) and thrashed the model-table cache.
    //   2nd attempt: reciprocal indexed by the COUNT SUM instead — den is fully
    //      determined by s = n0 + n1, and counts are uint8-capped, so the table
    //      is only 1021 entries (8 KB, L1-resident), avoiding the first
    //      attempt's flaw. recip[s] = floor(2^40/den)+1 was verified BIT-EXACT
    //      against this division over the entire reachable domain (n0,n1 in
    //      [0,510]; 2^32 and 2^34 are NOT exact, so precision matters). It was
    //      still not faster: measured 6.85/6.87 s against a 6.55/6.94 s
    //      baseline — inside run-to-run noise.
    // Reason: the 16 divisions per base (8 models x 2 nodes) are mutually
    // INDEPENDENT, so the out-of-order engine overlaps them and hides the
    // latency, while a table lookup adds load pressure competing with the model
    // tables. Profiling agrees: `code` is ~52% of coder time but that is the
    // mixer, APM and arithmetic coder, not this divide.
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
    bool     do_ir_pf = false;  // emit second prefetch for IR row (set by run_core/ARCS_IR_PF_N)
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
    // Five orders, not eight. Measured on three datasets, the three dropped
    // models (11, 18, 26) cost size AND time -- they were not paying for
    // themselves:
    //
    //   dataset  8 orders                 5 orders                 delta
    //   yeast    13,994,970 / 7.23s       13,978,158 / 5.27s       -16,812 B  -27%
    //   ecoli    46,130,519 / 6.70s       46,122,968 / 5.32s        -7,551 B  -21%
    //   pf       13,087,347 / 4.70s       13,080,256 / 4.17s        -7,091 B  -11%
    //
    // Smaller on every dataset, so this is not a size/speed trade. Two effects
    // compound: each model is a hash probe into a table far larger than L2, so
    // dropping three removes three cache misses per base on BOTH sides; and a
    // near-duplicate of an order already present adds correlated noise to the
    // mixer, which has to spend weight updates learning to discount it. Eight
    // orders was over-fitted -- more models is not more prediction.
    //
    // Compress time is unchanged (19.25 -> 19.31s on yeast): encode is
    // dominated by assembly, not by the coder.
    std::vector<int> orders = {2, 4, 8, 14, 22};
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
    // Upper bound on the table size. Default 24. --fast lowers it (ARCSDNA_HBITS_MAX)
    // to shrink the FCM tables on large/bacterial pseudogenomes where autosize would
    // otherwise pick 2^24 (~64 MB/model). This is a CAP, not a floor: small pgs still
    // size down naturally, so it never inflates a small input.
    // 16 bits: 65,536 entries x 4 B x 2 hashed models = 512 KB, which FITS IN L2.
    // The old default of 24 meant 134 MB of tables and a random probe into them
    // per model per base -- a guaranteed miss to main memory. Measured on all
    // three datasets, capping at 16 is smaller AND faster everywhere:
    //
    //   dataset  cap=24                cap=16                delta
    //   yeast    13,978,885 / 5.20s    13,969,650 / 3.78s    -9,235 B   -27%
    //   ecoli    46,125,657 / 5.54s    46,121,648 / 4.77s    -4,009 B   -14%
    //   pf       13,081,362 / 4.37s    13,077,934 / 2.55s    -3,428 B   -42%
    //
    // Smaller as well as faster, so there is nothing to trade off. A bigger
    // table is not better modelling here: high-order DNA contexts are sparse
    // enough that most entries are seen too few times to learn from, and the
    // collisions a small table introduces act as backoff rather than as loss.
    // This is the same effect that made five context models beat eight.
    //
    // An earlier comment here asserted the cap "cannot help speed: the coder's
    // cost is per-base model updates, not table size". That was wrong, and it
    // was wrong because it was reasoned rather than measured -- the profile puts
    // 51.5% of decode in model SELECT, which is exactly this lookup.
    int cap = 16;
    if (const char* e = getenv("ARCSDNA_HBITS_MAX")) { int v = atoi(e); if (v >= 12 && v <= 24) cap = v; }
    int b = 18; if (b > cap) b = cap;
    while (((size_t)1 << b) < 2 * N && b < cap) ++b;  // ≥2× N entries, capped at 2^cap
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

// ── TSC helper for ARCS_DNA_PROFILE instrumentation ──────────────────────────
// Reads the hardware time-stamp counter. Zero on non-x86 so the profiling path
// compiles but never fires there. rdtsc is not serializing; the minor reorder
// noise across phase boundaries is acceptable for cy/base averages over millions
// of iterations and does not affect the compressed output in any way.
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#  if defined(_MSC_VER)
#    include <intrin.h>
#  else
#    include <x86intrin.h>
#  endif
static inline uint64_t arcs_rdtsc() { return (uint64_t)__rdtsc(); }
#else
static inline uint64_t arcs_rdtsc() { return 0; }
#endif

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
void run_core(const std::string* pin, std::string* pout, size_t N, BitCoder& coder,
              const std::string& seed = "") {
    Engine e; e.build();
    bool use_fast_ir = true;
    for (auto& m : e.models) if (m.order > 31) use_fast_ir = false;
    uint64_t hist  = 0;
    uint64_t rcomp = 0;

    // Context seed: warm hist/rcomp shift-registers + freq tables from preceding block tail.
    // Both encoder and decoder run identically since seed = deterministic preceding bytes.
    // Must replicate the exact main-loop update order: roll rcomp first, then ir_update_fast,
    // then hist roll — matching the phase-3 sequence in the coded path.
    if (!seed.empty()) {
        for (char c : seed) {
            int b = b2i(c);
            for (auto& m : e.models) { m.select(hist); m.update(b); }
            for (auto& s : e.stcms)  { s.select();    s.update(b, hist); }
            rcomp = (rcomp >> 2) | ((uint64_t)(b ^ 3) << 62);  // roll before ir_update_fast
            if (g_cfg.ir) {
                if (use_fast_ir) for (auto& m : e.models) m.ir_update_fast(rcomp);
                else             for (auto& m : e.models) m.ir_update(hist, b);
            }
            hist  = (hist  << 2) | (uint64_t)b;
        }
    }

    // ── ARCS_IR_PF_N: selective IR-row prefetch ───────────────────────────────
    // Mark the top-N hashed models (highest order first) to emit an IR-row
    // prefetch in the select phase alongside the existing forward prefetch.
    //
    // For ENCODE the prefetch uses the EXACT rcomp_new = (rcomp>>2)|(comp(b)<<62)
    // because b is known at the top of the loop; index equality is 100%.
    // For DECODE b is unknown until code_node completes so we use the pre-roll
    // rcomp as the best approximation (1-base lag; index accuracy < 100%).
    //
    // ARCS_IR_PF_N controls how many models get the prefetch.
    // Default 2 (top-2 hashed orders: 26,22) — LFB budget 8+2=10/12, measured
    // −14.3% wall / −18.4% update on DS1; N=4 hits LFB limit and DEGRADES.
    // Set ARCS_IR_PF_N=0 to disable; 4+ not recommended (LFB saturation).
    {
        int n_ir = 2;
        if (const char* s = getenv("ARCS_IR_PF_N")) n_ir = atoi(s);
        if (n_ir > 0 && use_fast_ir && g_cfg.ir) {
            // Collect hashed-model indices sorted by order descending.
            int hmi[8]; int nhm = 0;
            for (int mi = 0; mi < (int)e.models.size(); ++mi)
                if (!e.models[mi].dense) hmi[nhm++] = mi;
            // Simple insertion sort (≤8 elements).
            for (int a = 1; a < nhm; ++a) {
                int key = hmi[a], kb = e.models[key].order;
                int b2 = a - 1;
                while (b2 >= 0 && e.models[hmi[b2]].order < kb) { hmi[b2+1] = hmi[b2]; --b2; }
                hmi[b2+1] = key;
            }
            int mark = n_ir < nhm ? n_ir : nhm;
            for (int k = 0; k < mark; ++k) e.models[hmi[k]].do_ir_pf = true;
        }
    }
    bool any_ir_pf = false;
    for (auto& m : e.models) if (m.do_ir_pf) { any_ir_pf = true; break; }
    // ARCS_SPEC_PF_N: how many of the widest models get speculative
    // next-position prefetches (default 2, 0 disables). See the issue site.
    int spec_pf_n = 4;   // measured optimum: 6.73s -> 5.07s at 4; 6 and 8 regress (bandwidth)
    if (const char* sp = getenv("ARCS_SPEC_PF_N")) { int v = atoi(sp); if (v >= 0 && v <= 8) spec_pf_n = v; }

    // ── ARCS_IR_VERIFY: index-equality audit ──────────────────────────────────
    // Compares the table index the IR prefetch would target (computed at select
    // time from rcomp_exact) against the index ir_update_fast() actually uses
    // (computed after rcomp is rolled). For ENCODE with exact rcomp_new, equality
    // must be 100%. Any deviation exposes a bug in the address computation.
    // Runs over all N bases (or VERIFY_LIMIT if env var specifies fewer).
    const bool do_verify = (getenv("ARCS_IR_VERIFY") != nullptr);
    uint64_t ver_match = 0, ver_total = 0;
    // Capture prefetch indices across select→update for up to 8 models.
    size_t pf_idx_ver[8] = {};

    // ── Per-base TSC profiling (ARCS_DNA_PROFILE=1) ───────────────────────────
    const bool do_prof = (getenv("ARCS_DNA_PROFILE") != nullptr);
    uint64_t cyc_sel = 0, cyc_code = 0, cyc_upd = 0, tsc0 = 0, prof_tsc_start = 0;
    using HRC = std::chrono::high_resolution_clock;
    HRC::time_point prof_wall_start{};
    if (do_prof) { prof_wall_start = HRC::now(); prof_tsc_start = arcs_rdtsc(); }

    for (size_t i = 0; i < N; ++i) {
        int b_known = ENCODE ? b2i((*pin)[i]) : 0;

        // Exact rcomp_new for ENCODE (b is known); 1-step lag for DECODE.
        // This is the value ir_update_fast() will receive this iteration.
        const uint64_t rcomp_ir = ENCODE
            ? ((rcomp >> 2) | ((uint64_t)(b_known ^ 3) << 62))
            : rcomp;

        // ── Phase 1: select + IR-row prefetch ────────────────────────────────
        if (do_prof) tsc0 = arcs_rdtsc();
        for (auto& m : e.models) m.select(hist);
        for (auto& s : e.stcms)  s.select();
        // Selective IR prefetch — issued immediately after forward prefetches so
        // both sets of fill-buffer requests are live before code_node starts.
        if (any_ir_pf) {
            for (auto& m : e.models) {
                if (m.do_ir_pf) {
                    size_t ir_idx = m.index(rcomp_ir >> m.ir_ctx_shift);
                    __builtin_prefetch(&m.tbl[ir_idx], 1, 1);
                    // Capture for verification in same pass (zero extra index call).
                    if (do_verify) {
                        size_t mi = (size_t)(&m - e.models.data());
                        pf_idx_ver[mi] = ir_idx;
                    }
                }
            }
        }
        // ── Speculative prefetch for the NEXT position ───────────────────────
        // Decode is a dependency chain: decode a base, fold it into hist,
        // compute the table index, wait for the load, decode the next. The
        // existing select() prefetch is issued immediately before the value is
        // used, so it has no time to land and every iteration pays a full miss.
        // pg_decode is 6.70s of an 8.88s decompress at 2.44 Mchar/s, which is
        // memory latency, not arithmetic.
        //
        // The next context is unknown -- but only four values are possible. So
        // issue all four candidate indices now, before code_node runs; whichever
        // base turns out to be right, its line is already in flight and the
        // other three are harmless. Prefetch is a hint: this cannot change a
        // single output byte, only when the data arrives.
        //
        // Restricted to the widest models by ARCS_SPEC_PF_N (default 2): the
        // low-order tables are small enough to stay cached, so speculating on
        // them would spend bandwidth for nothing.
        if (spec_pf_n > 0) {
            int done = 0;
            for (size_t mi = e.models.size(); mi-- > 0 && done < spec_pf_n; ) {
                auto& m = e.models[mi];
                if (m.dense) continue;   // small tables stay cached; no gain
                const uint64_t base_hist = hist << 2;
                __builtin_prefetch(&m.tbl[m.index(base_hist | 0ULL)], 1, 1);
                __builtin_prefetch(&m.tbl[m.index(base_hist | 1ULL)], 1, 1);
                __builtin_prefetch(&m.tbl[m.index(base_hist | 2ULL)], 1, 1);
                __builtin_prefetch(&m.tbl[m.index(base_hist | 3ULL)], 1, 1);
                ++done;
            }
        }
        if (do_prof) cyc_sel += arcs_rdtsc() - tsc0;

        // ── Phase 2: code_node × 2 ───────────────────────────────────────────
        if (do_prof) tsc0 = arcs_rdtsc();
        int hb = code_node(e, coder, 0, 0, hist, ENCODE ? ((b_known >> 1) & 1) : 0);
        int lb = code_node(e, coder, 1 + hb, hb, hist, ENCODE ? (b_known & 1) : 0);
        int b  = (hb << 1) | lb;
        if (do_prof) cyc_code += arcs_rdtsc() - tsc0;

        if (!ENCODE) (*pout)[i] = i2b(b);

        // ── Phase 3: update ───────────────────────────────────────────────────
        rcomp = (rcomp >> 2) | ((uint64_t)(b ^ 3) << 62);

        // Index-equality verification: rcomp is now the value ir_update_fast()
        // will receive; compare actual index against the prefetch index captured
        // in the select phase. For ENCODE with exact rcomp_ir this must be 100%.
        if (do_verify && any_ir_pf) {
            for (size_t mi = 0; mi < e.models.size(); ++mi) {
                if (e.models[mi].do_ir_pf) {
                    size_t actual = e.models[mi].index(rcomp >> e.models[mi].ir_ctx_shift);
                    if (actual == pf_idx_ver[mi]) ++ver_match;
                    ++ver_total;
                }
            }
        }

        if (do_prof) tsc0 = arcs_rdtsc();
        for (auto& m : e.models) {
            m.update(b);
            if (g_cfg.ir) { if (use_fast_ir) m.ir_update_fast(rcomp); else m.ir_update(hist, b); }
        }
        for (auto& s : e.stcms)  s.update(b, hist);
        if (do_prof) cyc_upd += arcs_rdtsc() - tsc0;

        hist = (hist << 2) | (uint64_t)b;
    }

    // ── Verification report ───────────────────────────────────────────────────
    if (do_verify && ver_total > 0) {
        int n_marked = 0; for (auto& m : e.models) if (m.do_ir_pf) ++n_marked;
        fprintf(stderr,
            "[IR_VERIFY] mode=%-6s  N=%zu  ir_pf_models=%d\n"
            "[IR_VERIFY] match=%zu / %zu  hit_rate=%.4f%%  %s\n",
            ENCODE ? "ENCODE" : "DECODE",
            N, n_marked,
            ver_match, ver_total,
            100.0 * ver_match / ver_total,
            (ver_match == ver_total) ? "EXACT (100%)" : "PARTIAL — address mismatch detected");
        fflush(stderr);
    }

    // ── Profiling report ──────────────────────────────────────────────────────
    if (do_prof && N > 0) {
        uint64_t tsc_total = arcs_rdtsc() - prof_tsc_start;
        double wall_s = std::chrono::duration<double>(HRC::now() - prof_wall_start).count();
        double ghz     = (wall_s > 1e-6) ? (double)tsc_total / (wall_s * 1e9) : 3.0;
        uint64_t cyc_other = (tsc_total > cyc_sel + cyc_code + cyc_upd)
                             ? tsc_total - cyc_sel - cyc_code - cyc_upd : 0;
        double inv = 1.0 / (double)N;
        size_t hashed_models = 0;
        for (auto& m : e.models) if (!m.dense) ++hashed_models;
        size_t dense_models  = e.models.size() - hashed_models;
        size_t tbl_bytes = hashed_models * ((size_t)1 << Model::hash_bits) * 4
                         + dense_models  * ((size_t)1 << (2 * 8)) * 4;
        int n_ir_pf = 0; for (auto& m : e.models) if (m.do_ir_pf) ++n_ir_pf;
        fprintf(stderr,
            "[DNA_PROFILE] ── run_core<%s>  N=%-9zu  models=%zu  ir_pf=%d\n"
            "[DNA_PROFILE]    hash_bits=%-2d  hashed_tbl=%zu entries × 4 B × %zu models\n"
            "[DNA_PROFILE]    approx FCM working set ~%.1f MB\n"
            "[DNA_PROFILE]    wall=%.3f s  tsc_total=%.3e cy  ghz_est=%.2f\n"
            "[DNA_PROFILE]\n"
            "[DNA_PROFILE]    cy/base  select=%7.1f  code=%7.1f  update=%7.1f  other=%7.1f  total=%7.1f\n"
            "[DNA_PROFILE]    %%total   select=%5.1f%%  code=%5.1f%%  update=%5.1f%%  other=%5.1f%%\n",
            ENCODE ? "ENCODE" : "DECODE",
            N, e.models.size(), n_ir_pf,
            Model::hash_bits, (size_t)1 << Model::hash_bits, hashed_models,
            (double)tbl_bytes / (1024.0 * 1024.0),
            wall_s, (double)tsc_total, ghz,
            cyc_sel  * inv, cyc_code * inv, cyc_upd  * inv, cyc_other * inv, tsc_total * inv,
            100.0 * cyc_sel  / tsc_total, 100.0 * cyc_code / tsc_total,
            100.0 * cyc_upd  / tsc_total, 100.0 * cyc_other / tsc_total);
        fflush(stderr);
    }
}

} // namespace

// ── Public API ────────────────────────────────────────────────────────────────
std::vector<uint8_t> dna_encode(
    const std::string&              pg,
    const std::vector<std::string>& /*reads*/,
    const std::vector<uint32_t>&    /*chain_order*/,
    const std::vector<uint32_t>&    /*pg_pos*/,
    const std::string&              seed)
{
    init_tables();
    load_cfg();
    const size_t N = pg.size();

    std::vector<uint8_t> out(8, 0);
    for (int i = 0; i < 8; ++i) out[i] = (uint8_t)((uint64_t)N >> ((7 - i) * 8));
    if (N == 0) return out;

    // Persist the FCM order list into the stream so decode reconstructs the exact
    // model set regardless of any ARCSDNA_ORDERS env in the (separate) decode
    // process. This makes the DNA stream fully self-describing and fixes a latent
    // losslessness hazard where an env override on only one side desynced the coder.
    // Order values are ≤ 63 and count ≤ 63, so each fits in one byte.
    autosize_hash_bits(N);   // fixes Model::hash_bits (from env override or auto-size)
    {
        size_t no = g_cfg.orders.size();
        if (no > 255) no = 255;
        out.push_back((uint8_t)no);
        for (size_t k = 0; k < no; ++k) {
            int o = g_cfg.orders[k]; if (o < 1) o = 1; if (o > 255) o = 255;
            out.push_back((uint8_t)o);
        }
        // Persist the hashed-model table size too, so a smaller ARCSDNA_HBITS cap (used
        // by --fast to shrink FCM RAM) is reproduced exactly at decode without the env.
        out.push_back((uint8_t)Model::hash_bits);
    }

    BinEncoder enc;
    enc.out.reserve(N / 3 + 16);
    EncBitCoder coder{&enc};
    run_core<true>(&pg, nullptr, N, coder, seed);
    enc.flush();

    out.insert(out.end(), enc.out.begin(), enc.out.end());
    return out;
}

std::string dna_decode(const std::vector<uint8_t>& data, const std::string& seed) {
    init_tables();
    load_cfg();
    if (data.size() < 8) return {};

    uint64_t N = 0;
    for (int i = 0; i < 8; ++i) N = (N << 8) | data[i];
    if (N == 0) return {};

    // Read the persisted FCM order list (written by dna_encode after the N header)
    // and make it authoritative — this is what the stream was actually coded with,
    // overriding any ARCSDNA_ORDERS env in this process so decode always matches.
    size_t off = 8;
    if (off < data.size()) {
        size_t no = data[off++];
        std::vector<int> ords; ords.reserve(no);
        for (size_t k = 0; k < no && off < data.size(); ++k) ords.push_back((int)data[off++]);
        if (!ords.empty()) g_cfg.orders = ords;
        // Read the persisted hashed-model table size and make it authoritative (matches
        // the encoder's actual size, incl. any --fast HBITS cap) so autosize can't diverge.
        if (off < data.size()) {
            int hb = data[off++];
            if (hb >= 1 && hb <= 30) { Model::hash_bits = hb; g_hbits_from_env = true; }
        }
    }

    autosize_hash_bits(N);   // early-returns: g_hbits_from_env set from the persisted value
    BinDecoder dec(data.data() + off, data.data() + data.size());
    DecBitCoder coder{&dec};
    std::string pg(N, 'A');
    run_core<false>(nullptr, &pg, N, coder, seed);
    return pg;
}

// ── VLE + LZMA pseudogenome backend ──────────────────────────────────────────
// 2-bit pack (A=0,C=1,G=2,T=3, 4 bases/byte, MSB first) + LZMA-9.
// Each block is independently decodable — no warm-up, no cross-block state.
// Decode is ~10-50ms per block vs ~1-2s for adaptive FCM (LZMA is streaming).
// Ratio cost vs FCM: ~5-20% larger on human pg (FCM models long-range repeats
// via high-order context; LZMA's sliding window is shorter).

// Plain LZMA-9 on ASCII pg text — no 2-bit packing.
// Zero regression on datasets where LZMA already beats FCM (bacterial/DS1).
// Some ratio cost on human WGS where FCM exploits long-range repeats.
std::vector<uint8_t> vle_encode_pg(const std::string& pg) {
    size_t N = pg.size();
    std::vector<uint8_t> out(8, 0);
    for (int i = 0; i < 8; i++) out[i] = (uint8_t)((uint64_t)N >> ((7 - i) * 8));
    if (N == 0) return out;
    std::vector<uint8_t> bytes(pg.begin(), pg.end());
    auto comp = arcs_compress(bytes, 9);
    out.insert(out.end(), comp.begin(), comp.end());
    return out;
}

std::string vle_decode_pg(const std::vector<uint8_t>& data) {
    if (data.size() < 8) return {};
    uint64_t N = 0;
    for (int i = 0; i < 8; i++) N = (N << 8) | data[i];
    if (N == 0) return {};
    auto comp = std::vector<uint8_t>(data.begin() + 8, data.end());
    auto raw  = arcs_decompress(comp.data(), comp.size());
    if (raw.empty()) return {};
    return std::string(raw.begin(), raw.end());
}

// ── 2-bit + LZMA pg codec (format 0x08) ──────────────────────────────────────
// Z4 encoding: A=0, C=1, T=2, G=3 (complement = +2 mod 4, complement = XOR high bit).
// Pack 4 bases per byte MSB-first, LZMA-9 compress. Decode ~10ms vs 26s FCM.
// M2 result: +2.2% vs FCM on DS2 (trial24). Non-ACTG chars → fall back to ASCII LZMA.

static const char S_I2B[4] = {'A','C','T','G'};
static int S_B2I[256];
static bool s_b2i_ready = false;
static void s_init_b2i() {
    if (s_b2i_ready) return;
    for (int& v : S_B2I) v = 0;
    S_B2I[(uint8_t)'A'] = S_B2I[(uint8_t)'a'] = 0;
    S_B2I[(uint8_t)'C'] = S_B2I[(uint8_t)'c'] = 1;
    S_B2I[(uint8_t)'T'] = S_B2I[(uint8_t)'t'] = 2;
    S_B2I[(uint8_t)'G'] = S_B2I[(uint8_t)'g'] = 3;
    s_b2i_ready = true;
}

static bool s_pg_is_pure_actg(const std::string& pg) {
    for (unsigned char c : pg)
        if (c != 'A' && c != 'C' && c != 'T' && c != 'G') return false;
    return true;
}

std::vector<uint8_t> pg_encode_2bit(const std::string& pg) {
    s_init_b2i();
    if (!s_pg_is_pure_actg(pg)) {
        // Non-ACTG chars present: fall back to LZMA on ASCII (same ratio, fully lossless).
        return vle_encode_pg(pg);
    }
    const size_t N = pg.size();
    // 8-byte big-endian N header
    std::vector<uint8_t> out(8);
    for (int k = 7; k >= 0; --k) out[7-k] = (uint8_t)(N >> (k * 8));
    // Pack: 4 bases/byte, MSB-first (base i → bits [7-6, 5-4, 3-2, 1-0])
    std::vector<uint8_t> packed((N + 3) / 4, 0);
    for (size_t i = 0; i < N; ++i)
        packed[i >> 2] |= (uint8_t)(S_B2I[(uint8_t)pg[i]] << (6 - (i & 3) * 2));
    // zstd-6: ~1 GB/s decompress vs LZMA's ~50 MB/s → 20× faster; ratio ~= LZMA-3.
    // arcs_decompress auto-detects zstd magic (at byte 8 of the zstd_compress output).
    auto comp = zstd_compress(packed, 6);
    out.insert(out.end(), comp.begin(), comp.end());
    return out;
}

std::string pg_decode_2bit(const std::vector<uint8_t>& data) {
    s_init_b2i();
    if (data.size() < 8) return {};
    uint64_t N = 0;
    for (int i = 0; i < 8; ++i) N = (N << 8) | data[i];
    if (N == 0) return {};
    std::vector<uint8_t> comp(data.begin() + 8, data.end());
    auto packed = arcs_decompress(comp.data(), comp.size());
    if (packed.empty()) return {};
    std::string out(N, 'A');
    for (size_t i = 0; i < N; ++i)
        out[i] = S_I2B[(packed[i >> 2] >> (6 - (i & 3) * 2)) & 3];
    return out;
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
