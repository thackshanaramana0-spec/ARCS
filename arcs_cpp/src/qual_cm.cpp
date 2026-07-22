// ── Adaptive Q-ary quality coder (own generative-model context) ───────────────
// See qual_cm.h and trial20_arcs_dna_crossstream/QUALITY_SOTA.md.
//
// DESIGN RATIONALE (not a transcription of fqzcomp):
// A Phred value is the sequencer's own estimate of per-base error probability. It
// is a smooth latent trajectory + small noise, driven by: (a) cycle/position
// drift, (b) autocorrelation with the previous value(s), (c) reversion to a local
// quality CEILING (the signal dips then returns to a recent best), (d) local
// volatility. We model those latent variables directly:
//   - q1, q2                : autocorrelation (order-2)
//   - ceil (decaying max)   : the local ceiling the signal reverts to. This is a
//                             DECAYING running-max — a general reversion estimator,
//                             of which fqzcomp's max(q2,q3) is only a crude 2-tap.
//   - vol (EWMA of |Δq|)    : how bumpy the recent run is (our own formulation)
//   - position bin          : cycle drift
//   - is_dev (mismatch flag): reference-free cross-stream signal (see chain-pg)
// Each quality symbol is coded WHOLE against a per-context ADAPTIVE frequency
// table with a classic carryless range coder — the model is never transmitted, so
// rich context is free (the lesson the static rANS model could not exploit).
// Integer-only → encode and decode are bit-exact and portable → lossless.
#include "qual_cm.h"
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <thread>

namespace {

// ── Classic carryless range coder (Subbotin form) ─────────────────────────────
constexpr uint32_t RC_TOP = 1u << 24;
constexpr uint32_t RC_BOT = 1u << 16;

struct RangeEnc {
    uint32_t low = 0, range = 0xFFFFFFFFu;
    std::vector<uint8_t> out;
    inline void encode(uint32_t cum, uint32_t freq, uint32_t tot) {
        range /= tot;
        low   += cum * range;
        range *= freq;
        while ((low ^ (low + range)) < RC_TOP ||
               (range < RC_BOT && ((range = (0u - low) & (RC_BOT - 1)), true))) {
            out.push_back((uint8_t)(low >> 24));
            low <<= 8; range <<= 8;
        }
    }
    inline void flush() { for (int i = 0; i < 4; ++i) { out.push_back((uint8_t)(low >> 24)); low <<= 8; } }
};

struct RangeDec {
    uint32_t low = 0, range = 0xFFFFFFFFu, code = 0;
    const uint8_t* ptr; const uint8_t* end;
    RangeDec(const uint8_t* p, const uint8_t* e) : ptr(p), end(e) {
        for (int i = 0; i < 4; ++i) code = (code << 8) | (ptr < end ? *ptr++ : 0u);
    }
    inline uint32_t getfreq(uint32_t tot) { range /= tot; return (code - low) / range; }
    inline void update(uint32_t cum, uint32_t freq) {
        low   += cum * range;
        range *= freq;
        while ((low ^ (low + range)) < RC_TOP ||
               (range < RC_BOT && ((range = (0u - low) & (RC_BOT - 1)), true))) {
            code = (code << 8) | (ptr < end ? *ptr++ : 0u);
            low <<= 8; range <<= 8;
        }
    }
};

// ── Adaptive per-context frequency table ──────────────────────────────────────
// Small dense array (quality alphabet ≤ 48). Increment on observe, halve on
// overflow — a standard adaptive model. Total kept < RC_BOT so range/tot is safe.
constexpr int   QA        = 94;      // max alphabet slots (full legal Phred 0..93)
constexpr int   FREQ_MAX  = 60000;   // rescale threshold (< 65536)
// Per-run quality-CM config. Passed EXPLICITLY into encode_block/decode_block (never
// via mutable globals) so parallel quality blocks AND concurrent chunk-decoders are
// race-free and byte-exact — a shared mutable USE_SEQ silently corrupted block workers
// (they saw the default, dropping the sequence-context model → +5 MB on DS7 quality).
struct QcmCfg { int freq_inc = 4; int use_seq = 0; int qcm_k = 20; };
static QcmCfg qcm_cfg_from_env() {
    QcmCfg c;
    if (const char* e = getenv("ARCS_QCM_INC")) { int v = atoi(e); if (v>=1 && v<=200) c.freq_inc = v; }
    if (const char* e = getenv("ARCS_QCM_K"))   { int v = atoi(e); if (v>=1 && v<=255) c.qcm_k   = v; }
    return c;
}
static inline int b2i5(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2;
                 case 'T': return 3; default: return 4; }  // N/other = 4
}
// Local 3-mer sequence context {prev,cur,next base} at position j → 0..124.
static inline int seq3(const std::string* s, int j, int L) {
    if (!s) return 0;
    int p = (j > 0)     ? b2i5((*s)[j - 1]) : 4;
    int c = (j < (int)s->size()) ? b2i5((*s)[j]) : 4;
    int n = (j + 1 < L) ? b2i5((*s)[j + 1]) : 4;
    return (p * 5 + c) * 5 + n;
}

struct Freq {
    uint16_t f[QA];
    uint32_t tot;
    uint32_t obs;                // real observation count (drives interpolation λ)
    Freq() { for (int i = 0; i < QA; ++i) f[i] = 1; tot = QA; obs = 0; }
    inline void update(int sym, int freq_inc) {
        f[sym] = (uint16_t)(f[sym] + freq_inc); tot += freq_inc; ++obs;
        if (tot >= FREQ_MAX) {
            uint32_t t = 0;
            for (int i = 0; i < QA; ++i) { f[i] = (uint16_t)((f[i] >> 1) | 1); t += f[i]; }
            tot = t;
        }
    }
};

// ── Flat open-addressing context table (drop-in for unordered_map<u64,Freq>) ──
// The CM math depends only on which Freq a key maps to, not on the container, so
// this is BIT-IDENTICAL to the unordered_map version (same Freq values, same
// update order) — zero ratio change, faster. Keys live in a contiguous 8-byte
// array (cache-friendly linear probing); the 196-byte Freq structs live in a
// separate pool allocated only for contexts that actually occur, so this also uses
// LESS memory than the node-based map (no per-node pointers, no empty 196 B slots).
struct FreqMap {
    static constexpr uint64_t EMPTY = ~0ULL;   // real keys are ≤ ~72M, never ~0ULL
    static constexpr size_t   CHUNK = 8192;    // Freqs per pool chunk
    std::vector<uint64_t> keys;
    std::vector<uint32_t> slot_idx;                 // slot → pool index
    std::vector<std::vector<Freq>> chunks;          // chunked pool: exactly count Freqs, no 2× slack
    size_t mask = 0, count = 0;

    explicit FreqMap(unsigned bits = 15) { init(bits); }
    void init(unsigned bits) {
        size_t cap = (size_t)1 << bits;
        keys.assign(cap, EMPTY);
        slot_idx.assign(cap, 0u);
        chunks.clear();
        mask = cap - 1; count = 0;
    }
    static inline size_t mix(uint64_t k) { return (size_t)(((k + 1) * 0x9E3779B97F4A7C15ULL) >> 32); }
    // pool index → Freq& (chunk buffers are heap-stable across chunks-vector growth).
    inline Freq& poolref(uint32_t pi) { return chunks[pi / CHUNK][pi % CHUNK]; }
    inline Freq& get(uint64_t k) {
        size_t i = mix(k) & mask;
        while (keys[i] != EMPTY) {
            if (keys[i] == k) return poolref(slot_idx[i]);
            i = (i + 1) & mask;
        }
        keys[i] = k;
        uint32_t pi = (uint32_t)count;
        if (pi % CHUNK == 0) chunks.emplace_back(CHUNK);   // never reallocs existing Freqs
        slot_idx[i] = pi;
        ++count;
        if (count * 10 >= (mask + 1) * 7) grow();           // grow rehashes keys only; pi stays valid
        return poolref(pi);
    }
    void grow() {
        size_t newcap = (mask + 1) << 1, nmask = newcap - 1;
        std::vector<uint64_t> nk(newcap, EMPTY);
        std::vector<uint32_t> ni(newcap, 0u);
        for (size_t s = 0; s <= mask; ++s) {
            if (keys[s] == EMPTY) continue;
            size_t i = mix(keys[s]) & nmask;
            while (nk[i] != EMPTY) i = (i + 1) & nmask;
            nk[i] = keys[s]; ni[i] = slot_idx[s];
        }
        keys.swap(nk); slot_idx.swap(ni); mask = nmask;
    }
};

// ── Interpolated coding distribution (PPM-style backoff) ──────────────────────
// Build an integer coding table cf[0..qmax] ∝ λ·P_child + (1−λ)·P_parent, with
// λ = obs_c/(obs_c+K). Derivation (common denominator (obs_c+K)·tot_c·tot_p):
//   weight(q) = obs_c·tot_p·f_child[q] + K·tot_c·f_parent[q]
// Rescaled to total ≈ TARGET so the range coder stays in bounds. Deterministic
// integer math → encoder and decoder build identical cf → lossless.
constexpr uint32_t CF_TARGET = 1u << 14;

// Exact unsigned 64-bit division by an invariant divisor (libdivide algorithm,
// Granlund-Montgomery). Precompute the reciprocal ONCE for a divisor, then each
// division is a mulhi + shift instead of a hardware idiv (~30→~5 cycles). Produces
// the EXACT floor(n/d) — bit-identical to `/` — so build_cf output is unchanged.
struct U64Div {
    uint64_t magic; uint8_t more;   // more: low 6 bits = shift, bit 6 = ADD marker
    explicit U64Div(uint64_t d) {
        uint32_t l = 63 - (uint32_t)__builtin_clzll(d);
        if ((d & (d - 1)) == 0) { magic = 0; more = (uint8_t)l; return; }   // power of two
        __uint128_t n = ((__uint128_t)1) << (64 + l);
        uint64_t pm  = (uint64_t)(n / d);
        uint64_t rem = (uint64_t)(n - (__uint128_t)pm * d);
        uint64_t e   = d - rem;
        if (e < ((uint64_t)1 << l)) { more = (uint8_t)l; }
        else {
            pm += pm; uint64_t tr = rem + rem;
            if (tr >= d || tr < rem) pm += 1;
            more = (uint8_t)(l | 0x40);
        }
        magic = pm + 1;
    }
    inline uint64_t operator()(uint64_t n) const {
        if (magic == 0) return n >> more;
        uint64_t q = (uint64_t)(((__uint128_t)magic * n) >> 64);
        if (more & 0x40) { uint64_t t = ((n - q) >> 1) + q; return t >> (more & 0x3F); }
        return q >> more;
    }
};

inline uint32_t build_cf(const Freq& child, const Freq& parent, int qmax,
                         uint32_t* cf, int qcm_k) {
    const uint64_t nc = child.obs;
    const uint64_t tp = parent.tot, tc = child.tot;
    const uint64_t A = nc * tp;                 // hoist constant products out of loop
    const uint64_t B = (uint64_t)qcm_k * tc;
    uint64_t w[QA]; uint64_t sum = 0;
    for (int q = 0; q <= qmax; ++q) {
        uint64_t wq = A * child.f[q] + B * parent.f[q];
        w[q] = wq; sum += wq;
    }
    if (sum == 0) { for (int q = 0; q <= qmax; ++q) cf[q] = 1; return (uint32_t)(qmax + 1); }
    U64Div dv(sum);                             // one reciprocal for all qmax+1 divisions
    uint32_t tot = 0;
    for (int q = 0; q <= qmax; ++q) {
        uint32_t v = (uint32_t)dv(w[q] * CF_TARGET);
        if (v == 0) v = 1;
        cf[q] = v; tot += v;
    }
    return tot;
}

// ── Context (our generative-model latent variables) ───────────────────────────
struct CtxState {
    int q1 = 0, q2 = 0, ceil = 0, vol = 0, run = 0;
    inline void advance(int q) {
        ceil = std::max(q, ceil - 1);            // decaying running-max
        vol  = (vol * 3 + 8 * std::abs(q - q1)) >> 2;  // EWMA of |Δq|
        run  = (q == q1) ? (run < 7 ? run + 1 : 7) : 0; // consecutive-repeat length of q1
        q2 = q1; q1 = q;
    }
    inline void reset() { q1 = q2 = ceil = vol = run = 0; }
};

inline int posbin(int j, int L) { int pb = (int)((int64_t)j * 8 / (L > 0 ? L : 1)); return pb > 7 ? 7 : pb; }

inline uint64_t ctx_key(const CtxState& s, int j, int L, int is_dev, int qmax) {
    int pb   = posbin(j, L);
    int ceilb= std::min(s.ceil / 3, 15);
    int volb = std::min(s.vol / 8, 3);
    int q1   = std::min(s.q1, qmax);
    int q2   = std::min(s.q2, qmax);
    uint64_t k = (uint32_t)q1;
    k = k * 94u + (uint32_t)q2;
    k = k * 16u + (uint32_t)ceilb;
    k = k * 4u  + (uint32_t)volb;
    k = k * 8u  + (uint32_t)pb;
    k = k * 8u  + (uint32_t)s.run;
    k = k * 2u  + (uint32_t)(is_dev ? 1 : 0);
    return k;
}
// Low-order PARENT context {q1, q2, pos} — dense, warms fast; used to seed a new
// rich context so it does not start cold (uniform). Both sides create rich
// contexts in the same symbol order with the same parent state → deterministic.
inline uint64_t lo_key(const CtxState& s, int j, int L, int qmax) {
    uint64_t k = (uint32_t)std::min(s.q1, qmax);
    k = k * 48u + (uint32_t)std::min(s.q2, qmax);
    k = k * 8u  + (uint32_t)posbin(j, L);
    k = k * 8u  + (uint32_t)s.run;   // run-length: dominant lever on binned quality
    return k;
}

// Compute the (parent, coded) context keys for one position; the coder codes
// each symbol against λ·P_coded + (1−λ)·P_parent (interpolation).
//   USE_SEQ : parent = dense {q1,q2,pos}; coded = parent ⊗ local 3-mer sequence.
//             This is the held-out-measured winner (2.204 bpq): the sequence
//             context provides the signal, and per-symbol backoff to the warm
//             {q1,q2,pos} parent removes its cold-start cost.
//   else    : parent = {q1,q2,pos}; coded = full quality context (q1,q2,ceil,vol,
//             pos,is_dev) — the quality-history model used when no sequence given.
// qmb: mate-quality bin (0=<10, 1=10-20, 2=20-30, 3=>=30, 4=no-mate/R1/SE).
// Added ONLY to coded_k so the dense parent (lo) stays warm without the extra
// dimension. The interpolation backoff absorbs cold-start for new qmb contexts.
inline void keys(const CtxState& s, int j, int L, int is_dev, int qmax,
                 const std::string* seq, uint64_t& parent_k, uint64_t& coded_k,
                 int use_seq, int qmb = 4) {
    const uint32_t q1c = (uint32_t)std::min(s.q1, qmax);
    const uint32_t q2c = (uint32_t)std::min(s.q2, qmax);
    const uint32_t pb  = (uint32_t)posbin(j, L);
    uint64_t lk = q1c;
    lk = lk * 48u + q2c;
    lk = lk * 8u  + pb;
    lk = lk * 8u  + (uint32_t)s.run;
    parent_k = lk;                               // parent unchanged — stays dense
    if (use_seq) {
        coded_k = (lk * 125u + (uint32_t)seq3(seq, j, L)) * 5u + (uint32_t)qmb;
    } else {
        const uint32_t ceilb = (uint32_t)std::min(s.ceil / 3, 15);
        const uint32_t volb  = (uint32_t)std::min(s.vol / 8, 3);
        uint64_t k = q1c;
        k = k * 94u + q2c;
        k = k * 16u + ceilb;
        k = k * 4u  + volb;
        k = k * 8u  + pb;
        k = k * 8u  + (uint32_t)s.run;
        k = k * 2u  + (uint32_t)(is_dev ? 1 : 0);
        coded_k = k * 5u + (uint32_t)qmb;
    }
}

// ── One block = an independent CM model + range coder over order[begin,end) ────
// Blocks share nothing (own tbl/lo maps, own coder, CtxState reset per read), so
// they encode/decode concurrently and are bit-exact independent of thread count. The
// config (freq_inc/use_seq/qcm_k) is passed BY VALUE (not shared globals) so it is
// race-free across blocks and across concurrent chunk-decoders.
static std::vector<uint8_t> encode_block(
    const std::vector<std::vector<uint8_t>>& rq,
    const std::vector<uint32_t>& order, size_t begin, size_t end,
    const std::vector<std::vector<bool>>& dev_sets, int L,
    const std::vector<std::string>* seqs, int qmax, QcmCfg cfg, bool is_pe) {
    FreqMap tbl(15), lo(12);
    RangeEnc enc;
    uint32_t cf[QA];
    // Track which R1 indices were processed before their R2 in this block.
    // Encoder and decoder make the identical decision (same order, same set) → lossless.
    std::unordered_set<uint32_t> seen;
    for (size_t oi = begin; oi < end; ++oi) {
        uint32_t idx = order[oi];
        if (idx >= rq.size()) { seen.insert(idx); continue; }
        const auto& r = rq[idx];
        const bool has_dev = idx < dev_sets.size();
        const std::string* seq = (seqs && idx < seqs->size()) ? &(*seqs)[idx] : nullptr;
        const int RL = (int)r.size();
        // R2 with mate already encoded in this block → mate quality available
        bool has_mate = is_pe && (idx & 1) && idx > 0 && seen.count(idx - 1) &&
                        (idx - 1) < rq.size();
        const std::vector<uint8_t>* mq = has_mate ? &rq[idx - 1] : nullptr;
        CtxState s;
        for (int j = 0; j < RL; ++j) {
            int q = r[j];
            if (q > qmax) q = qmax;
            int is_dev = (has_dev && j < (int)dev_sets[idx].size() && dev_sets[idx][j]) ? 1 : 0;
            int qmb = 4; // no-mate / R1 / SE sentinel
            if (mq && j < (int)mq->size()) {
                int qm = (*mq)[j];
                qmb = qm < 10 ? 0 : qm < 20 ? 1 : qm < 30 ? 2 : 3;
            }
            uint64_t pk, ck; keys(s, j, L, is_dev, qmax, seq, pk, ck, cfg.use_seq, qmb);
            Freq& FP = lo.get(pk);
            Freq& FC = tbl.get(ck);
            uint32_t tot = build_cf(FC, FP, qmax, cf, cfg.qcm_k);
            uint32_t cum = 0; for (int i = 0; i < q; ++i) cum += cf[i];
            enc.encode(cum, cf[q], tot);
            FC.update(q, cfg.freq_inc); FP.update(q, cfg.freq_inc);
            s.advance(q);
        }
        seen.insert(idx);
    }
    enc.flush();
    return std::move(enc.out);
}

static bool decode_block(
    const uint8_t* p, const uint8_t* e,
    const std::vector<uint32_t>& order, size_t begin, size_t end,
    const std::vector<std::vector<bool>>& dev_sets, int L,
    std::vector<std::vector<uint8_t>>& rq_out,
    const std::vector<std::string>* seqs, const std::vector<int>* rlens, int qmax, QcmCfg cfg,
    bool is_pe) {
    FreqMap tbl(15), lo(12);
    RangeDec dec(p, e);
    uint32_t cf[QA];
    // Tracks which indices in this block have been fully decoded (for mate-context lookup).
    std::vector<bool> decoded(rq_out.size(), false);
    for (size_t oi = begin; oi < end; ++oi) {
        uint32_t idx = order[oi];
        if (idx >= rq_out.size()) return false;
        auto& r = rq_out[idx];
        const int RL = (rlens && idx < rlens->size()) ? (*rlens)[idx] : L;
        r.assign((size_t)RL, 0);
        const bool has_dev = idx < dev_sets.size();
        const std::string* seq = (seqs && idx < seqs->size()) ? &(*seqs)[idx] : nullptr;
        // R2 with R1 already decoded in this block → mate quality available
        bool has_mate = is_pe && (idx & 1) && idx > 0 && idx - 1 < decoded.size() &&
                        decoded[idx - 1] && !rq_out[idx - 1].empty();
        const std::vector<uint8_t>* mq = has_mate ? &rq_out[idx - 1] : nullptr;
        CtxState s;
        for (int j = 0; j < RL; ++j) {
            int is_dev = (has_dev && j < (int)dev_sets[idx].size() && dev_sets[idx][j]) ? 1 : 0;
            int qmb = 4;
            if (mq && j < (int)mq->size()) {
                int qm = (*mq)[j];
                qmb = qm < 10 ? 0 : qm < 20 ? 1 : qm < 30 ? 2 : 3;
            }
            uint64_t pk, ck; keys(s, j, L, is_dev, qmax, seq, pk, ck, cfg.use_seq, qmb);
            Freq& FP = lo.get(pk);
            Freq& FC = tbl.get(ck);
            uint32_t tot = build_cf(FC, FP, qmax, cf, cfg.qcm_k);
            uint32_t target = dec.getfreq(tot);
            uint32_t cum = 0; int q = 0;
            while (q <= qmax && cum + cf[q] <= target) { cum += cf[q]; ++q; }
            if (q > qmax) return false;
            dec.update(cum, cf[q]);
            FC.update(q, cfg.freq_inc); FP.update(q, cfg.freq_inc);
            r[j] = (uint8_t)q;
            s.advance(q);
        }
        decoded[idx] = true;
    }
    return true;
}

// Choose block count from read count and hardware threads. More reads → more
// blocks (better parallel speedup); tiny inputs stay single-block (no overhead,
// no ratio loss). Each block resets its model, so N blocks cost N cold-starts —
// with ≥40k reads/block that is a negligible fraction of the stream.
static int choose_nblocks(size_t n) {
    if (getenv("ARCS_QUAL_NOPAR")) return 1;
    if (const char* e = getenv("ARCS_QUAL_BLOCKS")) { int v = atoi(e); if (v>=1 && v<=64) return v; }
    // Scale with available cores (competitors run multithreaded by default). Each
    // block re-warms its model, so we keep ≥18k reads/block to bound the ratio cost
    // (~0.3%/block); the speed win (parallel encode + parallel decode) dominates.
    // ARCS_QUAL_BLOCKS=1 (or ARCS_QUAL_NOPAR) restores the absolute-best-ratio path.
    unsigned hc = std::thread::hardware_concurrency(); if (!hc) hc = 4;
    size_t by_reads = n / 18000; if (by_reads < 1) by_reads = 1;
    int nb = (int)std::min<size_t>(hc, by_reads);
    if (nb < 1) nb = 1; if (nb > 64) nb = 64;
    return nb;
}
static inline void put_u32(std::vector<uint8_t>& o, uint32_t v) {
    o.push_back((uint8_t)(v & 0xFF));      o.push_back((uint8_t)((v >> 8) & 0xFF));
    o.push_back((uint8_t)((v >> 16) & 0xFF)); o.push_back((uint8_t)((v >> 24) & 0xFF));
}
static inline uint32_t get_u32(const uint8_t* p) {
    return (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24);
}

} // namespace

std::vector<uint8_t> qual_cm_encode(
    const std::vector<std::vector<uint8_t>>& rq,
    const std::vector<uint32_t>&             order,
    const std::vector<std::vector<bool>>&    dev_sets,
    int                                      L,
    const std::vector<std::string>*          seqs,
    bool                                     is_pe) {

    int qmax = 0;
    for (uint32_t idx : order) {
        if (idx >= rq.size()) continue;
        const auto& r = rq[idx];
        // Per-read length, NOT global L: variable-length data (esp. a short/empty
        // reads[0] that makes L tiny) must still scan every read's full quality, or
        // qmax is undercounted and high qualities get clamped away (data loss).
        for (int j = 0; j < (int)r.size(); ++j) if (r[j] > qmax) qmax = r[j];
    }
    if (qmax > QA - 1) qmax = QA - 1;
    QcmCfg cfg = qcm_cfg_from_env();
    cfg.use_seq = (seqs != nullptr) ? 1 : 0;

    const size_t n = order.size();
    int nb = choose_nblocks(n);
    std::vector<size_t> bnd(nb + 1);
    for (int b = 0; b <= nb; ++b) bnd[b] = n * (size_t)b / (size_t)nb;
    // For PE mate-context, block boundaries must fall on even indices so R1 and R2
    // from each pair stay in the same block (mate quality is local to the block).
    if (is_pe && nb > 1) {
        for (int b = 1; b < nb; ++b)
            if (bnd[b] & 1) ++bnd[b]; // round up to even
    }

    std::vector<std::vector<uint8_t>> parts((size_t)nb);
    if (nb == 1) {
        parts[0] = encode_block(rq, order, bnd[0], bnd[1], dev_sets, L, seqs, qmax, cfg, is_pe);
    } else {
        std::vector<std::thread> th;
        th.reserve((size_t)nb);
        for (int b = 0; b < nb; ++b)
            th.emplace_back([&, b] {
                parts[(size_t)b] = encode_block(rq, order, bnd[b], bnd[b+1], dev_sets, L, seqs, qmax, cfg, is_pe);
            });
        for (auto& t : th) t.join();
    }

    // Header: qmax, freq_inc, use_seq_flags, qcm_k, nblocks, [nblocks × uint32 length], payloads…
    // use_seq_flags: bit0 = use_seq, bit1 = is_pe (inter-mate quality context)
    std::vector<uint8_t> out;
    size_t tot = 5 + (size_t)nb * 4;
    for (auto& p : parts) tot += p.size();
    out.reserve(tot);
    out.push_back((uint8_t)qmax);
    out.push_back((uint8_t)cfg.freq_inc);
    out.push_back((uint8_t)((cfg.use_seq ? 0x01 : 0) | (is_pe ? 0x02 : 0)));
    out.push_back((uint8_t)cfg.qcm_k);
    out.push_back((uint8_t)nb);
    for (auto& p : parts) put_u32(out, (uint32_t)p.size());
    for (auto& p : parts) out.insert(out.end(), p.begin(), p.end());
    return out;
}

bool qual_cm_decode(
    const std::vector<uint8_t>&              data,
    const std::vector<uint32_t>&             order,
    const std::vector<std::vector<bool>>&    dev_sets,
    int                                      L,
    std::vector<std::vector<uint8_t>>&       rq_out,
    const std::vector<std::string>*          seqs,
    const std::vector<int>*                  rlens) {

    if (data.size() < 5) return false;
    const int qmax = data[0];
    QcmCfg cfg;
    cfg.freq_inc = data[1];
    const uint8_t use_seq_flags = data[2];
    cfg.use_seq  = (use_seq_flags & 0x01) ? 1 : 0;
    const bool is_pe = (use_seq_flags & 0x02) != 0;
    cfg.qcm_k    = data[3];
    const int nb = data[4];
    if (qmax < 0 || qmax > QA - 1 || cfg.freq_inc < 1 || cfg.freq_inc > 200 ||
        use_seq_flags > 3 || cfg.qcm_k < 1 || cfg.qcm_k > 255 ||
        nb < 1 || nb > 64) return false;
    if (cfg.use_seq && !seqs) return false;

    const size_t hdr = 5 + (size_t)nb * 4;
    if (data.size() < hdr) return false;
    std::vector<uint32_t> len((size_t)nb);
    size_t payload = hdr;
    for (int b = 0; b < nb; ++b) { len[(size_t)b] = get_u32(data.data() + 5 + (size_t)b*4); payload += len[(size_t)b]; }
    if (payload != data.size()) return false;

    const size_t n = order.size();
    std::vector<size_t> bnd(nb + 1);
    for (int b = 0; b <= nb; ++b) bnd[b] = n * (size_t)b / (size_t)nb;
    // Mirror the encoder's even-boundary rounding for PE blocks
    if (is_pe && nb > 1) {
        for (int b = 1; b < nb; ++b)
            if (bnd[b] & 1) ++bnd[b];
    }

    // Block byte offsets
    std::vector<size_t> off((size_t)nb + 1);
    off[0] = hdr;
    for (int b = 0; b < nb; ++b) off[(size_t)b+1] = off[(size_t)b] + len[(size_t)b];

    std::vector<uint8_t> ok((size_t)nb, 0);
    auto run = [&](int b) {
        const uint8_t* p = data.data() + off[(size_t)b];
        const uint8_t* e = data.data() + off[(size_t)b+1];
        ok[(size_t)b] = decode_block(p, e, order, bnd[b], bnd[b+1],
                                     dev_sets, L, rq_out, seqs, rlens, qmax, cfg, is_pe) ? 1 : 0;
    };
    if (nb == 1) {
        run(0);
    } else {
        std::vector<std::thread> th;
        th.reserve((size_t)nb);
        for (int b = 0; b < nb; ++b) th.emplace_back(run, b);
        for (auto& t : th) t.join();
    }
    for (int b = 0; b < nb; ++b) if (!ok[(size_t)b]) return false;
    return true;
}
