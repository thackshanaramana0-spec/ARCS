#include "caller.h"
#include "chain_encoder.h"     // CallData
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

// Faithful C++ port of the frozen reference-free bubble-caller (scripts/rf_frozen2.py):
//   HDMAX=2, MAF=0.20, DHI=2.5, KHI=1.1, MC=3, TRI=0.12, HALF=15.
// Pileup is built from ARCS's own read→consensus placements (CallData); k-mers are
// counted internally. Read-placement confidence uses MAPQ=60-6*mm (mm = mismatches
// vs contig); the mapq<20 filter of the reference pipeline == skip reads with mm>=7.

namespace {

constexpr int HDMAX = 2, MC = 3, HALF = 15;
constexpr double MAF = 0.20, DHI = 2.5, KHI = 1.1, TRI = 0.12;

// ── Indel calling parameters ──────────────────────────────────────────────────
// A het indel present in ~half the reads is invisible to the ungapped SNV pileup:
// the contig consensus carries one allele, and the reads of the OTHER haplotype
// frameshift against it and get dropped by the mm>=7 (mapq<20) filter. We recover
// exactly those reads with a banded fitting alignment, collect the indel they carry,
// and call it when enough reads agree.
constexpr int MAXINDEL     = 30;    // longest indel event we call
constexpr double INDEL_AF_HI = 0.85; // above this it is ~homozygous (contig would carry it)

inline int b2i(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; }
    return -1;
}
inline char compl_base(char c) {
    switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; }
    return 'N';
}
std::string rc_str(const std::string& s) {
    std::string r(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) r[s.size() - 1 - i] = compl_base(s[i]);
    return r;
}

// ── Packed canonical 31-mer counting (replaces DSK) ───────────────────────────
inline bool pack31(const char* s, uint64_t& out) {
    uint64_t v = 0;
    for (int i = 0; i < 31; ++i) { int b = b2i(s[i]); if (b < 0) return false; v = (v << 2) | (uint64_t)b; }
    out = v; return true;
}
inline uint64_t rc31(uint64_t v) {
    uint64_t r = 0;
    for (int i = 0; i < 31; ++i) { r = (r << 2) | (3u - (v & 3u)); v >>= 2; }
    return r;
}
inline uint64_t canon31(uint64_t v) { uint64_t r = rc31(v); return v < r ? v : r; }

inline uint64_t colkey(uint32_t cid, uint32_t pos) { return ((uint64_t)cid << 32) | pos; }

// ── STR / tandem-repeat detection ────────────────────────────────────────────
// Returns true when `seq` is composed of a short tandem unit (≤4 bp) that also
// appears in the left-flank context `flank`. Such bubbles are STR length changes:
// they require stronger anchor evidence before we emit a call.
static bool is_str_event(const std::string& seq, const std::string& flank) {
    if (seq.size() < 2) return false;
    for (int u = 1; u <= std::min((int)seq.size(), 4); ++u) {
        // check that seq is a perfect tandem of the unit
        bool tandem = true;
        for (size_t i = 0; i < seq.size(); i += (size_t)u) {
            size_t rem = std::min((size_t)u, seq.size() - i);
            if (seq.substr(i, rem) != seq.substr(0, rem)) { tandem = false; break; }
        }
        if (!tandem) continue;
        // check the same unit appears in the left flank (confirms STR context)
        std::string unit = seq.substr(0, (size_t)u);
        if (flank.find(unit) != std::string::npos) return true;
    }
    return false;
}

// ── Contig-bubble indel extraction ───────────────────────────────────────────
// Het indels segregate the two haplotypes into separate contigs (the frameshift
// breaks the k-mer chain), so they appear as a BUBBLE between two contigs that share
// flanking sequence, not as a within-read event. Given a shared anchor (A[pA]==B[qB]),
// walk both contigs from the anchor: they match, diverge by exactly one indel, then
// re-converge (a right flank of FLANK identical bases). A is the reference haplotype;
// the indel is expressed in A's coordinates.
struct Bubble { int type; uint32_t apos; int len; std::string ins; bool ok = false; };

inline bool flank_match(const std::string& A, size_t ai, const std::string& B, size_t bi, int F) {
    if (ai + (size_t)F > A.size() || bi + (size_t)F > B.size()) return false;
    for (int k = 0; k < F; ++k) if (A[ai + (size_t)k] != B[bi + (size_t)k]) return false;
    return true;
}

Bubble extract_bubble(const std::string& A, uint32_t pA, const std::string& B, uint32_t qB,
                      int MAXINDEL, int FLANK) {
    Bubble r;
    const size_t la = A.size(), lb = B.size();
    // extend the matched flank from the anchor
    uint32_t d = 0;
    while (pA + d < la && qB + d < lb && A[pA + d] == B[qB + d]) ++d;
    if (d == 0) return r;                       // no anchoring match (need a left flank)
    if (pA + d >= la || qB + d >= lb) return r;  // ran off an end: no bubble
    // insertion in B (alt haplotype gained bases): A continues after B skips g bases
    for (int g = 1; g <= MAXINDEL; ++g) {
        if (qB + d + (uint32_t)g + (uint32_t)FLANK > lb) break;
        if (flank_match(A, pA + d, B, qB + d + (uint32_t)g, FLANK)) {
            r.type = 1; r.apos = pA + d; r.len = g;
            r.ins = B.substr(qB + d, (size_t)g); r.ok = true; return r;
        }
    }
    // deletion in B (alt haplotype lost bases): B continues after A skips g bases
    for (int g = 1; g <= MAXINDEL; ++g) {
        if (pA + d + (uint32_t)g + (uint32_t)FLANK > la) break;
        if (flank_match(A, pA + d + (uint32_t)g, B, qB + d, FLANK)) {
            r.type = 0; r.apos = pA + d; r.len = g; r.ok = true; return r;
        }
    }
    return r;
}

// ── Cross-contig SNV bubble ───────────────────────────────────────────────────
// A 1-substitution bubble between two haplotype contigs: they share flanking
// sequence but differ at exactly one base. The assembler's correct haplotype
// separation makes this invisible to the pileup; it shows up here as a bubble
// where the two contigs diverge at one base then immediately re-converge.
struct SnvBubble { bool ok = false; uint32_t apos; char ref_base, alt_base; };

inline SnvBubble extract_snv_bubble(const std::string& A, uint32_t pA,
                                     const std::string& B, uint32_t qB, int FLANK) {
    SnvBubble r;
    const size_t la = A.size(), lb = B.size();
    uint32_t d = 0;
    while (pA + d < la && qB + d < lb && A[pA + d] == B[qB + d]) ++d;
    if (pA + d >= la || qB + d >= lb) return r;          // ran off an end
    char ra = A[pA + d], rb = B[qB + d];
    if (b2i(ra) < 0 || b2i(rb) < 0) return r;            // non-ACGT at bubble
    // right flank: next FLANK bases must be identical (ensures a clean single-SNV bubble)
    if (pA + d + 1u + (uint32_t)FLANK > (uint32_t)la) return r;
    if (qB + d + 1u + (uint32_t)FLANK > (uint32_t)lb) return r;
    if (!flank_match(A, pA + d + 1u, B, qB + d + 1u, FLANK)) return r;
    r.ok = true; r.apos = pA + d; r.ref_base = ra; r.alt_base = rb;
    return r;
}

// packed forward 25-mer (50 bits); returns false on non-ACGT
inline bool pack25(const char* s, uint64_t& out) {
    uint64_t v = 0;
    for (int i = 0; i < 25; ++i) { int b = b2i(s[i]); if (b < 0) return false; v = (v << 2) | (uint64_t)b; }
    out = v; return true;
}
inline uint64_t rc25(uint64_t v) {
    uint64_t r = 0;
    for (int i = 0; i < 25; ++i) { r = (r << 2) | (3u - (v & 3u)); v >>= 2; }
    return r;
}

double loglik(const std::vector<std::pair<int,int>>& rl, const int* al, int m) {
    double ll = 0.0;
    for (auto& bq : rl) {
        double e = std::pow(10.0, -bq.second / 10.0);
        double p = 0.0;
        for (int a = 0; a < m; ++a) p += (bq.first == al[a]) ? (1.0 - e) : (e / 3.0);
        ll += std::log(std::max(p / m, 1e-300));
    }
    return ll;
}

} // namespace

int run_variant_call(const std::vector<Read>& reads, const CallData& cd,
                     const std::string& out_vcf) {
    if (!cd.valid) { fprintf(stderr, "caller: no placement data\n"); return -1; }
    const size_t n = reads.size();

    // Ploidy: default 2 = diploid (the frozen het-SNV path, byte-identical). ARCS_PLOIDY=k
    // (3=triploid, 4=tetraploid) allows up to k co-occurring alleles at a bubble column and
    // emits multi-allelic SNV records. Only k>2 changes any behaviour.
    int PLOIDY = 2;
    if (const char* pe = std::getenv("ARCS_PLOIDY")) { int v = atoi(pe); if (v >= 2 && v <= 4) PLOIDY = v; }

    // Optional: dump contig sequences (for the self-contained lift test placer).
    if (const char* dp = std::getenv("ARCS_DUMP_CONTIGS")) {
        FILE* df = fopen(dp, "w");
        if (df) {
            for (size_t ci = 0; ci < cd.contigs.size(); ++ci)
                fprintf(df, ">contig_%zu\n%s\n", ci, cd.contigs[ci].c_str());
            fclose(df);
        }
    }

    // ── 1. Internal canonical-31-mer counts ──
    std::unordered_map<uint64_t, uint32_t> kc;
    kc.reserve(1u << 21);
    for (const auto& rd : reads) {
        const std::string& s = rd.seq;
        if (s.size() < 31) continue;
        for (size_t i = 0; i + 31 <= s.size(); ++i) {
            uint64_t v;
            if (pack31(s.data() + i, v)) ++kc[canon31(v)];
        }
    }
    // H = expected k-mer count at a single-copy het site (haploid depth estimate).
    // Count-frequency histogram mode: single-copy k-mers are the most numerous as
    // DISTINCT entries (~320k regardless of coverage), outnumbering repeat k-mers.
    // k²-weighting fails because repeat counts (e.g. 13000²) swamp signal counts
    // (240²).  Pileup-median fails when assembly fragments (low pileup → low H →
    // KHI threshold rejects all real k-mers at high coverage → 0 calls).
    uint32_t H = 30;
    {
        uint32_t cnt_max = 0;
        for (auto& [km, cnt] : kc) cnt_max = std::max(cnt_max, cnt);
        cnt_max = std::min(cnt_max, 5000u);
        if (cnt_max >= 4) {
            uint32_t bw = std::max(1u, cnt_max / 200u);
            uint32_t nb = cnt_max / bw + 2;
            std::vector<uint64_t> bkt(nb, 0);
            for (auto& [km, cnt] : kc)
                if (cnt >= 2 && cnt <= cnt_max) bkt[cnt / bw]++;
            // Find first local minimum in first third (valley between error and signal peaks).
            size_t valley = 2; // fallback: skip first 2 buckets (error-adjacent low counts)
            for (size_t b = 1; b + 1 < nb / 3 && b + 1 < bkt.size(); ++b) {
                if (bkt[b] < bkt[b - 1] && bkt[b] < bkt[b + 1]) { valley = b; break; }
            }
            // Highest peak after valley = single-copy k-mer count = H.
            size_t peak = (valley + 1 < bkt.size()) ? valley + 1 : valley;
            for (size_t b = valley + 1; b < bkt.size(); ++b)
                if (bkt[b] > bkt[peak]) peak = b;
            uint32_t H_hist = (uint32_t)((peak + 0.5) * bw);
            if (H_hist >= 5) H = H_hist;
        }
    }
    auto kcount = [&](const std::string& km) -> uint32_t {
        if (km.size() != 31) return 0;
        uint64_t v; if (!pack31(km.data(), v)) return 0;
        auto it = kc.find(canon31(v)); return it == kc.end() ? 0u : it->second;
    };

    // ── 2. Pileup from placements (contig frame). Skip reads with mm>=7 (mapq<20). ──
    // Store per read its contig-frame sequence + start, for the flank pass.
    struct Rec { uint32_t cid, pos; std::string seq, qual; };
    std::vector<Rec> recs; recs.reserve(n);
    std::unordered_map<uint64_t, std::vector<std::pair<int,int>>> col;  // key → [(base,qual)]
    col.reserve(1u << 20);

    for (size_t oi = 0; oi < n; ++oi) {
        uint32_t cid = cd.read_cid[oi], pos = cd.read_pos[oi];
        if (cid >= cd.contigs.size()) continue;
        const std::string& cc = cd.contigs[cid];
        bool rc = cd.read_rc[oi] != 0;
        std::string seq = rc ? rc_str(reads[oi].seq) : reads[oi].seq;     // contig frame
        std::string qual = reads[oi].qual;
        if (rc) std::reverse(qual.begin(), qual.end());
        const int rl = (int)seq.size();
        // placement confidence (mm vs contig) → mapq filter
        int mm = 0;
        for (int j = 0; j < rl; ++j) {
            uint32_t p = pos + (uint32_t)j;
            if (p >= cc.size()) { mm = rl; break; }
            char a = seq[(size_t)j];
            if (b2i(a) >= 0 && cc[p] != a) ++mm;
        }
        if (mm >= 7) continue;                                            // mapq<20
        recs.push_back({cid, pos, seq, qual});
        for (int j = 0; j < rl; ++j) {
            uint32_t p = pos + (uint32_t)j;
            if (p >= cc.size()) break;
            int b = b2i(seq[(size_t)j]);
            if (b < 0) continue;
            int q = (j < (int)qual.size() && qual[(size_t)j] >= 33) ? (qual[(size_t)j] - 33) : 0;
            col[colkey(cid, p)].push_back({b, q});
        }
    }

    // ── 3. Candidate columns: major/minor allele, depth, quality, likelihood ratio ──
    struct Cand { int M, mn; int cnt[4]; int d; };
    std::unordered_map<uint64_t, Cand> C;
    for (auto& kv : col) {
        auto& rl = kv.second;
        int d = (int)rl.size();
        if (d < 6) continue;
        int cnt[4] = {0,0,0,0};
        for (auto& bq : rl) cnt[bq.first]++;
        int o[4] = {0,1,2,3};
        std::sort(o, o + 4, [&](int a, int b){ return cnt[a] > cnt[b]; });
        int M = o[0], mn = o[1];
        if (cnt[mn] < 2) continue;
        std::vector<int> mq;
        for (auto& bq : rl) if (bq.first == mn) mq.push_back(bq.second);
        std::sort(mq.begin(), mq.end());
        double medq = mq.empty() ? 0 : (mq.size() % 2 ? (double)mq[mq.size()/2]
                                        : (mq[mq.size()/2 - 1] + mq[mq.size()/2]) / 2.0);
        if (medq < 20) continue;
        int alMn[2] = {M, mn}, alM[1] = {M};
        if (10.0 * (loglik(rl, alMn, 2) - loglik(rl, alM, 1)) / std::log(10.0) < 10) continue;
        Cand c; c.M = M; c.mn = mn; c.cnt[0]=cnt[0]; c.cnt[1]=cnt[1]; c.cnt[2]=cnt[2]; c.cnt[3]=cnt[3]; c.d = d;
        C[kv.first] = c;
    }
    if (C.empty()) { fprintf(stderr, "caller: 0 candidates\n"); }

    if (const char* he = std::getenv("ARCS_HAPLOID_COV")) { int v = atoi(he); if (v > 0) H = (uint32_t)v; }

    // ── 4. Flank pass: major/minor 31-mer flanks + per-offset base votes ──
    std::unordered_map<uint64_t, std::unordered_map<std::string,int>> FLmaj, FLmin;
    std::unordered_map<uint64_t, std::array<std::array<int,4>,31>> majoff, minoff;  // off idx = off+15
    for (auto& r : recs) {
        const int rl = (int)r.seq.size();
        for (int j = 0; j < rl; ++j) {
            uint64_t key = colkey(r.cid, r.pos + (uint32_t)j);
            auto ci = C.find(key);
            if (ci == C.end()) continue;
            int b = b2i(r.seq[(size_t)j]);
            if (b < 0) continue;
            bool isM = (b == ci->second.M), ismn = (b == ci->second.mn);
            if (!isM && !ismn) continue;
            if (j - HALF >= 0 && j + HALF + 1 <= rl)
                (isM ? FLmaj : FLmin)[key][r.seq.substr((size_t)(j - HALF), 2 * HALF + 1)]++;
            auto& grp = isM ? majoff[key] : minoff[key];
            for (int off = -HALF; off <= HALF; ++off) {
                if (off == 0) continue;
                int p = j + off;
                if (p >= 0 && p < rl) { int bb = b2i(r.seq[(size_t)p]); if (bb >= 0) grp[(size_t)(off + HALF)][(size_t)bb]++; }
            }
        }
    }

    auto most_common = [](const std::unordered_map<std::string,int>& m, int& cntout) -> std::string {
        std::string best; int bc = 0;
        for (auto& kv : m) if (kv.second > bc) { bc = kv.second; best = kv.first; }
        cntout = bc; return best;
    };

    // ── 5. Frozen filters → kept calls ──
    std::vector<std::pair<uint32_t,uint32_t>> kept;   // (cid, pos)
    for (auto& kv : C) {
        uint64_t key = kv.first; const Cand& c = kv.second;
        // bubble cleanliness: haplotype-flank Hamming
        int diffs = 0;
        auto mao = majoff.find(key); auto mio = minoff.find(key);
        if (mao != majoff.end() && mio != minoff.end()) {
            for (int oi2 = 0; oi2 < 31; ++oi2) {
                if (oi2 == HALF) continue;
                const auto& ma = mao->second[(size_t)oi2]; const auto& mi = mio->second[(size_t)oi2];
                int sma = ma[0]+ma[1]+ma[2]+ma[3], smi = mi[0]+mi[1]+mi[2]+mi[3];
                if (sma < 3 || smi < 3) continue;
                int am = 0, im = 0;
                for (int b = 1; b < 4; ++b) { if (ma[b] > ma[am]) am = b; if (mi[b] > mi[im]) im = b; }
                if (am != im) ++diffs;
            }
        }
        if (diffs > HDMAX) continue;
        if (c.d > (int)(DHI * H * 1.25 + 0.5)) continue;
        int cmajN = 0, cminN = 0;
        std::string fmaj, fmin;
        auto mj = FLmaj.find(key); auto mn2 = FLmin.find(key);
        if (mj != FLmaj.end()) { int cc; fmaj = most_common(mj->second, cc); }
        if (mn2 != FLmin.end()) { int cc; fmin = most_common(mn2->second, cc); }
        uint32_t cmaj = fmaj.empty() ? 0 : kcount(fmaj);
        uint32_t cmin = fmin.empty() ? 0 : kcount(fmin);
        (void)cmajN; (void)cminN;
        if ((double)std::max(cmaj, cmin) > KHI * H) continue;
        if ((double)c.cnt[c.mn] / c.d < MAF) continue;
        int o[4] = {0,1,2,3};
        std::sort(o, o + 4, [&](int a, int b){ return c.cnt[a] > c.cnt[b]; });
        // Allow up to PLOIDY co-occurring alleles: reject a column only if the
        // (PLOIDY+1)-th allele exceeds the noise fraction. PLOIDY=2 → o[2] (identical to
        // the frozen diploid triallelic-rejection); PLOIDY=3 → o[3]; PLOIDY=4 → no 5th base.
        if (PLOIDY < 4 && (double)c.cnt[o[PLOIDY]] > TRI * c.d) continue;
        int mcnt = 0;
        if (mn2 != FLmin.end()) { std::string s = most_common(mn2->second, mcnt); (void)s; }
        if (mcnt < MC) continue;
        kept.push_back({(uint32_t)(key >> 32), (uint32_t)(key & 0xFFFFFFFFu)});
    }

    // ── 6a. SNV records (contig coordinates; reference-free) ──
    struct OutRec { uint32_t cid, pos; std::string ref, alt, info; };
    std::vector<OutRec> orecs;
    std::sort(kept.begin(), kept.end());
    for (auto& kp : kept) {
        uint32_t cid = kp.first, pos = kp.second;
        const std::string& cc = cd.contigs[cid];
        char ref = pos < cc.size() ? cc[pos] : 'N';
        Cand& c = C[colkey(cid, pos)];
        int refb = b2i(ref);
        if (PLOIDY == 2) {
            // ── frozen diploid path (byte-identical): ALT = the M/mn allele ≠ REF ──
            int altb = (c.M != refb) ? c.M : c.mn;
            char alt = "ACGT"[altb & 3];
            char info[64];
            snprintf(info, sizeof info, "AF=%.3f;DP=%d", (double)c.cnt[c.mn] / c.d, c.d);
            orecs.push_back({cid, pos + 1u, std::string(1, ref), std::string(1, alt), info});
        } else {
            // ── polyploid: emit every top-PLOIDY allele (freq≥MAF) that differs from REF ──
            int o[4] = {0,1,2,3};
            std::sort(o, o + 4, [&](int a, int b){ return c.cnt[a] > c.cnt[b]; });
            std::string alts; int nalt = 0;
            for (int t = 0; t < PLOIDY; ++t) {
                int a = o[t];
                if (a == refb) continue;
                if ((double)c.cnt[a] / c.d < MAF) continue;
                if (!alts.empty()) alts += ",";
                alts += "ACGT"[a & 3];
                ++nalt;
            }
            if (nalt == 0) continue;             // all supported alleles == REF (not a variant)
            char info[80];
            snprintf(info, sizeof info, "AF=%.3f;DP=%d;PLOIDY=%d", (double)c.cnt[c.mn] / c.d, c.d, PLOIDY);
            orecs.push_back({cid, pos + 1u, std::string(1, ref), alts, info});
        }
    }
    size_t n_snv = orecs.size();

    // ── 6b. Indel pass: het indels are BUBBLES between haplotype contigs ──
    // A het indel frameshifts one haplotype's reads, so the assembler splits them into
    // a separate contig. The indel therefore shows up as two contigs that share flanking
    // sequence but differ by an inserted/deleted stretch. We find shared 25-mer anchors
    // between contig pairs, extract the clean single-indel bubble, and express it on the
    // longer (reference) haplotype. Support/allele-balance come from each contig's read
    // coverage. (default on; ARCS_NO_INDELS disables to reproduce the frozen SNV caller.)
    size_t n_indel = 0, n_xsnv = 0;
    if (!std::getenv("ARCS_NO_INDELS") && cd.contigs.size() >= 2) {
        constexpr int BK = 25, FLANK = 15;
        // per-contig read coverage (placed reads, ungapped span)
        std::vector<std::vector<uint16_t>> cov(cd.contigs.size());
        for (size_t ci = 0; ci < cd.contigs.size(); ++ci)
            cov[ci].assign(cd.contigs[ci].size(), 0);
        for (size_t oi = 0; oi < n; ++oi) {
            uint32_t cid = cd.read_cid[oi], pos = cd.read_pos[oi];
            if (cid >= cd.contigs.size()) continue;
            int rl = (int)reads[oi].seq.size();
            for (int j = 0; j < rl; ++j) {
                uint32_t p = pos + (uint32_t)j;
                if (p < cov[cid].size() && cov[cid][p] < 60000) ++cov[cid][p];
            }
        }
        // canonical 25-mer → occurrences (cid, pos, orient); orient=0 fwd==canonical else 1
        std::unordered_map<uint64_t, std::vector<std::tuple<uint32_t,uint32_t,uint8_t>>> kidx;
        kidx.reserve(1u << 20);
        for (size_t ci = 0; ci < cd.contigs.size(); ++ci) {
            const std::string& c = cd.contigs[ci];
            for (size_t i = 0; i + BK <= c.size(); ++i) {
                uint64_t v; if (!pack25(c.data() + i, v)) continue;
                uint64_t rcv = rc25(v), canon = v < rcv ? v : rcv;
                kidx[canon].push_back(std::make_tuple((uint32_t)ci, (uint32_t)i, (uint8_t)(v <= rcv ? 0 : 1)));
            }
        }
        // aggregate bubbles: key (refcid, apos, type, len, ins) → (anchors, altcid)
        struct Agg { int anchors = 0; uint32_t altcid = 0, altpos = 0; };
        std::map<std::tuple<uint32_t,uint32_t,int,int,std::string>, Agg> im;
        int dbg_pairs = 0;
        for (auto& kv : kidx) {
            auto& occ = kv.second;
            if (occ.size() != 2) continue;                       // unique cross-contig anchor
            uint32_t ca, pa, cb, pb; uint8_t oa, ob;
            std::tie(ca, pa, oa) = occ[0]; std::tie(cb, pb, ob) = occ[1];
            if (ca == cb) continue;                              // same contig (tandem) — skip
            ++dbg_pairs;
            // reference = longer contig (better lift), other = alt, oriented to match ref
            uint32_t rc_, rp, ac, ap; uint8_t ro, ao;
            if (cd.contigs[ca].size() >= cd.contigs[cb].size()) { rc_=ca; rp=pa; ro=oa; ac=cb; ap=pb; ao=ob; }
            else { rc_=cb; rp=pb; ro=ob; ac=ca; ap=pa; ao=oa; }
            const std::string& R = cd.contigs[rc_];
            bool opp = (ro != ao);
            std::string Aalt = opp ? rc_str(cd.contigs[ac]) : cd.contigs[ac];
            uint32_t qB = opp ? (uint32_t)(cd.contigs[ac].size() - ap - BK) : ap;
            Bubble bub = extract_bubble(R, rp, Aalt, qB, MAXINDEL, FLANK);
            if (!bub.ok || bub.apos == 0) continue;
            auto& a = im[std::make_tuple(rc_, bub.apos, bub.type, bub.len, bub.ins)];
            a.anchors++; a.altcid = ac; a.altpos = ap;
        }
        if (std::getenv("ARCS_INDEL_DEBUG"))
            fprintf(stderr, "[IDBG] contigs=%zu cross-anchors=%d bubbles=%zu\n",
                    cd.contigs.size(), dbg_pairs, im.size());
        // windowed mean coverage — robust to single-base dropouts near contig ends
        auto covwin = [&](uint32_t ci, uint32_t p) -> int {
            const auto& cv = cov[ci]; if (cv.empty()) return 0;
            int lo = (int)p - 15, hi = (int)p + 15, s = 0, cnt = 0;
            for (int x = std::max(0, lo); x <= std::min((int)cv.size() - 1, hi); ++x) { s += cv[x]; ++cnt; }
            return cnt ? s / cnt : 0;
        };
        // per-contig median coverage — the alt haplotype is a short fragment whose
        // per-position support is what a single anchor kmer can't localize; its overall
        // depth is the honest measure of how strongly reads support that allele.
        std::vector<int> medcov(cd.contigs.size(), 0);
        for (size_t ci = 0; ci < cd.contigs.size(); ++ci) {
            if (cov[ci].empty()) continue;
            std::vector<uint16_t> v = cov[ci];
            std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
            medcov[ci] = v[v.size() / 2];
        }
        for (auto& kv : im) {
            uint32_t cid, apos; int type, len; std::string ins;
            std::tie(cid, apos, type, len, ins) = kv.first;
            Agg& a = kv.second;
            const std::string& cc = cd.contigs[cid];
            if (apos == 0 || apos > cc.size()) continue;
            // allele balance: ref-contig coverage at the locus vs alt-contig overall depth
            int refd = covwin(cid, apos - 1);
            int altd = std::max(medcov[a.altcid], (int)covwin(a.altcid, a.altpos));
            double af = (double)altd / std::max(1, refd + altd);
            if (std::getenv("ARCS_INDEL_DEBUG"))
                fprintf(stderr, "[IDBG] bubble cid=%u apos=%u type=%d len=%d ins=%s refd=%d altd=%d af=%.3f anchors=%d altlen=%zu\n",
                        cid, apos, type, len, ins.c_str(), refd, altd, af, a.anchors, cd.contigs[a.altcid].size());
            // A bubble is heterozygous by construction (a homozygous indel yields one
            // contig, no bubble). The real requirements: both haplotypes carry read
            // support, and the shared flank is long enough to be a true bubble (not a
            // repeat coincidence). Alt-contig coverage under-counts the alt allele
            // (alt reads spill onto the shared ref contig), so no SNV-style AF band.
            if (refd < MC || altd < MC) continue;
            if (a.anchors < 3) continue;                          // solid flank, not a fluke
            // STR filter: tandem-repeat bubbles (homopolymer/STR length changes) are the
            // dominant false-positive class on real GIAB data. Require 5 anchors instead
            // of 3 when the indel sequence is a tandem unit that also appears in the left
            // flank — this eliminates collapsed-STR FPs without touching clean indels.
            {
                std::string left_flank = (apos >= 12) ? cc.substr(apos - 12, 12)
                                                       : cc.substr(0, apos);
                std::string indel_seq  = (type == 0)
                    ? ((apos + (uint32_t)len <= cc.size()) ? cc.substr(apos, (size_t)len) : "")
                    : ins;
                if (!indel_seq.empty() && is_str_event(indel_seq, left_flank) && a.anchors < 5)
                    continue;
            }
            char anchor = cc[apos - 1];
            std::string ref, alt;
            if (type == 0) {                                     // deletion in alt haplotype
                if (apos + (uint32_t)len > cc.size()) continue;
                ref = std::string(1, anchor) + cc.substr(apos, (size_t)len);
                alt = std::string(1, anchor);
            } else {                                             // insertion in alt haplotype
                ref = std::string(1, anchor);
                alt = std::string(1, anchor) + ins;
            }
            char info[96];
            snprintf(info, sizeof info, "SVTYPE=INDEL;AF=%.3f;DP=%d;ANCHORS=%d",
                     af, refd + altd, a.anchors);
            orecs.push_back({cid, apos, ref, alt, info});        // apos already 1-based anchor
            ++n_indel;
        }

        // ── 6b2. Cross-contig SNV pass (opt-in: ARCS_XSNV=1) ─────────────────
        // When the assembler correctly haplotype-separates reads, the pileup on
        // each contig shows only one allele → FN. We scan the same anchor pairs
        // used for indels, looking for a 1-substitution bubble.
        // Disabled by default: the 1-SNP bubble pattern is less specific than the
        // indel re-convergence pattern, causing net F1 regression on some windows.
        // Enable with ARCS_XSNV=1 for experimental use; tune parameters per dataset.
        if (!std::getenv("ARCS_XSNV")) goto skip_xsnv;
        //
        // Three quality gates enforce high precision:
        //  (1) ANCHORS ≥ 15: many independent 25-mer anchors in the left flank
        //      confirm the bubble is not a chance match.
        //  (2) Position-specific allele balance (AF ∈ [0.25,0.75]): both contigs
        //      must carry similar read depth at the bubble site. A collapsed contig
        //      (both haplotypes mis-assembled into one) appears as high ref-contig
        //      depth, while its spurious alt-pair also has high depth → AF ≈ 0.5
        //      but combined DP >> expected, so we also apply:
        //  (3) DP ≤ 1.8×H: for a genuine het pair each contig carries ~half the
        //      reads at that locus, so total ≈ H. Inflated DP reveals a collapsed
        //      contig falsely paired with an unrelated alt contig.
        //
        // Position-specific alt coverage (covwin at the bubble site on the alt
        // contig, not medcov across the full contig) gives a reliable AF estimate.
        {
            constexpr int XSNV_FLANK = 15;
            constexpr int XSNV_MIN_ANCHORS = 15;
            // collect pileup SNV keys (1-based VCF pos) so we don't double-emit
            std::unordered_set<uint64_t> pileup_keys;
            pileup_keys.reserve(n_snv);
            for (size_t i = 0; i < n_snv; ++i)
                pileup_keys.insert(colkey(orecs[i].cid, orecs[i].pos));
            struct XsnvAgg { int anchors = 0; uint32_t altcid = 0, altpos = 0; };
            std::map<std::tuple<uint32_t,uint32_t,char,char>, XsnvAgg> xim;
            for (auto& kv : kidx) {
                const auto& occ = kv.second;
                if (occ.size() != 2) continue;
                uint32_t cax, pax, cbx, pbx; uint8_t oax, obx;
                std::tie(cax, pax, oax) = occ[0]; std::tie(cbx, pbx, obx) = occ[1];
                if (cax == cbx) continue;
                uint32_t rcx, rpx, acx, apx; uint8_t rox, aox;
                if (cd.contigs[cax].size() >= cd.contigs[cbx].size()) {
                    rcx=cax; rpx=pax; rox=oax; acx=cbx; apx=pbx; aox=obx;
                } else {
                    rcx=cbx; rpx=pbx; rox=obx; acx=cax; apx=pax; aox=oax;
                }
                bool oppx = (rox != aox);
                std::string Bx = oppx ? rc_str(cd.contigs[acx]) : cd.contigs[acx];
                uint32_t qBx = oppx ? (uint32_t)(cd.contigs[acx].size() - apx - BK) : apx;
                SnvBubble sb = extract_snv_bubble(cd.contigs[rcx], rpx, Bx, qBx, XSNV_FLANK);
                if (!sb.ok) continue;
                // position of the bubble in alt contig's FORWARD strand
                uint32_t d = sb.apos - rpx;
                uint32_t alt_pos_fwd = oppx
                    ? (uint32_t)(cd.contigs[acx].size() - 1u) - (qBx + d)
                    : qBx + d;
                auto& xa = xim[std::make_tuple(rcx, sb.apos, sb.ref_base, sb.alt_base)];
                xa.anchors++;
                xa.altcid = acx;
                xa.altpos = alt_pos_fwd;   // last writer wins; all anchors for same bubble agree
            }
            for (auto& kv : xim) {
                uint32_t cid2, apos2; char rb2, ab2;
                std::tie(cid2, apos2, rb2, ab2) = kv.first;
                XsnvAgg& xa = kv.second;
                if (xa.anchors < XSNV_MIN_ANCHORS) continue;
                uint32_t vcf_pos2 = apos2 + 1u;              // 0-based apos → 1-based VCF POS
                if (pileup_keys.count(colkey(cid2, vcf_pos2))) continue;
                // position-specific coverage on both contigs (much more accurate than medcov)
                int refd2 = covwin(cid2, apos2);
                int altd2 = covwin(xa.altcid, xa.altpos);
                if (refd2 < MC || altd2 < MC) continue;
                int dp2 = refd2 + altd2;
                double af2 = (double)altd2 / std::max(1, dp2);
                // gate 2: allele balance — genuine het pair has ~equal depth on each contig
                if (af2 < 0.25 || af2 > 0.75) continue;
                // gate 3: total DP close to one haplotype coverage (H); inflated DP = collapsed contig
                if (dp2 > (int)(H * 1.8 + 0.5)) continue;
                char inf2[96];
                snprintf(inf2, sizeof inf2, "AF=%.3f;DP=%d;ANCHORS=%d;SOURCE=XCONTIG",
                         af2, dp2, xa.anchors);
                orecs.push_back({cid2, vcf_pos2, std::string(1, rb2), std::string(1, ab2), inf2});
                ++n_xsnv;
            }
        }
        skip_xsnv:;
    }

    // ── 6c. Write all records, sorted by (contig, pos) so the VCF is valid ──
    std::sort(orecs.begin(), orecs.end(), [](const OutRec& x, const OutRec& y){
        return x.cid != y.cid ? x.cid < y.cid : x.pos < y.pos;
    });
    FILE* f = fopen(out_vcf.c_str(), "w");
    if (!f) { fprintf(stderr, "caller: cannot open %s\n", out_vcf.c_str()); return -1; }
    fprintf(f, "##fileformat=VCFv4.2\n##source=ARCS-reffree-caller\n");
    // Declare each assembled contig so the VCF is self-contained (POS is 1-based per VCF v4.2).
    for (size_t ci = 0; ci < cd.contigs.size(); ++ci)
        fprintf(f, "##contig=<ID=contig_%zu,length=%zu>\n", ci, cd.contigs[ci].size());
    fprintf(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (auto& r : orecs)
        fprintf(f, "contig_%u\t%u\t.\t%s\t%s\t.\tPASS\t%s\n",
                r.cid, r.pos, r.ref.c_str(), r.alt.c_str(), r.info.c_str());
    fclose(f);
    fprintf(stderr, "[CALL] H=%u candidates=%zu SNVs=%zu+%zu indels=%zu -> %s\n",
            H, C.size(), n_snv, n_xsnv, n_indel, out_vcf.c_str());
    return (int)orecs.size();
}
