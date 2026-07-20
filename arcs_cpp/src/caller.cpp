#include "caller.h"
#include "chain_encoder.h"     // CallData
#include <unordered_map>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>

// Faithful C++ port of the frozen reference-free bubble-caller (scripts/rf_frozen2.py):
//   HDMAX=2, MAF=0.20, DHI=2.5, KHI=1.1, MC=3, TRI=0.12, HALF=15.
// Pileup is built from ARCS's own read→consensus placements (CallData); k-mers are
// counted internally. Read-placement confidence uses MAPQ=60-6*mm (mm = mismatches
// vs contig); the mapq<20 filter of the reference pipeline == skip reads with mm>=7.

namespace {

constexpr int HDMAX = 2, MC = 3, HALF = 15;
constexpr double MAF = 0.20, DHI = 2.5, KHI = 1.1, TRI = 0.12;

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

    // ── 1. Internal canonical-31-mer counts + coverage histogram (H = homozygous peak) ──
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
    std::unordered_map<uint32_t, uint64_t> hist;
    for (auto& kv : kc) ++hist[kv.second];
    uint32_t H = 30;
    { uint64_t bestc = 0; for (auto& hv : hist) if (hv.first >= 5 && hv.second > bestc) { bestc = hv.second; H = hv.first; } }
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

    // median candidate depth
    std::vector<int> depths; depths.reserve(C.size());
    for (auto& kv : C) depths.push_back(kv.second.d);
    std::sort(depths.begin(), depths.end());
    int med = depths.empty() ? 30 : depths[depths.size() / 2];

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
        if (c.d > DHI * med) continue;
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
        if ((double)c.cnt[o[2]] > TRI * c.d) continue;
        int mcnt = 0;
        if (mn2 != FLmin.end()) { std::string s = most_common(mn2->second, mcnt); (void)s; }
        if (mcnt < MC) continue;
        kept.push_back({(uint32_t)(key >> 32), (uint32_t)(key & 0xFFFFFFFFu)});
    }

    // ── 6. VCF output (contig coordinates; reference-free) ──
    std::sort(kept.begin(), kept.end());
    FILE* f = fopen(out_vcf.c_str(), "w");
    if (!f) { fprintf(stderr, "caller: cannot open %s\n", out_vcf.c_str()); return -1; }
    fprintf(f, "##fileformat=VCFv4.2\n##source=ARCS-reffree-caller\n");
    fprintf(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (auto& kp : kept) {
        uint32_t cid = kp.first, pos = kp.second;
        const std::string& cc = cd.contigs[cid];
        char ref = pos < cc.size() ? cc[pos] : 'N';
        // ALT = the observed allele (M or mn) that differs from the consensus REF.
        Cand& c = C[colkey(cid, pos)];
        int refb = b2i(ref);
        int altb = (c.M != refb) ? c.M : c.mn;
        char alt = "ACGT"[altb & 3];
        fprintf(f, "contig_%u\t%u\t.\t%c\t%c\t.\tPASS\tAF=%.3f;DP=%d\n",
                cid, pos + 1u, ref, alt, (double)c.cnt[c.mn] / c.d, c.d);
    }
    fclose(f);
    fprintf(stderr, "[CALL] H=%u candidates=%zu calls=%zu -> %s\n",
            H, C.size(), kept.size(), out_vcf.c_str());
    return (int)kept.size();
}
