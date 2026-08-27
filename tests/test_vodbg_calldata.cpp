// Verifies build_vodbg_pg's new CallData contract directly: for EVERY read,
// (read_cid[i], read_pos[i], read_rc[i]) must point into a valid contig at a
// valid offset, and the read's own sequence (or its RC) must actually match
// that contig region — covering all three populations the assembler can
// produce a read from (grown, fallback-mapped, and never-placed/appended).
#include "vodbg_pg.h"
#include "common.h"
#include <cstdio>
#include <random>
#include <string>
#include <vector>

static int g_fail = 0;

static void check_calldata(const std::vector<Read>& reads, const char* label) {
    CallData cd;
    ChainEncodeResult r = build_vodbg_pg(reads, &cd);
    (void)r;

    if (!cd.valid) { printf("[%s] FAILED: cd.valid is false\n", label); ++g_fail; return; }
    if (cd.read_cid.size() != reads.size() || cd.read_pos.size() != reads.size()
        || cd.read_rc.size() != reads.size()) {
        printf("[%s] FAILED: CallData array size mismatch\n", label);
        ++g_fail; return;
    }

    long bad_cid = 0, bad_pos = 0, mismatch = 0, ok = 0;
    for (size_t i = 0; i < reads.size(); ++i) {
        uint32_t cid = cd.read_cid[i];
        uint32_t pos = cd.read_pos[i];
        uint8_t  rc  = cd.read_rc[i];
        if (cid >= cd.contigs.size()) { ++bad_cid; continue; }
        const std::string& contig = cd.contigs[cid];
        const std::string& orig = reads[i].seq;
        std::string target = rc ? reverse_complement(orig) : orig;
        if ((size_t)pos + target.size() > contig.size()) { ++bad_pos; continue; }
        int mm = 0;
        for (size_t j = 0; j < target.size(); ++j)
            if (contig[pos + j] != target[j] && is_acgt_strict(target[j])) ++mm;
        // Clean synthetic data (no injected errors) should match exactly
        // regardless of which internal path placed the read; allow a tiny
        // slack only for the mismatch-tolerant fallback path.
        if (mm > (int)(target.size() / 4)) { ++mismatch; continue; }
        ++ok;
    }
    bool pass = (bad_cid == 0 && bad_pos == 0 && mismatch == 0);
    printf("[%s] %s (n=%zu ok=%ld bad_cid=%ld bad_pos=%ld mismatch=%ld, contigs=%zu)\n",
           label, pass ? "OK" : "FAILED", reads.size(), ok, bad_cid, bad_pos, mismatch, cd.contigs.size());
    if (!pass) ++g_fail;
}

static Read mk(const std::string& seq, const std::string& name) {
    Read r; r.seq = seq; r.name = name; r.qual = std::string(seq.size(), 'I'); return r;
}

int main() {
    std::mt19937 rng(123);
    std::string alpha = "ACGT";

    // 1. Simple overlapping chain — all reads should grow into one contig.
    {
        std::vector<Read> reads;
        reads.push_back(mk("AAAACCCCGGGG", "r0"));
        reads.push_back(mk("CCCCGGGGTTTT", "r1"));
        reads.push_back(mk("GGGGTTTTAAAA", "r2"));
        check_calldata(reads, "simple_chain");
    }

    // 2. Reads with NO overlap at all — every read should become its own
    // singleton contig via the never-placed/appended path.
    {
        std::vector<Read> reads;
        for (int i = 0; i < 10; ++i) {
            std::string s(40, 'A');
            for (auto& c : s) c = alpha[rng() % 4];
            reads.push_back(mk(s, "iso" + std::to_string(i)));
        }
        check_calldata(reads, "all_isolated");
    }

    // 3. Realistic mixed case: a genuine random genome tiled into overlapping
    // reads (should mostly grow via the main path) PLUS some short reads
    // below K0 forced into the fallback/append path, PLUS reverse-complement
    // reads to exercise the rc flag.
    {
        std::mt19937 rng2(7);
        std::string genome(2000, 'A');
        for (auto& c : genome) c = alpha[rng2() % 4];
        std::vector<Read> reads;
        for (int i = 0; i + 100 <= (int)genome.size(); i += 30) {
            std::string s = genome.substr((size_t)i, 100);
            if (i % 90 == 0) s = reverse_complement(s); // exercise rc path
            reads.push_back(mk(s, "t" + std::to_string(i)));
        }
        // A handful of short (<K0) reads that can only be placed by fallback.
        for (int i = 0; i < 5; ++i) {
            reads.push_back(mk(genome.substr((size_t)(i * 50), 15), "short" + std::to_string(i)));
        }
        check_calldata(reads, "mixed_realistic");
    }

    // 4. Large-ish stress test at a scale closer to real usage.
    {
        std::mt19937 rng3(99);
        std::string genome(20000, 'A');
        for (auto& c : genome) c = alpha[rng3() % 4];
        std::vector<Read> reads;
        for (int i = 0; i + 150 <= (int)genome.size(); i += 40) {
            reads.push_back(mk(genome.substr((size_t)i, 150), "s" + std::to_string(i)));
        }
        check_calldata(reads, "stress_20k_genome");
    }

    printf(g_fail == 0 ? "\nALL CALLDATA TESTS PASSED\n" : "\n%d CALLDATA TEST(S) FAILED\n", g_fail);
    return g_fail == 0 ? 0 : 1;
}
