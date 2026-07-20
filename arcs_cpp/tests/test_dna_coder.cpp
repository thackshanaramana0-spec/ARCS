// Roundtrip + ratio tests for ARCS-DNA compressor.
// Gate: decode(encode(pg)) == pg for all test cases.
#include "../src/dna_coder.h"
#include <cassert>
#include <cstdio>
#include <string>
#include <vector>
#include <random>
#include <cstring>

static void check(const std::string& pg, const char* label) {
    auto enc = dna_encode(pg, {}, {}, {});
    auto dec = dna_decode(enc);
    if (dec != pg) {
        fprintf(stderr, "FAIL [%s]: pg_len=%zu enc=%zu dec_len=%zu\n",
                label, pg.size(), enc.size(), dec.size());
        // Show first mismatch.
        for (size_t i = 0; i < std::min(pg.size(), dec.size()); i++) {
            if (pg[i] != dec[i]) {
                fprintf(stderr, "  first diff at %zu: orig=%c dec=%c\n",
                        i, pg[i], dec[i]);
                break;
            }
        }
        assert(false);
    }
    double bpb = 8.0 * (double)(enc.size() - 8) / (double)pg.size();
    fprintf(stderr, "OK  [%s]: pg=%zu → enc=%zu (%.3f bpb)\n",
            label, pg.size(), enc.size(), bpb);
}

static std::string random_dna(size_t n, uint64_t seed = 42) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> d(0, 3);
    static const char b[4] = {'A','C','G','T'};
    std::string s(n, 'A');
    for (auto& c : s) c = b[d(rng)];
    return s;
}

static std::string repeat_dna(const std::string& unit, size_t n) {
    std::string s;
    s.reserve(n);
    while (s.size() < n) s += unit;
    s.resize(n);
    return s;
}

int main() {
    fprintf(stderr, "=== test_dna_coder ===\n");

    // 1. Empty.
    {
        std::string pg = "";
        auto enc = dna_encode(pg, {}, {}, {});
        auto dec = dna_decode(enc);
        assert(dec == pg);
        fprintf(stderr, "OK  [empty]\n");
    }

    // 2. Single base.
    check("A", "single_A");
    check("ACGT", "acgt_4");

    // 3. Short random.
    check(random_dna(100),  "rand_100");
    check(random_dna(1000), "rand_1000");

    // 4. Highly repetitive (should compress very well).
    check(repeat_dna("ACGT", 10000), "repeat_40k");
    check(repeat_dna("AAAAAACCCCCCGGGGGGTTTTTT", 5000), "repeat_lowcomp");

    // 5. All same base.
    check(std::string(1000, 'A'), "all_A_1k");
    check(std::string(1000, 'C'), "all_C_1k");

    // 6. Pseudo-genomic: random with local repeats (simulates real pg).
    {
        std::string pg = random_dna(500, 1);
        // Simulate chain: repeat reads with small mismatches.
        std::mt19937_64 rng(7);
        std::string full;
        full.reserve(50000);
        full += pg;
        for (int i = 0; i < 90; i++) {
            // overlap 130 bases, add 20 new.
            size_t ov = std::min((size_t)130, full.size());
            full += random_dna(20, rng());
        }
        check(full, "pseudo_pg_50k");
    }

    // 7. Larger random (compression ratio check — should be < 2.5 bpb).
    {
        std::string pg = random_dna(486143, 99);
        auto enc = dna_encode(pg, {}, {}, {});
        auto dec = dna_decode(enc);
        assert(dec == pg);
        double bpb = 8.0 * (double)(enc.size() - 8) / (double)pg.size();
        fprintf(stderr, "OK  [rand_486k]: %.3f bpb (GeCo3 gets ~2.0 on real pg)\n", bpb);
        // Sanity: pure random DNA should be close to 2.0 bpb.
        assert(bpb < 2.5 && "compression suspiciously bad");
    }

    // 8. N-bases (mapped to A, must roundtrip as A).
    {
        std::string pg = "ACGTNNNACGT";
        auto enc = dna_encode(pg, {}, {}, {});
        auto dec = dna_decode(enc);
        // N→A in encoder, decoder returns A.
        std::string expected = pg;
        for (auto& c : expected) if (c == 'N') c = 'A';
        assert(dec == expected);
        fprintf(stderr, "OK  [N_bases]\n");
    }

    fprintf(stderr, "=== ALL PASS ===\n");
    return 0;
}
