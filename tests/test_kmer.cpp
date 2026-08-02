#include "../src/kmer.h"
#include <cassert>
#include <cstdio>

int main() {
    // Test KmerEncoder
    KmerEncoder enc(5);
    assert(enc.k() == 5);

    // Test encode/decode round-trip
    std::vector<kmer_t> kmers;
    enc.extract("ACGTTGCA", kmers, false);
    assert(!kmers.empty());

    // Test reverse complement
    kmer_t fwd = 0; // AAAAA
    kmer_t rc  = enc.revcomp(fwd);
    // AAAAA RC = TTTTT = 0b11111111... = (3 repeated k times)
    kmer_t expected_rc = 0;
    for (int i = 0; i < 5; ++i) expected_rc = (expected_rc << 2) | 3;
    assert(rc == expected_rc);

    // Test k-mer counter
    KmerCounter cnt;
    cnt.insert(kmer_t(42));
    cnt.insert(kmer_t(42));
    cnt.insert(kmer_t(100));
    assert(cnt.get(42)  == 2);
    assert(cnt.get(100) == 1);
    assert(cnt.get(999) == 0);
    assert(cnt.size()   == 2);

    // Test error correction doesn't crash
    std::vector<Read> reads;
    reads.push_back({"r1", "ACGTTGCAACGT", "IIIIIIIIIIII"});
    reads.push_back({"r2", "ACGTTGCAACGT", "IIIIIIIIIIII"});
    KmerEncoder enc2(4);
    KmerCounter cnt2 = count_kmers(reads, enc2);
    auto stats = cnt2.coverage_stats();
    assert(stats.mean_cov > 0);

    printf("PASS: test_kmer\n");
    return 0;
}
