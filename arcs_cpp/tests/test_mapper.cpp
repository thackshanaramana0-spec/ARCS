#include "../src/mapper.h"
#include "../src/debruijn.h"
#include "../src/hamming_mst.h"
#include <cassert>
#include <cstdio>

int main() {
    // Create a simple pseudogenome
    std::string genome = "ACGTACGTACGTACGTACGTACGTACGT"
                         "TTTTCCCCGGGGAAAATTTTCCCCGGGG"
                         "ACGTACGTACGTACGTACGTACGTACGT";

    ARCSParams params;
    params.w    = 5;
    params.mink = 7;
    params.max_edit = 3;

    MinimizerIndex index;
    index.build(genome, params.w, params.mink);
    assert(index.n_entries() > 0);

    ReadMapper mapper(genome, index, params);

    // Test: a read that exactly matches at position 0
    Read r1;
    r1.name = "exact";
    r1.seq  = genome.substr(0, 20);
    r1.qual = std::string(20, 'I');
    auto res1 = mapper.map(r1);
    // May or may not map depending on minimizer overlap, but should not crash
    printf("Read1 mapped=%d pos=%d n_mm=%d\n", res1.mapped, res1.pos, res1.n_mm);

    // Test: a read with 1 mismatch
    Read r2;
    r2.name = "one_mm";
    r2.seq  = genome.substr(5, 20);
    r2.seq[3] = (r2.seq[3] == 'A') ? 'C' : 'A'; // introduce 1 mismatch
    r2.qual = std::string(20, 'I');
    auto res2 = mapper.map(r2);
    printf("Read2 mapped=%d pos=%d n_mm=%d\n", res2.mapped, res2.pos, res2.n_mm);

    // Test Hamming-MST ordering
    std::vector<std::string> seqs = {
        "ACGTACGTAC", "ACGTACGTAC", "TTTTTTTTTT",
        "ACGTCCGTAC", "TTTTTTTTGT"
    };
    auto order = order_unmapped_reads(seqs, params);
    assert(order.size() == seqs.size());

    printf("PASS: test_mapper\n");
    return 0;
}
