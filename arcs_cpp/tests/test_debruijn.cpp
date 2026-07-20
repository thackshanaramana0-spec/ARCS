#include "../src/debruijn.h"
#include <cassert>
#include <cstdio>

int main() {
    // Test with a known tiny genome
    // Genome: "ACGTACGT" (8 bases)
    // 4-mers: ACGT, CGTA, GTAC, TACG, ACGT (wraps)
    std::string genome = "ACGTACGTTT";
    KmerEncoder enc(4);
    std::vector<Read> reads;
    // Simulate reads from genome
    for (int i = 0; i <= (int)genome.size() - 4; ++i)
        reads.push_back({"r" + std::to_string(i), genome.substr(i, 4), "IIII"});

    KmerCounter cnt = count_kmers(reads, enc, true);
    assert(cnt.size() > 0);

    DeBruijnGraph graph;
    graph.build(cnt, enc, 1);
    assert(graph.n_nodes() > 0);
    assert(graph.n_edges() > 0);

    // Test contig extraction
    auto contigs = graph.extract_contigs(enc);
    assert(!contigs.empty());

    // Test that all contigs are at least k bases long
    for (const auto& c : contigs)
        assert((int)c.size() >= enc.k());

    printf("PASS: test_debruijn (nodes=%zu, edges=%zu, contigs=%zu)\n",
           graph.n_nodes(), graph.n_edges(), contigs.size());
    return 0;
}
