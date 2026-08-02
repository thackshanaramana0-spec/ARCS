#pragma once
#include "common.h"
#include <vector>
#include <string>

// ── Hamming-Shift MST for unmapped reads ─────────────────────────────────────
// Based on: "Hamming-shifting graph of genomic short reads" (PLOS Comp Bio, 2021)
// Computes approximate Minimum Spanning Tree of reads by Hamming+shift distance.
// Reads ordered by DFS traversal of MST achieve 10-30% better BWT compressibility
// than random ordering (measured in the 2021 paper).
//
// For n unmapped reads, full pairwise distance is O(n²×L) — infeasible.
// We approximate using LSH (Locality Sensitive Hashing) to find k-nearest neighbors.
// Then run Kruskal's MST on the k-NN graph.

// ── Hamming+shift distance ────────────────────────────────────────────────────
// Distance(r1, r2) = min over shifts s in [-max_shift, +max_shift] of
//                   Hamming(r1, r2[s:s+len]) / len
// Normalized to [0, 1] for MST edge weights.
float hamming_shift_dist(const std::string& a, const std::string& b,
                         int max_shift = 5);

// ── LSH for Hamming distance ──────────────────────────────────────────────────
// Random bit-sampling: hash a read to an L-bit integer by sampling L random
// positions from its 2-bit encoding. Two reads at Hamming distance d have
// collision probability (1 - d/len)^L.
// Use multiple (n_bands) hash functions for k-NN approximation.
class HammingLSH {
public:
    // seed: random seed for reproducibility
    // n_hashes: number of independent hash functions
    // bits_per_hash: bits per hash (affects recall/precision tradeoff)
    HammingLSH(int read_len, int n_hashes, int bits_per_hash, uint64_t seed = 42);

    // Compute all hash values for a read sequence.
    // Returns vector of n_hashes uint64_t values.
    std::vector<uint64_t> hash(const std::string& seq) const;

    // Check if two reads are candidates (share >= 1 hash bucket)
    bool is_candidate(const std::vector<uint64_t>& h1,
                      const std::vector<uint64_t>& h2) const;

private:
    int read_len_;
    int n_hashes_;
    int bits_;
    // positions_[i][j] = which sequence position to sample for hash i, bit j
    std::vector<std::vector<int>> positions_;
};

// ── k-NN graph ────────────────────────────────────────────────────────────────
struct KNNEdge {
    uint32_t u, v;   // read indices
    float    dist;   // Hamming-shift distance
};

// Build k-nearest-neighbor graph for unmapped reads using LSH.
// Returns edges sorted by distance (ascending).
std::vector<KNNEdge> build_knn_graph(
    const std::vector<std::string>& seqs,
    int k_neighbors,
    int n_lsh_hashes,
    int max_shift = 5,
    uint64_t seed = 42
);

// ── Kruskal's MST ─────────────────────────────────────────────────────────────
// Union-Find for Kruskal's algorithm.
class UnionFind {
public:
    explicit UnionFind(size_t n);
    int  find(int x);
    bool unite(int x, int y);
private:
    std::vector<int> parent_, rank_;
};

// Run Kruskal's MST on edges, return MST edge list.
std::vector<KNNEdge> kruskal_mst(std::vector<KNNEdge> edges, size_t n_nodes);

// ── DFS order from MST ────────────────────────────────────────────────────────
// Returns indices of reads in DFS traversal order of the MST.
// Reads visited in this order have high pairwise similarity → good BWT compression.
std::vector<uint32_t> mst_dfs_order(
    const std::vector<KNNEdge>& mst,
    size_t n_nodes
);

// ── Full pipeline ─────────────────────────────────────────────────────────────
// Input: sequences of unmapped reads.
// Output: reordering indices (result[i] = index in input of the i-th read).
std::vector<uint32_t> order_unmapped_reads(
    const std::vector<std::string>& seqs,
    const ARCSParams& params
);
