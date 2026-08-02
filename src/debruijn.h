#pragma once
#include "kmer.h"
#include <unordered_map>
#include <vector>
#include <string>
#include <functional>

// ── De Bruijn Graph ───────────────────────────────────────────────────────────
// Nodes = (k-1)-mers  (encoded as kmer_t with k-1 bases)
// Edges = k-mers     (encoded as kmer_t with k bases)
//
// The graph is stored as adjacency lists: for each node (k-1)-mer,
// we store its out-edges (the base appended to form a k-mer).
//
// Canonical edges: we store k-mers in canonical form to halve storage,
// but track orientation separately for Eulerian path construction.
//
// Design: separate construction (KmerCounter input) from traversal (Eulerian path).

class DeBruijnGraph {
public:
    struct EdgeInfo {
        uint8_t base;    // base extending this node (0-3)
        bool    rc;      // edge is from reverse-complement direction
        count_t count;   // k-mer count (multiplicity)
    };

    // Build graph from k-mer counts.
    // Only solid k-mers (count >= min_count) are included.
    void build(const KmerCounter& counts,
               const KmerEncoder& enc,
               count_t min_count = 1);

    // Remove tips: dead-end branches shorter than max_tip_len bases.
    void remove_tips(int max_tip_len);

    // Extract contigs: maximal paths with no branching (in_deg=out_deg=1).
    // Returns list of contig strings.
    std::vector<std::string> extract_contigs(const KmerEncoder& enc) const;

    // Eulerian path (Hierholzer's algorithm).
    // Returns the pseudogenome string G* as a single sequence.
    // Concatenates all linear contigs; for real Eulerian circuits,
    // decomposes the graph into paths and joins them.
    std::string eulerian_pseudogenome(const KmerEncoder& enc) const;

    // Stats
    size_t n_nodes() const { return adj_.size(); }
    size_t n_edges() const { return n_edges_; }
    int    k()       const { return k_; }

    // In/out degrees for Eulerian condition check
    int out_degree(kmer_t node) const;
    int in_degree(kmer_t node) const;

private:
    struct NodeData {
        std::vector<EdgeInfo> out;     // out-edges
        int                   in_deg = 0;
    };

    std::unordered_map<kmer_t, NodeData> adj_;  // node → adjacency
    size_t n_edges_ = 0;
    int    k_       = 0;

    // Hierholzer implementation
    std::vector<kmer_t> hierholzer(kmer_t start) const;

    // Convert node path to sequence string
    std::string path_to_string(const std::vector<kmer_t>& path,
                                const KmerEncoder& enc) const;

    // Find all connected components (for multi-path genomes)
    std::vector<kmer_t> find_start_nodes() const;
};

// ── Full pipeline: reads → pseudogenome string ────────────────────────────────
struct PseudogenomeResult {
    std::string           genome;      // G* string
    int                   k;           // k used
    size_t                n_kmers;     // solid k-mers in graph
    double                est_genome_size;
    std::string           method;      // "debruijn" or "greedy_scs"
    // Exact positions from SCS construction (empty if method=="debruijn")
    // If populated: read i maps to genome[positions[i]..positions[i]+read_len]
    // with 0 mismatches — no re-mapping needed.
    std::vector<pos_t>    scs_positions;
};

PseudogenomeResult build_pseudogenome(
    const std::vector<Read>& reads,
    const ARCSParams& params
);

// ── Greedy SCS builder (fallback for low coverage) ────────────────────────────
// Chains reads by k-length suffix/prefix overlap using a hash table.
// Works at any coverage. Returns both the pseudogenome string AND the
// exact position of each read in the pseudogenome (0 mismatches guaranteed
// for all chained reads).
struct SCSResult {
    std::string          genome;
    std::vector<pos_t>   positions;   // genome position of read i
    std::vector<bool>    rc;          // reverse complement flag (always false for SCS)
    std::vector<int>     n_mm;        // mismatch count (0 for all SCS reads)
};

SCSResult build_greedy_scs_with_positions(const std::vector<std::string>& seqs, int k);
std::string build_greedy_scs(const std::vector<std::string>& seqs, int k); // legacy
