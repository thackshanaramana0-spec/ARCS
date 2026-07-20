#pragma once
#include "kmer.h"
#include <unordered_map>
#include <vector>

// ── (w,k)-Minimizer ────────────────────────────────────────────────────────────
// For a window of w consecutive k-mers, the minimizer is the lexicographically
// smallest canonical k-mer. This gives a sparse index of the genome.
// Roberts et al. 2004: (w,k)-minimizers have expected density 2/(w+1).

struct Minimizer {
    kmer_t  value;  // minimizer k-mer (canonical)
    pos_t   pos;    // position in source string (0-based, leftmost base of window)
    bool    rc;     // true if minimizer was on reverse strand
};

// Extract all minimizers from a string (genome or read).
std::vector<Minimizer> extract_minimizers(const std::string& seq,
                                          int w, int k);

// ── Minimizer index ────────────────────────────────────────────────────────────
// Maps minimizer k-mer → sorted list of positions in the pseudogenome.
// Used by the mapper to seed read alignment.

class MinimizerIndex {
public:
    // Build index from pseudogenome.
    void build(const std::string& genome, int w, int k);

    // Look up positions for a minimizer.
    // Returns empty span if not found.
    const std::vector<pos_t>* lookup(kmer_t minimizer) const;

    // Parameters
    int w() const { return w_; }
    int k() const { return k_; }
    size_t genome_len() const { return genome_len_; }
    size_t n_entries()  const { return n_entries_; }

private:
    std::unordered_map<kmer_t, std::vector<pos_t>> index_;
    int    w_          = 0;
    int    k_          = 0;
    size_t genome_len_ = 0;
    size_t n_entries_  = 0;
};
