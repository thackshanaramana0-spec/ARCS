#pragma once
#include "common.h"
#include <unordered_map>
#include <vector>
#include <cmath>
#include <functional>

// ── K-mer encoding ────────────────────────────────────────────────────────────
// Encodes k-mers as uint64_t (2 bits per base, A=00, C=01, G=10, T=11).
// Maximum k = 31 (62 bits). Canonical form: min(fwd, revcomp).

class KmerEncoder {
public:
    explicit KmerEncoder(int k) : k_(k) {
        ARCS_CHECK(k >= 1 && k <= 31, "k must be in [1,31]");
        mask_ = (k < 32) ? ((kmer_t(1) << (2*k)) - 1) : ~kmer_t(0);
    }

    int k() const { return k_; }
    kmer_t mask() const { return mask_; }

    // Encode a single base. Returns KMER_INVALID if base is 'N' or invalid.
    kmer_t encode_base(char c) const {
        uint8_t b = ::encode_base(c);
        return (b < 4) ? kmer_t(b) : KMER_INVALID;
    }

    // Canonical k-mer: min(forward, reverse_complement)
    kmer_t canonical(kmer_t fwd) const {
        return std::min(fwd, revcomp(fwd));
    }

    // Reverse complement of a k-mer
    kmer_t revcomp(kmer_t km) const {
        kmer_t rc = 0;
        kmer_t x  = km;
        for (int i = 0; i < k_; ++i) {
            rc = (rc << 2) | (kmer_t(3) - (x & 3)); // complement of 2-bit base
            x >>= 2;
        }
        return rc;
    }

    // Extract all valid canonical k-mers from a sequence.
    // Skips k-mers containing 'N'.
    void extract(const std::string& seq,
                 std::vector<kmer_t>& out,
                 bool canonical_form = true) const;

    // Slide-extend: given current k-mer km and a new base b (0-3),
    // return new k-mer by shifting left and appending b.
    kmer_t slide(kmer_t km, uint8_t b) const {
        return ((km << 2) | b) & mask_;
    }

    // Decode k-mer back to string (for debugging)
    std::string decode(kmer_t km) const {
        std::string s(k_, 'N');
        for (int i = k_-1; i >= 0; --i) {
            s[i] = BASE_TO_CHAR[km & 3];
            km >>= 2;
        }
        return s;
    }

private:
    int    k_;
    kmer_t mask_;
};

// ── K-mer counter (robin-hood open-addressing hash table) ────────────────────
// Stores kmer_t → count_t. Uses linear probing with tombstone deletion.
// Load factor target: 70% before resize.

class KmerCounter {
public:
    explicit KmerCounter(size_t capacity_hint = 1 << 20);

    // Increment count for k-mer. Returns new count.
    count_t insert(kmer_t km);

    // Get count (0 if absent)
    count_t get(kmer_t km) const;

    // Iterate all (k-mer, count) pairs
    void foreach(const std::function<void(kmer_t, count_t)>& fn) const;

    // Remove all entries with count < threshold
    void filter_below(count_t threshold);

    size_t size()     const { return size_; }
    size_t capacity() const { return cap_; }

    // Coverage statistics from k-mer histogram
    struct CovStats {
        double mean_cov;   // estimated mean coverage
        count_t threshold; // error threshold = max(2, floor(sqrt(mean_cov)))
        size_t n_solid;    // k-mers above threshold
        size_t n_error;    // k-mers below threshold
    };
    CovStats coverage_stats() const;

private:
    static constexpr kmer_t EMPTY    = ~kmer_t(0);
    static constexpr kmer_t DELETED  = ~kmer_t(0) - 1;
    static constexpr float  MAX_LOAD = 0.70f;

    struct Slot { kmer_t key = EMPTY; count_t val = 0; };

    std::vector<Slot> table_;
    size_t size_ = 0;
    size_t cap_  = 0;
    size_t dels_ = 0;

    size_t probe(kmer_t km) const {
        // Fibonacci hashing for better distribution
        return (size_t)(km * 11400714819323198485ULL) & (cap_ - 1);
    }

    void rehash(size_t new_cap);
    void maybe_grow();
};

// ── Error corrector ────────────────────────────────────────────────────────────
// Uses k-mer frequency to correct low-frequency (likely erroneous) k-mers.
// Strategy: for each k-mer in a read with count < threshold, try all single
// substitutions; accept the one with the highest count (must be >= threshold).
// Only corrects if exactly one high-frequency neighbor exists (unambiguous).

class ErrorCorrector {
public:
    ErrorCorrector(const KmerCounter& counts,
                   const KmerEncoder& enc,
                   count_t threshold)
        : counts_(counts), enc_(enc), threshold_(threshold) {}

    // Correct a single read sequence in-place.
    // Returns number of bases corrected.
    int correct(std::string& seq) const;

    // Correct a batch of reads.
    size_t correct_batch(std::vector<Read>& reads) const;

private:
    const KmerCounter& counts_;
    const KmerEncoder& enc_;
    count_t threshold_;

    // Try to correct one k-mer at position i in seq.
    // Returns true if correction applied.
    bool correct_at(std::string& seq, int i) const;
};

// ── Auto k-selection ───────────────────────────────────────────────────────────
// Select k such that a random k-mer is genome-unique with high probability.
// Formula: P(unique) ≈ 1 - genome_size/4^k
// We want P(unique) > 0.99 → k > log4(100 × genome_size)
// genome_size is estimated from the k-mer histogram valley.
int select_k(size_t n_reads, size_t read_len,
             int k_min = 15, int k_max = 31);

// Build k-mer counter from a set of reads (canonical k-mers)
KmerCounter count_kmers(const std::vector<Read>& reads,
                        const KmerEncoder& enc,
                        bool canonical = true);
