#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <array>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cmath>

// ── Nucleotide encoding ───────────────────────────────────────────────────────
// A=0, C=1, G=2, T=3, N=4  (2-bit for ACGT, 3-bit for ACGTN)
static constexpr uint8_t BASE_A = 0;
static constexpr uint8_t BASE_C = 1;
static constexpr uint8_t BASE_G = 2;
static constexpr uint8_t BASE_T = 3;
static constexpr uint8_t BASE_N = 4;
static constexpr uint8_t BASE_INVALID = 255;

// ASCII → 2-bit encoding table (size 256)
inline const std::array<uint8_t,256>& base_table() {
    static std::array<uint8_t,256> tbl = [](){
        std::array<uint8_t,256> t;
        t.fill(BASE_INVALID);
        t['A'] = t['a'] = BASE_A;
        t['C'] = t['c'] = BASE_C;
        t['G'] = t['g'] = BASE_G;
        t['T'] = t['t'] = BASE_T;
        t['U'] = t['u'] = BASE_T; // RNA
        t['N'] = t['n'] = BASE_N;
        return t;
    }();
    return tbl;
}

// 2-bit → ASCII
static constexpr char BASE_TO_CHAR[5] = {'A','C','G','T','N'};

// Complement: A↔T, C↔G, N↔N
static constexpr uint8_t COMPLEMENT[5] = {BASE_T, BASE_G, BASE_C, BASE_A, BASE_N};

// ── k-mer types ───────────────────────────────────────────────────────────────
using kmer_t = uint64_t;   // up to k=31 (62 bits used)
using count_t = uint32_t;
using pos_t   = uint32_t;  // position in pseudogenome
using readid_t = uint32_t;

static constexpr kmer_t KMER_INVALID = ~kmer_t(0);

// ── Compression parameters (runtime-configurable) ─────────────────────────────
struct ARCSParams {
    int    k        = 0;    // k-mer size (0 = auto-select)
    int    k_min    = 15;
    int    k_max    = 31;
    int    w        = 10;   // minimizer window
    int    mink     = 15;   // minimizer k
    int    max_edit = -1;   // -1 = auto (read_len / 10)
    int    knn      = 40;   // Hamming-MST k-nearest-neighbors
    int    n_lsh    = 8;    // LSH hash functions
    int    pilot    = 5000; // reads for model training
    float  akc_amplicon  = 0.15f;
    float  akc_wgs       = 0.35f;
    bool   paired_end    = false;
    bool   store_names   = true;
    bool   gutensor      = false;
    bool   use_mst       = true;  // use MST tree encoder (beats SCS by 22-50%)
    bool   use_mst_dfs   = false; // MST with DFS-position parent encoding (better at low coverage)
    bool   use_chain     = false; // use greedy k-NN chain encoder (trial11 — beats MST+SPRING)
    bool   use_chain_pg  = false; // chain-pg: build pseudogenome from chain walk, compress with GeCo3
    int    threads       = 1;
    int    quality_bins  = 0;     // 0=lossless, 4=Novaseq 4-bin, 8=uniform 8-bin
    int    chunk_size    = 0;     // 0=load all reads (default); >0=chunked MST mode
};

// ── Quality binning (lossy pre-processing) ────────────────────────────────────
// Maps raw Phred [0,42] → a small set of representative values.
// Applied before compression; decoder outputs binned values as-is (no inverse).
inline char bin_quality_char(char qchar, int n_bins) {
    int q = (unsigned char)qchar - 33; // ASCII Phred+33 → raw score
    if (q < 0) q = 0; if (q > 42) q = 42;
    int binned;
    if (n_bins == 4) {
        // Novaseq-like 4-bin: breakpoints at 3,15,27
        if      (q <  3) binned = 2;
        else if (q < 15) binned = 12;
        else if (q < 27) binned = 23;
        else             binned = 37;
    } else if (n_bins == 2) {
        binned = (q < 20) ? 10 : 35;
    } else {
        // Uniform n_bins bands over [0,42]
        int band = q * n_bins / 43;
        if (band >= n_bins) band = n_bins - 1;
        binned = (band * 42) / (n_bins - 1);
    }
    return (char)(binned + 33);
}

// ── Read record ───────────────────────────────────────────────────────────────
struct Read {
    std::string name;
    std::string seq;
    std::string qual;
    std::string plus;                         // '+' line content AFTER the '+' (empty = bare "+")
    bool        reverse_complemented = false; // set after mapping
};

// ── Mapping result ────────────────────────────────────────────────────────────
struct MapResult {
    pos_t  pos    = 0;
    bool   mapped = false;
    bool   rc     = false;
    int    n_mm   = 0;                  // mismatch count
    std::vector<uint32_t> mm_offsets;  // within-read offsets
    std::vector<uint8_t>  mm_bases;    // original bases at mismatch positions
};

// ── AKC result ────────────────────────────────────────────────────────────────
enum class DataRegime : uint8_t {
    AMPLICON    = 0,
    TARGETED    = 1,
    WGS         = 2,
    METAGENOMIC = 3,
};

struct AKCResult {
    float      score  = 0.0f;
    DataRegime regime = DataRegime::WGS;
    bool       use_bsc = false;  // true for AMPLICON/TARGETED
};

// ── Error handling ────────────────────────────────────────────────────────────
class ARCSError : public std::runtime_error {
public:
    explicit ARCSError(const std::string& msg) : std::runtime_error(msg) {}
};

#define ARCS_CHECK(cond, msg) \
    do { if(!(cond)) throw ARCSError(std::string(__FILE__) + ":" + std::to_string(__LINE__) + " " + (msg)); } while(0)

// ── Utility ───────────────────────────────────────────────────────────────────
inline uint8_t encode_base(char c) {
    return base_table()[static_cast<uint8_t>(c)];
}

inline std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.size(), 'N');
    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t b = encode_base(seq[seq.size() - 1 - i]);
        rc[i] = (b <= BASE_T) ? BASE_TO_CHAR[COMPLEMENT[b]] : 'N';
    }
    return rc;
}

inline double log2_safe(double x) {
    return (x > 0.0) ? std::log2(x) : 0.0;
}

// CRC32 (Castagnoli) — for container integrity
uint32_t crc32c(const uint8_t* data, size_t len, uint32_t crc = 0xFFFFFFFFu);
