#pragma once
#include "minimizer.h"
#include "common.h"
#include <vector>
#include <string>

// ── Seed chain ────────────────────────────────────────────────────────────────
// A chain is a sorted sequence of matching seeds (query_pos, ref_pos) pairs.
// Seeds in a chain are collinear (same orientation, monotonically increasing).
struct Seed {
    int query_pos;  // position in read (0-based)
    int ref_pos;    // position in pseudogenome
};

struct Chain {
    std::vector<Seed> seeds;
    int   score   = 0;    // number of seeds × weight
    bool  rc      = false;
};

// ── Edit-distance alignment ────────────────────────────────────────────────────
// Returns number of mismatches between query[0..qlen) and ref[rpos..rpos+qlen)
// Also fills mm_offsets (within-query positions) and mm_bases (ref bases there).
// Uses banded DP if max_edit > 0, otherwise exact.
int align_edit(
    const std::string& query,
    const std::string& ref,
    int                ref_pos,
    int                max_edit,
    std::vector<uint32_t>& mm_offsets,
    std::vector<uint8_t>&  mm_bases
);

// ── Read mapper ────────────────────────────────────────────────────────────────
// Maps each read to the pseudogenome using minimizer seeds + chain scoring
// + edit distance extension.
//
// Mapping rate target: >85% on standard Illumina WGS data.

class ReadMapper {
public:
    ReadMapper(const std::string& genome,
               const MinimizerIndex& index,
               const ARCSParams& params);

    // Map a single read. Sets MapResult::mapped=true if mapping found.
    MapResult map(const Read& r) const;

    // Map a batch, returns results in order.
    std::vector<MapResult> map_batch(const std::vector<Read>& reads) const;

    // For PE: map R2 using R1 mapping to narrow search window
    // (insert_mean, insert_std: estimated from first mapped pairs)
    MapResult map_pair_r2(const Read& r2,
                          const MapResult& r1_result,
                          int insert_mean,
                          int insert_std) const;

    const std::string& genome() const { return genome_; }

private:
    const std::string&    genome_;
    const MinimizerIndex& index_;
    const ARCSParams&     params_;
    int max_edit_;

    // Collect candidate positions from minimizer index
    std::vector<Seed> collect_seeds(const std::string& seq, bool rc) const;

    // Chain seeds: keep collinear seeds, score = count
    std::vector<Chain> chain_seeds(const std::vector<Seed>& seeds, bool rc) const;

    // Score chains and pick best
    const Chain* best_chain(const std::vector<Chain>& chains) const;

    // Try to align read at ref position rpos (full edit distance)
    MapResult try_align(const std::string& seq, int rpos, bool rc) const;
};
