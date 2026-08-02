#include "mapper.h"
#include <algorithm>
#include <cstring>
#include <limits>

// ── Edit distance (banded DP) ─────────────────────────────────────────────────
// Aligns query against ref[ref_pos .. ref_pos + qlen - 1].
// Returns n_mismatches. Also fills mm_offsets and mm_bases.
// Uses exact match (no indels) for speed: FASTQ reads have mostly substitutions.
// Indel support can be added as an extension.
int align_edit(const std::string& query,
               const std::string& ref,
               int ref_pos,
               int max_edit,
               std::vector<uint32_t>& mm_offsets,
               std::vector<uint8_t>&  mm_bases) {
    mm_offsets.clear();
    mm_bases.clear();

    int qlen = (int)query.size();
    if (ref_pos < 0 || ref_pos + qlen > (int)ref.size())
        return max_edit + 1; // out of bounds

    int n_mm = 0;
    for (int i = 0; i < qlen; ++i) {
        char q = query[i];
        char r = ref[ref_pos + i];
        if (q != r) {
            if (max_edit >= 0 && n_mm > max_edit) return n_mm;
            mm_offsets.push_back((uint32_t)i);
            uint8_t rb = encode_base(r);
            mm_bases.push_back(rb < 4 ? rb : BASE_N);
            ++n_mm;
        }
    }
    return n_mm;
}

// ── ReadMapper ────────────────────────────────────────────────────────────────
ReadMapper::ReadMapper(const std::string& genome,
                       const MinimizerIndex& index,
                       const ARCSParams& params)
    : genome_(genome), index_(index), params_(params) {
    max_edit_ = (params.max_edit < 0)
                ? std::max(1, (int)(genome.size() > 0 ?
                    params_.k_max / 10 : 5))
                : params.max_edit;
}

std::vector<Seed> ReadMapper::collect_seeds(const std::string& seq, bool rc) const {
    std::string query = rc ? reverse_complement(seq) : seq;
    auto mins = extract_minimizers(query, params_.w, params_.mink);

    std::vector<Seed> seeds;
    seeds.reserve(mins.size() * 3); // average bucket size

    for (const auto& m : mins) {
        if (m.value == KMER_INVALID) continue;
        const auto* positions = index_.lookup(m.value);
        if (!positions) continue;

        for (pos_t ref_pos : *positions) {
            // ref_pos is position of minimizer in genome.
            // query_pos is position of minimizer in query (m.pos).
            // Reference alignment start = ref_pos - m.pos
            int aln_start = (int)ref_pos - (int)m.pos;
            if (aln_start < 0 || aln_start + (int)query.size() > (int)genome_.size())
                continue;
            seeds.push_back({(int)m.pos, aln_start});
        }
    }
    return seeds;
}

std::vector<Chain> ReadMapper::chain_seeds(const std::vector<Seed>& seeds, bool rc) const {
    if (seeds.empty()) return {};

    // Group seeds by their implied alignment start position (ref_pos - query_pos)
    // Seeds with the same alignment start are in the same chain.
    std::unordered_map<int, Chain> chains_map;

    for (const auto& s : seeds) {
        int aln = s.ref_pos; // we already computed aln_start = ref_pos
        auto& chain = chains_map[aln];
        chain.seeds.push_back(s);
        chain.score++;
        chain.rc = rc;
    }

    // Convert to vector and sort by score
    std::vector<Chain> chains;
    chains.reserve(chains_map.size());
    for (auto& [k, c] : chains_map)
        chains.push_back(std::move(c));

    std::sort(chains.begin(), chains.end(),
              [](const Chain& a, const Chain& b) { return a.score > b.score; });

    // Return top-K chains
    static constexpr int TOP_K = 10;
    if ((int)chains.size() > TOP_K) chains.resize(TOP_K);
    return chains;
}

const Chain* ReadMapper::best_chain(const std::vector<Chain>& chains) const {
    if (chains.empty()) return nullptr;
    return &chains[0]; // already sorted by score
}

MapResult ReadMapper::try_align(const std::string& seq, int rpos, bool rc) const {
    std::vector<uint32_t> mm_offsets;
    std::vector<uint8_t>  mm_bases;

    int n_mm = align_edit(seq, genome_, rpos, max_edit_, mm_offsets, mm_bases);

    MapResult res;
    res.pos    = (pos_t)rpos;
    res.mapped = (n_mm <= max_edit_);
    res.rc     = rc;
    res.n_mm   = n_mm;
    if (res.mapped) {
        res.mm_offsets = std::move(mm_offsets);
        res.mm_bases   = std::move(mm_bases);
    }
    return res;
}

MapResult ReadMapper::map(const Read& r) const {
    int max_edit_for_read = std::max(1, (int)r.seq.size() / 10);

    // Try forward strand
    for (bool rc : {false, true}) {
        std::string query = rc ? reverse_complement(r.seq) : r.seq;
        auto seeds  = collect_seeds(r.seq, rc);
        auto chains = chain_seeds(seeds, rc);
        const Chain* best = best_chain(chains);
        if (!best) continue;

        // Try top chains in order
        for (const auto& chain : chains) {
            int aln_start = chain.seeds.empty() ? -1 : chain.seeds[0].ref_pos;
            if (aln_start < 0) continue;

            MapResult res = try_align(query, aln_start, rc);
            if (res.mapped && res.n_mm <= max_edit_for_read) {
                return res;
            }
        }
    }

    MapResult unmapped;
    unmapped.mapped = false;
    return unmapped;
}

std::vector<MapResult> ReadMapper::map_batch(const std::vector<Read>& reads) const {
    std::vector<MapResult> results(reads.size());
    for (size_t i = 0; i < reads.size(); ++i)
        results[i] = map(reads[i]);
    return results;
}

MapResult ReadMapper::map_pair_r2(const Read& r2,
                                   const MapResult& r1_res,
                                   int insert_mean,
                                   int insert_std) const {
    if (!r1_res.mapped) return map(r2); // fallback to full search

    // Expected R2 position: R1_pos + insert_mean (R2 maps opposite strand)
    int expected = (int)r1_res.pos + insert_mean;
    int window   = insert_std * 3; // ±3σ covers 99.7% of inserts

    // Try positions in the window
    int rlen = (int)r2.seq.size();
    std::string rc_seq = reverse_complement(r2.seq); // R2 is typically on RC strand

    int best_mm = max_edit_ + 1;
    MapResult best;
    best.mapped = false;

    for (int delta = -window; delta <= window; delta += 5) {
        int rpos = expected + delta;
        if (rpos < 0 || rpos + rlen > (int)genome_.size()) continue;

        std::vector<uint32_t> mm_off;
        std::vector<uint8_t>  mm_bas;
        int n_mm = align_edit(rc_seq, genome_, rpos, max_edit_, mm_off, mm_bas);
        if (n_mm < best_mm) {
            best_mm     = n_mm;
            best.pos    = (pos_t)rpos;
            best.mapped = (n_mm <= max_edit_);
            best.rc     = true;
            best.n_mm   = n_mm;
            best.mm_offsets = mm_off;
            best.mm_bases   = mm_bas;
        }
    }

    if (!best.mapped) return map(r2); // fallback
    return best;
}
