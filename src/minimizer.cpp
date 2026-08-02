#include "minimizer.h"
#include <deque>

// ── extract_minimizers ─────────────────────────────────────────────────────────
// Sliding window minimum using a monotonic deque (O(n) total).
std::vector<Minimizer> extract_minimizers(const std::string& seq, int w, int k) {
    ARCS_CHECK(w >= 1 && k >= 1, "w and k must be >= 1");
    if ((int)seq.size() < k) return {};

    KmerEncoder enc(k);
    std::vector<Minimizer> result;
    result.reserve(seq.size() / w + 1);

    // Pre-compute all k-mers
    struct KmerPos { kmer_t km; bool rc; };
    std::vector<KmerPos> kmers;
    kmers.reserve(seq.size() - k + 1);

    kmer_t fwd = 0, rev = 0;
    int valid = 0;
    kmer_t rc_mask = (k < 32) ? ((kmer_t(1) << (2*k)) - 1) : ~kmer_t(0);

    for (int i = 0; i < (int)seq.size(); ++i) {
        uint8_t b = encode_base(seq[i]);
        if (b >= 4) {
            fwd = rev = 0; valid = 0;
            kmers.push_back({KMER_INVALID, false});
            continue;
        }
        // Forward strand
        fwd = enc.slide(fwd, b);
        // Reverse complement: new base goes at MSB
        rev = ((kmer_t(COMPLEMENT[b]) << (2*(k-1))) | (rev >> 2)) & rc_mask;
        ++valid;

        if (valid >= k) {
            bool is_rc = rev < fwd;
            kmers.push_back({std::min(fwd, rev), is_rc});
        } else {
            kmers.push_back({KMER_INVALID, false});
        }
    }

    // Sliding window minimum over w k-mers (positions [i, i+w-1])
    // Deque stores indices of candidate minimizers in increasing order of value
    struct Cand { kmer_t km; int pos; };
    std::deque<Cand> deq;

    kmer_t prev_min = KMER_INVALID;
    int    prev_pos = -1;

    auto emit = [&](kmer_t km, int pos, bool rc) {
        if (km != KMER_INVALID && (km != prev_min || pos != prev_pos)) {
            result.push_back({km, (pos_t)pos, rc});
            prev_min = km;
            prev_pos = pos;
        }
    };

    for (int i = 0; i < (int)kmers.size(); ++i) {
        kmer_t km = kmers[i].km;

        // Remove elements outside window
        while (!deq.empty() && deq.front().pos < i - w + 1)
            deq.pop_front();

        // Remove elements larger than current (maintain monotonic increasing)
        if (km != KMER_INVALID) {
            while (!deq.empty() && deq.back().km >= km)
                deq.pop_back();
            deq.push_back({km, i});
        }

        // Window is complete when i >= w - 1 (using k-mer indices, not base indices)
        if (i >= w - 1 && !deq.empty()) {
            int pos_in_kmers = deq.front().pos;
            emit(deq.front().km, pos_in_kmers, kmers[pos_in_kmers].rc);
        }
    }

    return result;
}

// ── MinimizerIndex::build ─────────────────────────────────────────────────────
void MinimizerIndex::build(const std::string& genome, int w, int k) {
    w_ = w; k_ = k;
    genome_len_ = genome.size();
    index_.clear();
    n_entries_ = 0;

    auto mins = extract_minimizers(genome, w, k);
    for (const auto& m : mins) {
        if (m.value != KMER_INVALID) {
            index_[m.value].push_back(m.pos);
            ++n_entries_;
        }
    }

    // Sort each bucket for binary search
    for (auto& [km, positions] : index_)
        std::sort(positions.begin(), positions.end());
}

const std::vector<pos_t>* MinimizerIndex::lookup(kmer_t minimizer) const {
    auto it = index_.find(minimizer);
    return (it != index_.end()) ? &it->second : nullptr;
}
