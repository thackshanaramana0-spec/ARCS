#include "kmer.h"
#include <cmath>
#include <numeric>
#include <functional>

// ── KmerEncoder ───────────────────────────────────────────────────────────────
void KmerEncoder::extract(const std::string& seq,
                          std::vector<kmer_t>& out,
                          bool canonical_form) const {
    if ((int)seq.size() < k_) return;
    out.reserve(out.size() + seq.size() - k_ + 1);

    kmer_t fwd = 0;
    int    valid = 0; // consecutive valid (non-N) bases

    for (int i = 0; i < (int)seq.size(); ++i) {
        uint8_t b = ::encode_base(seq[i]);
        if (b >= 4) {        // N or invalid
            fwd   = 0;
            valid = 0;
            continue;
        }
        fwd   = slide(fwd, b);
        valid++;
        if (valid >= k_) {
            kmer_t km = canonical_form ? canonical(fwd) : fwd;
            out.push_back(km);
        }
    }
}

// ── KmerCounter ───────────────────────────────────────────────────────────────
KmerCounter::KmerCounter(size_t capacity_hint) {
    // Round up to next power of 2
    cap_ = 1;
    while (cap_ < capacity_hint * 2) cap_ <<= 1;
    table_.assign(cap_, {EMPTY, 0});
}

count_t KmerCounter::insert(kmer_t km) {
    maybe_grow();
    size_t h = probe(km);
    size_t first_del = SIZE_MAX;

    while (true) {
        Slot& s = table_[h];
        if (s.key == EMPTY) {
            size_t idx = (first_del != SIZE_MAX) ? first_del : h;
            if (first_del != SIZE_MAX) --dels_;
            table_[idx] = {km, 1};
            ++size_;
            return 1;
        }
        if (s.key == DELETED) {
            if (first_del == SIZE_MAX) first_del = h;
        } else if (s.key == km) {
            return ++s.val;
        }
        h = (h + 1) & (cap_ - 1);
    }
}

count_t KmerCounter::get(kmer_t km) const {
    size_t h = probe(km);
    while (true) {
        const Slot& s = table_[h];
        if (s.key == EMPTY)   return 0;
        if (s.key == km)      return s.val;
        h = (h + 1) & (cap_ - 1);
    }
}

void KmerCounter::foreach(const std::function<void(kmer_t, count_t)>& fn) const {
    for (size_t i = 0; i < cap_; ++i) {
        const Slot& s = table_[i];
        if (s.key != EMPTY && s.key != DELETED)
            fn(s.key, s.val);
    }
}

void KmerCounter::filter_below(count_t threshold) {
    for (auto& s : table_) {
        if (s.key != EMPTY && s.key != DELETED && s.val < threshold) {
            s.key = DELETED;
            --size_;
            ++dels_;
        }
    }
}

void KmerCounter::rehash(size_t new_cap) {
    std::vector<Slot> old = std::move(table_);
    cap_  = new_cap;
    size_ = 0;
    dels_ = 0;
    table_.assign(cap_, {EMPTY, 0});
    for (const auto& s : old)
        if (s.key != EMPTY && s.key != DELETED) {
            // Re-insert without calling maybe_grow
            size_t h = probe(s.key);
            while (table_[h].key != EMPTY && table_[h].key != DELETED)
                h = (h + 1) & (cap_ - 1);
            table_[h] = s;
            ++size_;
        }
}

void KmerCounter::maybe_grow() {
    if ((float)(size_ + dels_) / cap_ > MAX_LOAD)
        rehash(cap_ * 2);
}

KmerCounter::CovStats KmerCounter::coverage_stats() const {
    // Build histogram of counts (cap at 1000 to avoid huge arrays)
    static constexpr count_t MAX_HIST = 1000;
    std::vector<size_t> hist(MAX_HIST + 1, 0);
    size_t total_solid = 0;
    double total_count = 0.0;

    auto hist_fn = [&](kmer_t, count_t c) {
        hist[std::min(c, MAX_HIST)]++;
        total_count += c;
    };
    foreach(hist_fn);

    // Find the valley between error peak (count=1-3) and true k-mer peak
    // Simple heuristic: threshold at sqrt(peak_count)
    size_t peak_val = 1;
    size_t peak_freq = 0;
    for (size_t c = 2; c <= MAX_HIST; ++c) {
        if (hist[c] > peak_freq) { peak_freq = hist[c]; peak_val = c; }
    }

    double mean_cov = (size_ > 0) ? total_count / size_ : 1.0;
    count_t thr = std::max(count_t(2), (count_t)std::floor(std::sqrt(mean_cov)));

    for (size_t c = thr; c <= MAX_HIST; ++c) total_solid += hist[c];

    return { mean_cov, thr, total_solid, size_ - total_solid };
}

// ── ErrorCorrector ────────────────────────────────────────────────────────────
bool ErrorCorrector::correct_at(std::string& seq, int pos) const {
    // Build the k-mer containing position pos (last base of k-mer starting at pos-k+1)
    // Try all 3 alternative bases at position pos.
    char original = seq[pos];
    uint8_t orig_b = ::encode_base(original);
    if (orig_b >= 4) return false; // N cannot be corrected

    kmer_t best_km   = KMER_INVALID;
    count_t best_cnt = 0;
    char   best_char = original;
    int    n_good    = 0;

    // We need to find all k-mers covering position pos and check each substitution.
    // For efficiency, we re-extract one k-mer at a time around pos.
    int k = enc_.k();
    int start = std::max(0, pos - k + 1);
    int end   = std::min((int)seq.size() - k, pos);

    for (int s = start; s <= end; ++s) {
        // Extract k-mer starting at s
        kmer_t fwd = 0;
        bool valid = true;
        for (int j = s; j < s + k; ++j) {
            uint8_t b = ::encode_base(seq[j]);
            if (b >= 4) { valid = false; break; }
            fwd = enc_.slide(fwd, b);
        }
        if (!valid) continue;

        // Check if this k-mer is solid
        kmer_t can = enc_.canonical(fwd);
        if (counts_.get(can) >= threshold_) continue; // already solid, no need

        // Try all 3 substitutions at local position (pos - s)
        int local = pos - s;
        uint8_t orig = ::encode_base(seq[pos]);

        for (uint8_t alt = 0; alt < 4; ++alt) {
            if (alt == orig) continue;
            // Rebuild k-mer with substitution at local position
            kmer_t alt_fwd = 0;
            for (int j = s; j < s + k; ++j) {
                uint8_t b = (j == pos) ? alt : ::encode_base(seq[j]);
                if (b >= 4) { alt_fwd = KMER_INVALID; break; }
                alt_fwd = enc_.slide(alt_fwd, b);
            }
            if (alt_fwd == KMER_INVALID) continue;
            count_t cnt = counts_.get(enc_.canonical(alt_fwd));
            if (cnt >= threshold_) {
                ++n_good;
                if (cnt > best_cnt) {
                    best_cnt  = cnt;
                    best_char = BASE_TO_CHAR[alt];
                    best_km   = alt_fwd;
                }
            }
        }
    }

    // Only correct if exactly one high-frequency neighbor (unambiguous)
    if (n_good == 1 && best_char != original) {
        seq[pos] = best_char;
        return true;
    }
    return false;
}

int ErrorCorrector::correct(std::string& seq) const {
    int corrected = 0;
    for (int i = 0; i < (int)seq.size(); ++i)
        if (correct_at(seq, i)) ++corrected;
    return corrected;
}

size_t ErrorCorrector::correct_batch(std::vector<Read>& reads) const {
    size_t total = 0;
    for (auto& r : reads)
        total += correct(r.seq);
    return total;
}

// ── Auto k-selection ───────────────────────────────────────────────────────────
int select_k(size_t n_reads, size_t read_len, int k_min, int k_max) {
    // Genome size estimate (very rough): total_bases / coverage
    // Coverage estimate: assume ~30x as default if unknown
    // Use Lander-Waterman: coverage = n_reads * read_len / genome_size
    // We don't know genome_size, so use heuristic based on read_len.
    // For k_opt: k-mer should be genome-unique with P > 0.99
    // → 4^k > 100 * genome_size_estimate
    // Conservative estimate: genome ~500 MB (mammalian), so k_opt ≈ 15 for bacteria, 19 for human
    // We return a fixed default of 21 (genome-assembler standard)

    // Better: use n_reads and read_len to estimate coverage at different k values
    // Optimal k maximises n_solid_kmers * (read_len - k + 1) / n_reads
    // Without actual counting, return standard 21 (works well empirically)

    int k = 21;
    if (read_len <= 50)  k = 15;  // short reads (M.tb 50bp)
    if (read_len >= 150) k = 21;  // standard Illumina
    if (read_len >= 250) k = 25;  // long Illumina

    return std::clamp(k, k_min, k_max);
}

// ── Build k-mer counter from reads ────────────────────────────────────────────
KmerCounter count_kmers(const std::vector<Read>& reads,
                        const KmerEncoder& enc,
                        bool canonical) {
    KmerCounter cnt(reads.size() * 50); // rough initial capacity
    std::vector<kmer_t> kmers;
    for (const auto& r : reads) {
        kmers.clear();
        enc.extract(r.seq, kmers, canonical);
        for (kmer_t km : kmers)
            cnt.insert(km);
    }
    return cnt;
}
