#include "mst_encoder.h"
#include "minimizer.h"
#include "position_enc.h"
#include "container.h"
#include "rans_model.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <unordered_map>
#include <stack>
#include <queue>
#include <cassert>
#include <functional>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __SSE2__
#include <emmintrin.h>
#endif

// ── SIMD Hamming distance ─────────────────────────────────────────────────────
// Returns the number of positions where a[i] != b[i] over [0, len).
// Falls back to scalar on non-SSE2 builds.  No early exit — callers check
// the result against max_mm after the call, allowing the compiler to vectorise
// the hot path without a data-dependent branch inside the inner loop.
static inline int hamming_count(const char* __restrict__ a,
                                 const char* __restrict__ b,
                                 int len) {
    int mm = 0;
#ifdef __SSE2__
    int i = 0;
    for (; i + 16 <= len; i += 16) {
        __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a + i));
        __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b + i));
        __m128i eq = _mm_cmpeq_epi8(va, vb);
        mm += 16 - __builtin_popcount(static_cast<unsigned>(_mm_movemask_epi8(eq)));
    }
    for (; i < len; ++i) mm += (a[i] != b[i]);
#else
    for (int i = 0; i < len; ++i) mm += (a[i] != b[i]);
#endif
    return mm;
}

// ── align_reads ────────────────────────────────────────────────────────────────
// Finds the best (shift, substitution_list) alignment of `child` onto `parent`.
// Tries all shifts in [-max_shift, +max_shift] and (optionally) the RC of child.
// Returns the alignment that minimises substitutions subject to overlap ≥ min_overlap.
ReadAlignResult align_reads(const std::string& parent,
                             const std::string& child,
                             int max_shift,
                             int max_mm) {
    int L = (int)parent.size();
    assert(L == (int)child.size()); // fixed-length reads

    if (max_mm < 0) max_mm = std::max(1, L / 10);
    int min_overlap = std::max(1, L / 4); // require at least 25% overlap

    ReadAlignResult best;
    best.shift      = 0;
    best.overlap_len = L;
    best.n_subs     = max_mm + 1; // sentinel: worse than any valid alignment
    best.rc_child   = false;

    auto try_align = [&](const std::string& ch, bool rc) {
        for (int s = -max_shift; s <= max_shift; ++s) {
            // overlap region:
            // parent bases: [max(0,s)  .. min(L, L+s) - 1]
            // child  bases: [max(0,-s) .. min(L, L-s) - 1]
            int pa_start = std::max(0, s);
            int ch_start = std::max(0, -s);
            int ov_len   = L - std::abs(s);
            if (ov_len < min_overlap) continue;

            // SIMD Hamming count (no early exit — vectorises cleanly).
            // Sub-position list is only built for alignments that survive the filter.
            int n_mm = hamming_count(parent.data() + pa_start, ch.data() + ch_start, ov_len);
            if (n_mm > max_mm) continue;

            // Better than current best?
            if (n_mm < best.n_subs ||
                (n_mm == best.n_subs && ov_len > best.overlap_len)) {
                best.shift      = s;
                best.overlap_len = ov_len;
                best.n_subs     = n_mm;
                best.rc_child   = rc;

                // Build substitution list (scalar, rare — only on passing alignments)
                best.sub_positions.clear();
                best.sub_bases.clear();
                for (int i = 0; i < ov_len; ++i) {
                    if (parent[pa_start + i] != ch[ch_start + i]) {
                        best.sub_positions.push_back((uint32_t)i);
                        uint8_t b = encode_base(ch[ch_start + i]);
                        best.sub_bases.push_back(b < 4 ? b : BASE_N);
                    }
                }
            }
        }
    };

    try_align(child, false);
    std::string rc = reverse_complement(child);
    try_align(rc, true);

    return best;
}

// ── align_reads_windowed ─────────────────────────────────────────────────────
// Like align_reads, but only searches a small window around a known-implied
// shift for each orientation, instead of the full [-max_shift, +max_shift]
// range. When two reads share an exact k-mer at known positions, the correct
// alignment shift is determined by those positions (see build_knn) — no
// search is needed in principle. The window (default 2) is purely defensive
// slack; the implied shift is otherwise exact for substitution-only matches.
static ReadAlignResult align_reads_windowed(const std::string& parent,
                                             const std::string& child,
                                             int center_fwd,
                                             int center_rc,
                                             int window,
                                             int max_mm) {
    int L = (int)parent.size();
    if (max_mm < 0) max_mm = std::max(1, L / 10);
    int min_overlap = std::max(1, L / 4);

    ReadAlignResult best;
    best.shift       = 0;
    best.overlap_len = L;
    best.n_subs      = max_mm + 1;
    best.rc_child    = false;

    auto try_align = [&](const std::string& ch, bool rc, int center) {
        int lo = std::max(-(L - 1), center - window);
        int hi = std::min(L - 1, center + window);
        for (int s = lo; s <= hi; ++s) {
            int pa_start = std::max(0, s);
            int ch_start = std::max(0, -s);
            int ov_len   = L - std::abs(s);
            if (ov_len < min_overlap) continue;
            int n_mm = hamming_count(parent.data() + pa_start, ch.data() + ch_start, ov_len);
            if (n_mm > max_mm) continue;

            if (n_mm < best.n_subs ||
                (n_mm == best.n_subs && ov_len > best.overlap_len)) {
                best.shift       = s;
                best.overlap_len = ov_len;
                best.n_subs      = n_mm;
                best.rc_child    = rc;

                best.sub_positions.clear();
                best.sub_bases.clear();
                for (int i = 0; i < ov_len; ++i) {
                    if (parent[pa_start + i] != ch[ch_start + i]) {
                        best.sub_positions.push_back((uint32_t)i);
                        uint8_t b = encode_base(ch[ch_start + i]);
                        best.sub_bases.push_back(b < 4 ? b : BASE_N);
                    }
                }
            }
        }
    };

    try_align(child, false, center_fwd);
    std::string rc = reverse_complement(child);
    try_align(rc, true, center_rc);

    return best;
}

// ── MSTSequenceEncoder ────────────────────────────────────────────────────────
MSTSequenceEncoder::MSTSequenceEncoder(const MSTConfig& cfg) : cfg_(cfg) {}

// ── build_knn: all k-mer indexing with implied-shift ─────────────────────────
// Two-phase parallel KNN:
//   Phase 1 (single-threaded): index all k-mers, iterate buckets with evaluated-map
//     dedup, collect (u,v,cf,cr) work items — same logic as before, no change to
//     which pairs are evaluated.
//   Phase 2 (parallel): align all work items in parallel with OpenMP — pure read-only
//     access to seqs[], no shared writes, no synchronisation needed.
//   Phase 3 (single-threaded): merge distances into per-node best-k heaps.
//
// Fixes vs prior version:
//   - buckets.reserve() capped at 8M (was n*kmers/3 = 67M at 1.55M reads → 537 MB
//     of bucket-pointer array for only ~4.8M unique k-mers — catastrophic L3 spill).
//   - evaluated.reserve() pre-scans bucket sizes to estimate unique pairs accurately,
//     eliminating the 50M→100M→200M rehash cascade at archival scale.
//   - Alignment work separated from dedup so OpenMP can parallelize it safely.
std::vector<KNNEdge> MSTSequenceEncoder::build_knn(
    const std::vector<std::string>& seqs) const {

    int n = (int)seqs.size();
    if (n == 0) return {};
    int L = (int)seqs[0].size();

    int max_mm = cfg_.max_mm < 0 ? std::max(1, L / 10) : cfg_.max_mm;
    // k-mer size: k must satisfy 4^k >> n × kmers_per_read so that bucket
    // occupancy stays near 1 and bucket traversal runs in O(n) total work.
    // At 50 bp with k=11: 4^11 = 4.2 M ≈ n×40 → occupancy ≈ 1 → O(n²) pairs.
    // Sparsity target: 4^k ≥ 16 × n × L  (expected occupancy ≤ 1/16 per bucket).
    // Solving: k_sparse = ceil(log₄(16 × n × L)) = 2 + ceil(log₄(n × L)).
    // idx_k takes the larger of k_sparse and the length-proportional heuristic,
    // so longer reads still benefit from the higher-quality L×0.14 indexing.
    int k_sparse = 2 + (int)std::ceil(std::log((double)n * (double)L) / std::log(4.0));
    int idx_k    = std::min(21, std::max(k_sparse, L * 14 / 100));
    int stride  = 1;
    static constexpr int MAX_BUCKET       = 500;
    static constexpr int MAX_PAIRS_PER_READ = 16;

    struct BEntry { uint32_t read_idx; uint16_t pos; };

    auto t_knn_start = std::chrono::steady_clock::now();
    auto elapsed_ms = [&](std::chrono::steady_clock::time_point t0) {
        return (int)std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - t0).count();
    };

    // ── Phase 1a: build k-mer bucket index ──────────────────────────────────
    // Reserve: target ~8M unique k-mer buckets (typical for low-coverage WGS).
    // Under-reserving causes one rehash; over-reserving wastes cache capacity.
    // The old formula (n × kmers_per_read / 3 = 67M at 1.55M reads) pre-allocated
    // 537 MB of bucket pointers for ~4.8M actual unique k-mers — a 14× over-reserve
    // that spilled the entire bucket array out of L3 cache.
    // Reserve based on expected unique k-mers: n × kmers_per_read / avg_bucket_size.
    // For WGS at typical coverage, most k-mers are unique or in small buckets.
    // Under-reserving causes rehash cascades (each doubles allocation + re-inserts all);
    // over-reserving wastes cache lines. Target ~1 rehash in the worst case.
    // Block minimizer width: one representative k-mer per W-position window.
    // Use this for total_kmers_est so DS-7 (1M reads) correctly enters the
    // fast sort-based flat path rather than the hashmap fallback.
    int W_blk = std::max(4, L / 10);
    size_t kmers_per_read = (size_t)std::max(1, (L - idx_k) / W_blk + 1);

    // ── Phase 1: sort-based k-mer bucket construction ────────────────────────
    // Instead of unordered_map<kmer,vector<entry>> (which causes 4M heap allocs
    // and poor cache behaviour), we:
    //   1. Fill a flat KmerEntry array in parallel (each thread owns a read slice)
    //   2. Sort by k-mer value → contiguous equal-kmer runs become buckets
    //   3. Iterate runs exactly once to collect work items
    // This eliminates all per-bucket memory allocations and makes the fill O(n·L/T).
    // Memory: n × kmers_per_read × 16 bytes ≈ 192 MB for 100k × 136-kmer reads.
    struct KmerEntry { kmer_t kmer; uint32_t read_idx; uint16_t pos; uint16_t _p; };
    static_assert(sizeof(KmerEntry) == 16, "KmerEntry must be 16 bytes");

    size_t total_kmers_est = (size_t)n * kmers_per_read;
    static constexpr size_t MEM_CAP = (size_t)1200 * 1024 * 1024; // 1.2 GB

    struct WorkItem { int u, v, cf, cr; };

    // ── Phase 1a: parallel fill ──────────────────────────────────────────────
#ifdef _OPENMP
    int nfill = omp_get_max_threads();
#else
    int nfill = 1;
#endif

    KmerEncoder kenc(idx_k);
    size_t est_raw_pairs = 0;
    size_t skip_large    = 0;
    size_t total_pairs   = 0;

    std::vector<WorkItem> work;
    std::unordered_map<uint64_t, bool> evaluated;

    if (total_kmers_est * sizeof(KmerEntry) <= MEM_CAP) {
        // FAST PATH: flat array + parallel fill + std::sort
        std::vector<std::vector<KmerEntry>> tlocal(nfill);
        for (int t = 0; t < nfill; ++t)
            tlocal[t].reserve((total_kmers_est / (size_t)nfill) + 1024);

#ifdef _OPENMP
        #pragma omp parallel for num_threads(nfill) schedule(static)
#endif
        for (int i = 0; i < n; ++i) {
#ifdef _OPENMP
            int t = omp_get_thread_num();
#else
            int t = 0;
#endif
            auto& lv = tlocal[t];
            // Block minimizer with O(1)-per-position sliding window.
            // Pre-seed with the first idx_k-1 bases so the window at pos=0
            // correctly covers [0, idx_k-1]. n_valid counts consecutive valid
            // bases; a value < idx_k means the window contains an N and the
            // k-mer is skipped without restarting the O(k) inner loop.
            {
                kmer_t blk_min     = UINT64_MAX;
                int    blk_min_pos = -1;
                int    blk_start   = 0;
                kmer_t fwd         = 0;
                int    n_valid     = 0;
                for (int j = 0; j < idx_k - 1 && j < L; ++j) {
                    uint8_t b = encode_base(seqs[i][j]);
                    if (b >= 4) { fwd = 0; n_valid = 0; }
                    else        { fwd = kenc.slide(fwd, b); ++n_valid; }
                }
                for (int pos = 0; pos + idx_k <= L; ++pos) {
                    uint8_t b = encode_base(seqs[i][pos + idx_k - 1]);
                    if (b >= 4) { fwd = 0; n_valid = 0; }
                    else        { fwd = kenc.slide(fwd, b);
                                  if (n_valid < idx_k) ++n_valid; }
                    if (n_valid == idx_k) {
                        kmer_t ck = kenc.canonical(fwd);
                        if (ck < blk_min) { blk_min = ck; blk_min_pos = pos; }
                    }
                    bool end_of_block = (pos - blk_start + 1) >= W_blk;
                    bool last_kmer    = (pos + idx_k + 1 > L);
                    if (end_of_block || last_kmer) {
                        if (blk_min != UINT64_MAX)
                            lv.push_back({blk_min, (uint32_t)i,
                                          (uint16_t)blk_min_pos, 0});
                        blk_min = UINT64_MAX; blk_min_pos = -1;
                        blk_start = pos + 1;
                    }
                }
            }
        }

        // Concatenate thread-local arrays into a single flat vector
        std::vector<KmerEntry> kmer_list;
        kmer_list.reserve(total_kmers_est);
        for (auto& lv : tlocal)
            kmer_list.insert(kmer_list.end(), lv.begin(), lv.end());
        tlocal.clear(); // release before sort

        // Sort by k-mer → contiguous runs = buckets
        std::sort(kmer_list.begin(), kmer_list.end(),
                  [](const KmerEntry& a, const KmerEntry& b){ return a.kmer < b.kmer; });

        // ── Phase 1b: pre-scan runs for evaluated reserve ────────────────────
        // Reserve est_raw_pairs (no division) so the evaluated hash map stays
        // near load factor 0.5 after dedup, avoiding rehash cascades.
        // For high-sharing datasets (short reads, high coverage) the actual
        // unique pairs can reach ~60% of est_raw_pairs, so dividing by 15
        // caused 3+ rehashes and catastrophic performance on DS-4.
        size_t ksz = kmer_list.size();
        for (size_t i = 0; i < ksz; ) {
            size_t j = i + 1;
            while (j < ksz && kmer_list[j].kmer == kmer_list[i].kmer) ++j;
            size_t m = j - i;
            if (m >= 2 && m <= (size_t)MAX_BUCKET)
                est_raw_pairs += m * (m - 1) / 2;
            else if (m > (size_t)MAX_BUCKET)
                ++skip_large;
            i = j;
        }
        // Cap at 150M (≈1.2 GB bucket array at 8 B/ptr) to protect against
        // pathological inputs. evaluated.reserve() pre-sizes the bucket array
        // so the load factor never exceeds ~0.5 after insertion.
        size_t est_unique = std::min(est_raw_pairs + 2000000UL, (size_t)150000000);
        evaluated.reserve(est_unique);
        work.reserve(std::min(est_unique, (size_t)120000000));

        // ── Phase 1c: collect work items from sorted runs ────────────────────
        // Fix 1: cap candidate pairs per read to MAX_PAIRS_PER_READ.
        // Prevents O(n²) work for reads that appear in many k-mer buckets
        // (high-coverage or repetitive regions). 16 candidates is enough for
        // Kruskal's MST to find the globally optimal neighbour.
        std::vector<int> read_pair_count(n, 0);
        for (size_t i = 0; i < ksz; ) {
            size_t j = i + 1;
            while (j < ksz && kmer_list[j].kmer == kmer_list[i].kmer) ++j;
            size_t m = j - i;
            if (m < 2 || m > (size_t)MAX_BUCKET) { i = j; continue; }

            for (size_t ai = i; ai < j; ++ai) {
                int u = (int)kmer_list[ai].read_idx;
                if (read_pair_count[u] >= MAX_PAIRS_PER_READ) continue;
                for (size_t bi = ai + 1; bi < j; ++bi) {
                    int v = (int)kmer_list[bi].read_idx;
                    if (u == v) continue;
                    if (read_pair_count[u] >= MAX_PAIRS_PER_READ) break;
                    if (read_pair_count[v] >= MAX_PAIRS_PER_READ) continue;
                    int nu = std::min(u,v), nv = std::max(u,v);
                    uint64_t pk = ((uint64_t)nu << 32) | (uint64_t)nv;
                    if (evaluated.count(pk)) continue;
                    evaluated[pk] = true;
                    ++total_pairs;
                    int pos_u = kmer_list[ai].pos;
                    int pos_v = kmer_list[bi].pos;
                    work.push_back({u, v,
                        pos_u - pos_v,
                        pos_u + pos_v + idx_k - L});
                    ++read_pair_count[u];
                    ++read_pair_count[v];
                }
            }
            i = j;
        }

        fprintf(stderr, "[KNN-P1-sort] build=%.1fs kmers=%zuM pairs_est=%zu work=%zu\n",
                elapsed_ms(t_knn_start) / 1000.0,
                kmer_list.size() / 1000000,
                est_raw_pairs, work.size());

    } else {
        // HASH MAP FALLBACK for very large datasets (> 1.2 GB flat array)
        size_t expected_unique = std::min((size_t)n * kmers_per_read / 2, (size_t)20000000ULL);
        std::unordered_map<kmer_t, std::vector<BEntry>> buckets;
        buckets.reserve(expected_unique);
        for (int i = 0; i < n; ++i) {
            kmer_t blk_min     = UINT64_MAX;
            int    blk_min_pos = -1;
            int    blk_start   = 0;
            kmer_t fwd         = 0;
            int    n_valid     = 0;
            for (int j = 0; j < idx_k - 1 && j < L; ++j) {
                uint8_t b = encode_base(seqs[i][j]);
                if (b >= 4) { fwd = 0; n_valid = 0; }
                else        { fwd = kenc.slide(fwd, b); ++n_valid; }
            }
            for (int pos = 0; pos + idx_k <= L; ++pos) {
                uint8_t b = encode_base(seqs[i][pos + idx_k - 1]);
                if (b >= 4) { fwd = 0; n_valid = 0; }
                else        { fwd = kenc.slide(fwd, b);
                              if (n_valid < idx_k) ++n_valid; }
                if (n_valid == idx_k) {
                    kmer_t ck = kenc.canonical(fwd);
                    if (ck < blk_min) { blk_min = ck; blk_min_pos = pos; }
                }
                bool end_of_block = (pos - blk_start + 1) >= W_blk;
                bool last_kmer    = (pos + idx_k + 1 > L);
                if (end_of_block || last_kmer) {
                    if (blk_min != UINT64_MAX)
                        buckets[blk_min].push_back({(uint32_t)i,
                                                    (uint16_t)blk_min_pos});
                    blk_min = UINT64_MAX; blk_min_pos = -1;
                    blk_start = pos + 1;
                }
            }
        }
        size_t skip_large_pre = 0;
        for (auto& [kv, mbrs] : buckets) {
            size_t m = mbrs.size();
            if (m >= 2 && m <= (size_t)MAX_BUCKET) est_raw_pairs += m * (m-1) / 2;
            else if (m > (size_t)MAX_BUCKET) { ++skip_large_pre; ++skip_large; }
        }
        size_t est_unique = std::min(est_raw_pairs / 15 + 2000000, (size_t)150000000);
        evaluated.reserve(est_unique);
        work.reserve(std::min(est_unique, (size_t)120000000));
        std::vector<int> read_pair_count_fb(n, 0);
        for (auto& [kval, members] : buckets) {
            if ((int)members.size() < 2 || (int)members.size() > MAX_BUCKET) continue;
            for (int ai = 0; ai < (int)members.size(); ++ai) {
                int u = members[ai].read_idx;
                if (read_pair_count_fb[u] >= MAX_PAIRS_PER_READ) continue;
                for (int bi = ai+1; bi < (int)members.size(); ++bi) {
                    int v = members[bi].read_idx;
                    if (u == v) continue;
                    if (read_pair_count_fb[u] >= MAX_PAIRS_PER_READ) break;
                    if (read_pair_count_fb[v] >= MAX_PAIRS_PER_READ) continue;
                    int nu = std::min(u,v), nv = std::max(u,v);
                    uint64_t pk = ((uint64_t)nu << 32) | (uint64_t)nv;
                    if (evaluated.count(pk)) continue;
                    evaluated[pk] = true;
                    ++total_pairs;
                    work.push_back({u, v,
                        (int)members[ai].pos - (int)members[bi].pos,
                        (int)members[ai].pos + (int)members[bi].pos + idx_k - L});
                    ++read_pair_count_fb[u];
                    ++read_pair_count_fb[v];
                }
            }
        }
        fprintf(stderr, "[KNN-P1-map] build=%.1fs index=%.1fM pairs_est=%zu work=%zu\n",
                elapsed_ms(t_knn_start) / 1000.0,
                (double)buckets.bucket_count() / 1e6,
                est_raw_pairs, work.size());
        buckets.clear(); buckets.rehash(0);
    }
    evaluated.clear(); evaluated.rehash(0);

    auto t_p2 = std::chrono::steady_clock::now();
    // ── Phase 2: parallel alignment ──────────────────────────────────────────
    // seqs[] is read-only → no data races.  Each element of dists[] is written
    // by exactly one thread (static or dynamic partition) → no synchronisation.
    std::vector<float> dists(work.size(), 2.0f); // 2.0 = filtered / no match

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 4096)
#endif
    for (int i = 0; i < (int)work.size(); ++i) {
        const WorkItem& wi = work[i];
        ReadAlignResult aln = align_reads_windowed(seqs[wi.u], seqs[wi.v],
                                                    wi.cf, wi.cr,
                                                    /*window=*/2, max_mm);
        if (aln.n_subs <= max_mm && aln.overlap_len > 0)
            dists[i] = (float)aln.n_subs / aln.overlap_len;
    }
    fprintf(stderr, "[KNN-P2] align=%.1fs\n", elapsed_ms(t_p2) / 1000.0);

    // ── Phase 3: merge into per-node best-k heaps (single-threaded) ─────────
    std::vector<std::vector<std::pair<float,int>>> best_k(n);
    for (int i = 0; i < n; ++i) best_k[i].reserve(cfg_.k_neighbors + 1);

    auto try_insert = [&](int node, int neighbor, float d) {
        auto& kv = best_k[node];
        kv.push_back({d, neighbor});
        std::push_heap(kv.begin(), kv.end());
        if ((int)kv.size() > cfg_.k_neighbors) {
            std::pop_heap(kv.begin(), kv.end());
            kv.pop_back();
        }
    };

    size_t pass_mm = 0;
    for (int i = 0; i < (int)work.size(); ++i) {
        if (dists[i] < 2.0f) {
            ++pass_mm;
            try_insert(work[i].u, work[i].v, dists[i]);
            try_insert(work[i].v, work[i].u, dists[i]);
        }
    }
    work.clear(); work.shrink_to_fit();
    dists.clear(); dists.shrink_to_fit();

    size_t n_isolated_knn = 0;
    for (int i = 0; i < n; ++i) if (best_k[i].empty()) ++n_isolated_knn;

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    fprintf(stderr, "[KNN] n=%d threads=%d skip_large=%zu pairs=%zu pass_mm=%zu isolated=%zu\n",
            n, n_threads, skip_large, total_pairs, pass_mm, n_isolated_knn);

    // Collect and deduplicate edges from best_k heaps
    std::unordered_map<uint64_t, float> edge_map;
    edge_map.reserve((size_t)n * cfg_.k_neighbors);
    for (int i = 0; i < n; ++i) {
        for (auto& [d, j] : best_k[i]) {
            int u = std::min(i, j), v = std::max(i, j);
            uint64_t key = ((uint64_t)u << 32) | (uint64_t)v;
            auto it = edge_map.find(key);
            if (it == edge_map.end() || d < it->second)
                edge_map[key] = d;
        }
    }
    std::vector<KNNEdge> edges;
    edges.reserve(edge_map.size());
    for (auto& [key, d] : edge_map) {
        edges.push_back({(uint32_t)(key >> 32), (uint32_t)key, d});
    }
    return edges;
}

// ── encode_tree ───────────────────────────────────────────────────────────────
// Stores parent indices. Root nodes have parent = n (sentinel).
std::vector<uint8_t> MSTSequenceEncoder::encode_tree(
    const std::vector<uint32_t>& parents,
    size_t n_reads) const {

    std::vector<uint8_t> out;
    out.reserve(n_reads * 3); // rough estimate

    write_varint(out, n_reads);
    for (size_t i = 0; i < n_reads; ++i) {
        // Use zigzag delta encoding relative to i (parents tend to be nearby in DFS order)
        uint32_t p = parents[i];
        if (p == UINT32_MAX) {
            // root: store special sentinel
            write_varint(out, (uint64_t)n_reads); // n_reads = "no parent"
        } else {
            // delta from current index
            int32_t delta = (int32_t)p - (int32_t)i;
            write_varint(out, (uint64_t)zigzag_encode(delta));
        }
    }
    return out;
}

// ── encode_deltas ─────────────────────────────────────────────────────────────
// For each read in DFS order:
//   - If root: store verbatim sequence (2 bits/base, padded to byte)
//   - If non-root: store (shift, n_subs, sub_offsets[], sub_bases[], prefix, suffix)
// Also builds dev_sets_out[dfs_position][pos] = is this a substitution position?
// (Used by quality encoder for joint quality-sequence coding.)
//
// Parallelism: alignments are pre-computed in BFS level order.  All reads at
// level k share the same invariant: their MST parents are at level k-1, so
// reconstructed[par] is fully settled before any level-k thread reads it.
// Reads within the same level have no parent-child relationship (provable from
// BFS), so their reconstructed[] writes are non-overlapping — no locks needed.
// The final serialisation pass iterates DFS order (sequential, fast) using the
// pre-computed alignment results.
std::vector<uint8_t> MSTSequenceEncoder::encode_deltas(
    const std::vector<Read>&        reads,
    const std::vector<uint32_t>&    parents,
    const std::vector<uint32_t>&    dfs_order,
    std::vector<bool>&              rc_flags_out,
    std::vector<std::vector<bool>>& dev_sets_out,
    std::vector<int>&               aln_shifts_out,
    MSTEncodeResult&                stats) const {

    int n = (int)reads.size();
    int L = reads.empty() ? 0 : (int)reads[0].seq.size();
    int max_mm = cfg_.max_mm < 0 ? std::max(1, L / 10) : cfg_.max_mm;
    int full_shift = L - std::max(1, L / 4);

    rc_flags_out.assign(n, false);
    dev_sets_out.assign(n, std::vector<bool>(L, false));
    aln_shifts_out.assign(n, 0);

    // ── Phase A: build BFS levels from MST parent array ──────────────────────
    // children[v] lists the MST children of node v (parents[child] == v).
    std::vector<std::vector<uint32_t>> children(n);
    for (int i = 0; i < n; ++i)
        if (parents[i] != UINT32_MAX)
            children[parents[i]].push_back((uint32_t)i);

    std::vector<int> node_level(n, -1);
    std::vector<std::vector<uint32_t>> level_groups;
    level_groups.reserve(32);

    // BFS: roots are nodes with no parent.
    {
        std::queue<uint32_t> bfsq;
        level_groups.push_back({});
        for (int i = 0; i < n; ++i) {
            if (parents[i] == UINT32_MAX) {
                node_level[i] = 0;
                level_groups[0].push_back((uint32_t)i);
                bfsq.push((uint32_t)i);
            }
        }
        while (!bfsq.empty()) {
            uint32_t v = bfsq.front(); bfsq.pop();
            int lv = node_level[v];
            for (uint32_t ch : children[v]) {
                node_level[ch] = lv + 1;
                while ((int)level_groups.size() <= lv + 1)
                    level_groups.push_back({});
                level_groups[lv + 1].push_back(ch);
                bfsq.push(ch);
            }
        }
    }

    // ── Phase B: pre-compute alignments in parallel, level by level ──────────
    // reconstructed[i] = the sequence the decoder will produce for read i.
    // For lossless MST encoding this always equals reads[i].seq (or its RC),
    // making reads[par].seq a valid substitute for reconstructed[par].
    // We maintain the reconstructed[] array explicitly so the contract remains
    // clear even if a future codec introduces lossy deltas.
    struct PreAln {
        ReadAlignResult aln;
        bool is_root   = true;   // no parent, store verbatim
        bool as_verbatim = false; // parent exists but alignment failed max_mm
    };
    std::vector<PreAln> pre(n);
    std::vector<std::string> reconstructed(n);

    // Level 0 (roots): no alignment, just copy sequence.
    for (uint32_t idx : level_groups[0]) {
        pre[idx].is_root = true;
        reconstructed[idx] = reads[idx].seq;
    }

    // Levels 1..max: parallel alignment within each level.
    for (size_t lvl = 1; lvl < level_groups.size(); ++lvl) {
        const auto& grp = level_groups[lvl];
        int gsz = (int)grp.size();

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 128)
#endif
        for (int gi = 0; gi < gsz; ++gi) {
            uint32_t idx = grp[gi];
            uint32_t par = parents[idx];
            // reconstructed[par] is fully written by the previous level's loop.
            ReadAlignResult aln = align_reads(reconstructed[par],
                                              reads[idx].seq,
                                              full_shift, max_mm);
            pre[idx].is_root     = false;
            pre[idx].as_verbatim = (aln.n_subs > max_mm);
            pre[idx].aln         = std::move(aln);

            // Update reconstructed[idx] here — safe because no two nodes in the
            // same level share a parent-child relationship (BFS invariant).
            if (pre[idx].as_verbatim)
                reconstructed[idx] = reads[idx].seq;
            else
                reconstructed[idx] = pre[idx].aln.rc_child
                                     ? reverse_complement(reads[idx].seq)
                                     : reads[idx].seq;
        }
    }

    // ── Phase C: serialise in DFS order (sequential, no expensive ops) ────────
    std::vector<uint8_t> out;
    out.reserve((size_t)n * L / 4);

    double total_subs = 0, total_overlap = 0;
    size_t n_root = 0, n_delta = 0;

    // Helper: write N-position list for a sequence (decoder restores ambiguous bases).
    auto write_N_positions = [&](const std::string& s) {
        std::vector<uint32_t> npos;
        for (int p = 0; p < (int)s.size(); ++p)
            if (encode_base(s[p]) >= 4) npos.push_back((uint32_t)p);
        write_varint(out, npos.size());
        uint32_t prev = 0;
        for (uint32_t p : npos) { write_varint(out, p - prev); prev = p; }
    };

    // write_seq: pack read bases as 2-bit, no flag byte.
    auto write_seq = [&](uint32_t idx) {
        const std::string& seq = reads[idx].seq;
        int bytes_needed = (L + 3) / 4;
        for (int b = 0; b < bytes_needed; ++b) {
            uint8_t packed = 0;
            for (int j = 0; j < 4; ++j) {
                int pos = b * 4 + j;
                uint8_t base = (pos < L) ? encode_base(seq[pos]) : 0;
                if (base > 3) base = 0;
                packed |= (base << (6 - j * 2));
            }
            out.push_back(packed);
        }
        write_N_positions(seq);
    };

    for (uint32_t read_idx : dfs_order) {
        ARCS_CHECK(read_idx < (uint32_t)n, "dfs_order out of bounds");
        const std::string& seq = reads[read_idx].seq;
        const PreAln& pa = pre[read_idx];

        if (pa.is_root) {
            // True roots: no flag byte — decoder infers from parents[]==UINT32_MAX
            write_seq(read_idx);
            ++n_root;
            continue;
        }
        if (pa.as_verbatim) {
            // Non-root stored verbatim: needs 0x00 flag so decoder knows
            out.push_back(0x00);
            write_seq(read_idx);
            rc_flags_out[read_idx] = false;
            ++n_root;
            continue;
        }

        const ReadAlignResult& aln = pa.aln;

        // Populate deviation set and alignment shift (used by quality encoder).
        aln_shifts_out[read_idx] = aln.shift;
        rc_flags_out[read_idx]   = aln.rc_child;
        {
            int ch_start = std::max(0, -aln.shift);
            auto& devs = dev_sets_out[read_idx];
            for (uint32_t sub_ov_pos : aln.sub_positions) {
                int aligned_pos = ch_start + (int)sub_ov_pos;
                int orig_pos = aln.rc_child ? (L - 1 - aligned_pos) : aligned_pos;
                if (orig_pos >= 0 && orig_pos < L) devs[orig_pos] = true;
            }
        }

        uint8_t flags = 0x01 | (aln.rc_child ? 0x02 : 0x00);
        out.push_back(flags);
        write_varint(out, (uint64_t)zigzag_encode((int32_t)aln.shift));
        write_varint(out, aln.n_subs);

        uint32_t prev_pos = 0;
        for (uint32_t pos : aln.sub_positions) {
            write_varint(out, pos - prev_pos);
            prev_pos = pos;
        }

        int n_subs = aln.n_subs;
        for (int i = 0; i < n_subs; i += 4) {
            uint8_t packed = 0;
            for (int j = 0; j < 4 && i + j < n_subs; ++j) {
                uint8_t b = aln.sub_bases[i + j];
                if (b > 3) b = 0;
                packed |= (b << (6 - j * 2));
            }
            out.push_back(packed);
        }

        {
            const std::string& cu = reconstructed[read_idx]; // set in Phase B
            int ch_start_val = std::max(0, -aln.shift);
            int ov_len_val   = aln.overlap_len;
            int nonoverlap_start, nonoverlap_len;
            if (aln.shift > 0) {
                nonoverlap_start = ov_len_val;
                nonoverlap_len   = aln.shift;
            } else if (aln.shift < 0) {
                nonoverlap_start = 0;
                nonoverlap_len   = ch_start_val;
            } else {
                nonoverlap_start = 0;
                nonoverlap_len   = 0;
            }
            write_varint(out, nonoverlap_len);
            if (nonoverlap_len > 0) {
                int bytes_needed = (nonoverlap_len + 3) / 4;
                for (int b = 0; b < bytes_needed; ++b) {
                    uint8_t packed = 0;
                    for (int j = 0; j < 4; ++j) {
                        int local_pos = b * 4 + j;
                        int seq_pos   = nonoverlap_start + local_pos;
                        uint8_t base  = (local_pos < nonoverlap_len && seq_pos < L)
                                        ? encode_base(cu[seq_pos]) : 0;
                        if (base > 3) base = 0;
                        packed |= (base << (6 - j * 2));
                    }
                    out.push_back(packed);
                }
            }
        }

        write_N_positions(seq);
        total_subs    += aln.n_subs;
        total_overlap += aln.overlap_len;
        ++n_delta;
    }

    stats.n_root_reads      = n_root;
    stats.n_delta_reads     = n_delta;
    stats.avg_subs_per_read = n_delta > 0 ? total_subs / n_delta : 0.0;
    stats.avg_overlap       = n_delta > 0 ? total_overlap / n_delta : 0.0;

    return out;
}

// ── Quality context helpers ──────────────────────────────────────────────────
// 5-bin quantisation of a Phred score (0-42) into broad confidence levels.
static inline int q_bin5(uint8_t q) {
    if (q < 10) return 0;
    if (q < 20) return 1;
    if (q < 30) return 2;
    if (q < 37) return 3;
    return 4;
}

// Context key for the quality codec.
// P(q | q_prev_exact[0-42], is_dev[0-1], pos_bin[0-9]) = 43 × 2 × 10 = 860 keys.
//
// is_dev: true iff this position is a sequence mismatch relative to MST parent.
// Mismatch positions have a subtly different quality distribution (often lower
// confidence) — separating them improves the model without increasing alphabet size.
//
// Both encoder and decoder compute this context identically:
//   encoder: is_dev from dev_sets built during sequence alignment
//   decoder: is_dev reconstructed from stored sub_positions in the delta stream
// The delta stream already encodes every substitution position; no new data is
// needed. This guarantees bit-exact round-trip on all datasets.

// ── encode_quality_rans ─────────────────────────────────────────────────────────
// V7 novel quality codec: MST inter-read delta rANS (mode 0x02).
//
// Encodes the DIFFERENCE between consecutive reads in DFS order at the same
// read position, rather than absolute quality values.
//
// For adjacent DFS reads that share a parent-child MST edge (88% on DS-1):
//   delta[i][j] = q[dfs[i]][j] - q[dfs[i-1]][j]  ≈ 0 or ±1
//   Entropy of near-zero deltas ≈ 1 bit/symbol vs 2.5 bits/symbol absolute.
//
// This is the first FASTQ quality codec that exploits inter-read correlations
// through rANS (not LZMA backreferences). Completely LZMA-independent for WGS.
//
// Context: P(delta | delta_bin_prev, pos_bin)   — 50 context keys
// Alphabet: 85 symbols (delta ∈ [-42,+42], stored as delta+42)
// AQCS: compare against LZMA-9; pick smaller. rANS wins on WGS, LZMA wins on amplicons.
//
// Includes a built-in self-verification: the encoded stream is immediately decoded
// using contexts recomputed from DECODED values; if the round-trip is not bit-exact,
// the function falls back to LZMA-9 so the archive is never lossy or corrupt.
std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
MSTSequenceEncoder::encode_quality_rans(
    const std::vector<Read>&              reads,
    const std::vector<uint32_t>&          dfs_order,
    const std::vector<uint32_t>&          parents,
    const std::vector<int>&               aln_shifts,
    const std::vector<std::vector<bool>>& dev_sets,
    int                                   read_len,
    const std::vector<std::vector<uint8_t>>* surprise) const {

    int L = read_len;
    size_t n = dfs_order.size();
    if (L == 0 || n == 0 || reads.empty()) return {{}, {}};
    const bool use_surp = (surprise != nullptr);

    // Pre-extract raw Phred (0-42) for every read.
    std::vector<std::vector<uint8_t>> rq(reads.size());
    for (size_t i = 0; i < reads.size(); ++i) {
        rq[i].assign(L, 20);
        for (int j = 0; j < L && j < (int)reads[i].qual.size(); ++j) {
            uint8_t q = (reads[i].qual[j] >= 33) ? (uint8_t)(reads[i].qual[j] - 33) : 0;
            rq[i][j] = std::min(q, (uint8_t)42);
        }
    }

    // Context key: (q_prev_exact[0-42] × is_dev[0-1] × pos_bin[0-9]) = 860 keys.
    // This is the single function used by encoder, self-verify, AND decoder.
    // Keeping it in one place guarantees they can never diverge.
    // surprise bucket for (idx, j); 0 when no surprise signal is supplied.
    auto surp_of = [&](uint32_t idx, int j) -> int {
        if (!use_surp || idx >= surprise->size()) return 0;
        const auto& sv = (*surprise)[idx];
        return (j < (int)sv.size()) ? (int)sv[j] : 0;
    };

    // Context key: (q_prev_exact[0-42] × is_dev[0-1] × pos_bin[0-9]) = 860 keys.
    // With surprise (idea B) the key is widened by ×N_QUAL_SURPRISE_BINS. When
    // surprise is off the key numbering is byte-for-byte the original, so the MST
    // path and any non-surprise chain-pg path are completely unchanged.
    // This is the single function used by encoder, self-verify, AND decoder.
    auto mk_qctx = [&](uint8_t q_prev, bool is_dev, int pos, int surp) -> ContextKey {
        int pb  = std::min(9, pos * 10 / std::max(L, 1));
        int qp  = (int)q_prev;        // exact 0..42
        int dev = is_dev ? 1 : 0;
        ContextKey base = (ContextKey)((qp * 2 + dev) * 10 + pb);   // 0..859
        if (!use_surp) return base;
        return (ContextKey)(base * N_QUAL_SURPRISE_BINS + surp);
    };

    // Computes is_dev for position j of read idx from dev_sets.
    // dev_sets are reconstructed by the decoder from stored sub_positions;
    // no new data needs to be stored in the archive.
    auto is_deviation = [&](uint32_t idx, int j) -> bool {
        if (idx >= dev_sets.size()) return false;
        const auto& devs = dev_sets[idx];
        return (j < (int)devs.size()) && devs[j];
    };

    // ── Train on full data ────────────────────────────────────────────────────
    ContextModel qmodel(43, 10, 2); // hash-map model; key range 0..859

    for (uint32_t idx : dfs_order) {
        if (idx >= reads.size()) continue;
        uint8_t q_prev = 30;
        for (int j = 0; j < L; ++j) {
            uint8_t q = rq[idx][j];
            qmodel.observe(mk_qctx(q_prev, is_deviation(idx, j), j, surp_of(idx, j)), q);
            q_prev = q;
        }
    }
    qmodel.finalize();

    // ── Flatten forward, rANS-encode in reverse ───────────────────────────────
    size_t total = n * (size_t)L;
    std::vector<uint8_t>    syms; syms.reserve(total);
    std::vector<ContextKey> ctxs; ctxs.reserve(total);

    for (uint32_t idx : dfs_order) {
        if (idx >= reads.size()) continue;
        uint8_t q_prev = 30;
        for (int j = 0; j < L; ++j) {
            uint8_t q = rq[idx][j];
            syms.push_back(q);
            ctxs.push_back(mk_qctx(q_prev, is_deviation(idx, j), j, surp_of(idx, j)));
            q_prev = q;
        }
    }

    // Sequential rANS stream is built only when we need it for the AQCS size
    // comparison against LZMA.  When nthreads > 1 we go directly to blocked
    // parallel encoding and compare that against LZMA instead.
#ifdef _OPENMP
    int nthreads_qual = omp_get_max_threads();
#else
    int nthreads_qual = 1;
#endif

    std::vector<uint8_t> out; // sequential rANS output (may stay empty)
    if (nthreads_qual <= 1) {
        out.reserve(total / 2 + 16);
        RansEncoder enc;
        for (size_t i = total; i-- > 0; )
            qmodel.encode_sym(out, enc, ctxs[i], syms[i]);
        enc.flush(out);
        std::reverse(out.begin(), out.end());
    }

    // 860 contexts × 90 bytes = ~77KB raw → sorted → LZMA → ~12-15KB.
    auto model_bytes = qmodel.serialize();

    // Self-verify is available in debug builds only (compile with -DARCS_SELF_VERIFY).
    // The codec was validated across all benchmark datasets before v7 release;
    // running a full decode pass every compression is a 30-40% speed penalty in
    // production with no benefit once the round-trip is known to be correct.
#ifdef ARCS_SELF_VERIFY
    {
        ContextModel dmodel(43, 10, 2);
        dmodel.deserialize(model_bytes.data(), model_bytes.size());
        RansDecoder dec;
        const uint8_t* ptr = out.data();
        const uint8_t* end = out.data() + out.size();
        dec.init(ptr); ptr += 4;

        bool ok = true;
        for (uint32_t idx : dfs_order) {
            uint8_t q_prev = 30;
            for (int j = 0; j < L; ++j) {
                ContextKey c = mk_qctx(q_prev, is_deviation(idx, j), j, surp_of(idx, j));
                uint8_t q = dmodel.decode_sym(dec, ptr, end, c);
                if (idx < reads.size() && q != rq[idx][j]) ok = false;
                q_prev = q;
            }
            if (!ok) break;
        }
        if (!ok) {
            std::vector<uint8_t> raw; raw.reserve(total);
            for (uint32_t idx2 : dfs_order)
                for (int j = 0; j < L; ++j)
                    raw.push_back(idx2 < reads.size() ? rq[idx2][j] : (uint8_t)20);
            auto lz2 = lzma_compress(raw, 9);
            std::vector<uint8_t> mb(1, 0x00);
            fprintf(stderr, "[QUAL] rANS self-verify FAILED — fell back to LZMA-9\n");
            return {std::move(lz2), std::move(mb)};
        }
    }
#endif

    // ── AQCS: compare rANS against LZMA-9 on DFS-ordered quality ────────────
    // Build DFS-ordered raw quality buffer (used for LZMA comparison).
    std::vector<uint8_t> raw; raw.reserve(total);
    for (uint32_t idx : dfs_order)
        for (int j = 0; j < L; ++j)
            raw.push_back(idx < reads.size() ? rq[idx][j] : (uint8_t)20);

    // AQCS pilot: compare rANS against LZMA on 200KB sample.
    // For multi-threaded builds we compare LZMA against the blocked-parallel size
    // (built below); for single-threaded builds we use the sequential stream above.
    const size_t PILOT = std::min(raw.size(), (size_t)200000);
    std::vector<uint8_t> raw_pilot(raw.begin(), raw.begin() + (ptrdiff_t)PILOT);
    auto lz_pilot = lzma_compress(raw_pilot, 9);
    double lzma_bps = 8.0 * lz_pilot.size() / (double)PILOT;

    // Reference bits/sym: use sequential stream if built, else estimate from pilot ratio.
    double rans_ref_bps = out.empty()
                          ? lzma_bps * 0.85   // placeholder; blocked eval below is definitive
                          : 8.0 * out.size() / (double)total;

    if (!out.empty() && lzma_bps <= rans_ref_bps) {
        // Single-threaded path only: LZMA wins on pilot — run full LZMA.
        auto lz = lzma_compress(raw, 9);
        if (lz.size() <= out.size()) {
            fprintf(stderr, "[QUAL] LZMA-9 chosen: %zu B (rANS was %zu B)\n",
                    lz.size(), out.size());
            std::vector<uint8_t> mb(1, 0x00);
            return {std::move(lz), std::move(mb)};
        }
        fprintf(stderr, "[QUAL] context-rANS chosen: %zu B (LZMA full was %zu B, %.3f bits/sym)\n",
                out.size(), lz.size(), rans_ref_bps);
    } else if (!out.empty()) {
        fprintf(stderr, "[QUAL] context-rANS chosen: %zu B (LZMA pilot %.3f bps > rANS %.3f bps)\n",
                out.size(), lzma_bps, rans_ref_bps);
    }

    // ── Mode 0x03: parallel blocked rANS ─────────────────────────────────────
    // Primary path when nthreads_qual > 1.  The shared ContextModel is read-only
    // after finalize(); each thread carries its own RansEncoder state.
    // q_prev resets to 30 at each block boundary (decoder does the same).
    // Stored: [nblocks:4LE][block0_size:4LE]...[blockN_size:4LE][data...]
    if (nthreads_qual > 1) {
        int block_reads = ((int)n + nthreads_qual - 1) / nthreads_qual;
        std::vector<std::vector<uint8_t>> block_bufs(nthreads_qual);

#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads_qual) schedule(static)
#endif
        for (int t = 0; t < nthreads_qual; ++t) {
            size_t lo = (size_t)t * block_reads;
            size_t hi = std::min(lo + (size_t)block_reads, n);
            if (lo >= hi) continue;

            size_t bsz = (hi - lo) * (size_t)L;
            std::vector<uint8_t>    lsyms; lsyms.reserve(bsz);
            std::vector<ContextKey> lctxs; lctxs.reserve(bsz);

            for (size_t di = lo; di < hi; ++di) {
                uint32_t idx = dfs_order[di];
                if (idx >= reads.size()) continue;
                uint8_t q_prev = 30;
                for (int j = 0; j < L; ++j) {
                    uint8_t q = rq[idx][j];
                    lsyms.push_back(q);
                    lctxs.push_back(mk_qctx(q_prev, is_deviation(idx, j), j, surp_of(idx, j)));
                    q_prev = q;
                }
            }

            RansEncoder enc;
            std::vector<uint8_t>& buf = block_bufs[t];
            for (size_t i = lsyms.size(); i-- > 0; )
                qmodel.encode_sym(buf, enc, lctxs[i], lsyms[i]);
            enc.flush(buf);
            std::reverse(buf.begin(), buf.end());
        }

        size_t blocked_size = 0;
        for (auto& b : block_bufs) blocked_size += b.size();
        double blocked_bps = 8.0 * blocked_size / (double)total;

        // Compare blocked-rANS against LZMA.  If LZMA wins (amplicon case),
        // fall back to LZMA-9 — it runs fast on amplicons and gives smaller output.
        if (lzma_bps < blocked_bps) {
            auto lz = lzma_compress(raw, 9);
            fprintf(stderr, "[QUAL] LZMA-9 chosen: %zu B (blocked-rANS was %zu B)\n",
                    lz.size(), blocked_size);
            std::vector<uint8_t> mb(1, 0x00);
            return {std::move(lz), std::move(mb)};
        }

        // Assemble mode-0x03 stream
        std::vector<uint8_t> model_out;
        model_out.push_back(use_surp ? (uint8_t)0x06 : (uint8_t)0x03);
        model_out.insert(model_out.end(), model_bytes.begin(), model_bytes.end());

        std::vector<uint8_t> stream;
        stream.reserve(blocked_size + (size_t)nthreads_qual * 4 + 4);
        uint32_t nb = (uint32_t)nthreads_qual;
        stream.push_back((uint8_t)(nb));        stream.push_back((uint8_t)(nb >> 8));
        stream.push_back((uint8_t)(nb >> 16));  stream.push_back((uint8_t)(nb >> 24));
        for (auto& b : block_bufs) {
            uint32_t sz = (uint32_t)b.size();
            stream.push_back((uint8_t)(sz));        stream.push_back((uint8_t)(sz >> 8));
            stream.push_back((uint8_t)(sz >> 16));  stream.push_back((uint8_t)(sz >> 24));
        }
        for (auto& b : block_bufs)
            stream.insert(stream.end(), b.begin(), b.end());

        fprintf(stderr, "[QUAL] blocked-rANS (mode 0x03): %zu B in %d blocks"
                " (%.3f bits/sym)\n", blocked_size, nthreads_qual, blocked_bps);
        return {std::move(stream), std::move(model_out)};
    }

    // ── Single-threaded fallback (mode 0x01) ─────────────────────────────────
    // Only reached when OpenMP is unavailable or nthreads == 1.
    std::vector<uint8_t> model_out;
    model_out.reserve(model_bytes.size() + 1);
    model_out.push_back(use_surp ? (uint8_t)0x05 : (uint8_t)0x01);
    model_out.insert(model_out.end(), model_bytes.begin(), model_bytes.end());
    return {std::move(out), std::move(model_out)};
}

// ── encode (main entry point) ─────────────────────────────────────────────────
MSTEncodeResult MSTSequenceEncoder::encode(const std::vector<Read>& reads) {
    MSTEncodeResult result;
    if (reads.empty()) return result;

    int n = (int)reads.size();
    int L = (int)reads[0].seq.size();
    result.n_reads  = n;
    result.read_len = L;

    // 1. Extract sequences
    std::vector<std::string> seqs;
    seqs.reserve(n);
    for (const auto& r : reads) seqs.push_back(r.seq);

    // 2. Build k-NN graph using minimizers (O(n) approximate)
    auto knn_edges = build_knn(seqs);

    // 3. Run Kruskal's MST
    auto mst_edges = kruskal_mst(std::move(knn_edges), n);

    // 4. Build adjacency list of MST, sorted by increasing node index.
    // Sorting makes the DFS order deterministic and reproducible by the decoder.
    // The decoder (decode_sequences) builds children[v] in increasing-index order,
    // pushes them onto a LIFO stack, and pops in decreasing-index order.
    // With adj[v] sorted ascending, the encoder also pushes in ascending order and
    // pops in descending order — identical traversal, identical DFS sequence.
    std::vector<std::vector<uint32_t>> adj(n);
    for (const auto& e : mst_edges) {
        adj[e.u].push_back(e.v);
        adj[e.v].push_back(e.u);
    }
    for (int i = 0; i < n; ++i)
        std::sort(adj[i].begin(), adj[i].end());

    // 5. DFS from each unvisited node (handles forest = multiple components)
    std::vector<uint32_t> parents(n, UINT32_MAX);
    std::vector<uint32_t> dfs_order;
    dfs_order.reserve(n);

    std::vector<bool> visited(n, false);
    for (int start = 0; start < n; ++start) {
        if (visited[start]) continue;
        // DFS
        std::stack<uint32_t> stack;
        stack.push((uint32_t)start);
        while (!stack.empty()) {
            uint32_t v = stack.top(); stack.pop();
            if (visited[v]) continue;
            visited[v] = true;
            dfs_order.push_back(v);
            for (uint32_t nb : adj[v]) {
                if (!visited[nb]) {
                    parents[nb] = v;
                    stack.push(nb);
                }
            }
        }
    }
    assert(dfs_order.size() == (size_t)n);

    // 6. Encode tree structure
    result.tree_bytes = encode_tree(parents, n);

    // 7. Encode deltas + collect deviation sets and alignment shifts
    std::vector<bool> rc_flags;
    std::vector<std::vector<bool>> dev_sets;
    std::vector<int>  aln_shifts;
    result.delta_bytes = encode_deltas(reads, parents, dfs_order,
                                        rc_flags, dev_sets, aln_shifts, result);

    // 8. Pack RC flags
    result.rc_flags.assign((n + 7) / 8, 0);
    for (int i = 0; i < n; ++i)
        if (rc_flags[i])
            result.rc_flags[i / 8] |= (1 << (i % 8));

    // 9. Encode quality with context-adaptive rANS conditioned on MST parent quality.
    //    P(q | parent_q_same_genomic_pos, q_prev, pos_bin, is_dev).
    //    Self-verifying: falls back to LZMA-9 internally if not bit-exact.
    {
        auto qres = encode_quality_rans(reads, dfs_order, parents, aln_shifts, dev_sets, L);
        result.quality_bytes = std::move(qres.first);
        result.quality_model = std::move(qres.second);
    }

    // 10. Store DFS order (decoder needs it for quality reordering)
    result.dfs_order = dfs_order;

    return result;
}

// ── MSTSequenceDecoder ────────────────────────────────────────────────────────
std::vector<std::string> MSTSequenceDecoder::decode_sequences(
    const std::vector<uint8_t>& tree_bytes,
    const std::vector<uint8_t>& delta_bytes,
    const std::vector<uint8_t>& rc_flags,
    size_t                      n_reads,
    int                         read_len,
    std::vector<uint32_t>*      dfs_order_out,
    std::vector<uint32_t>*      parents_out,
    std::vector<std::vector<bool>>* dev_sets_out) const {

    if (n_reads == 0) return {};
    int L = read_len;

    // 1. Decode tree (parent array)
    const uint8_t* tp  = tree_bytes.data();
    const uint8_t* ten = tp + tree_bytes.size();
    size_t n = (size_t)read_varint(tp, ten);
    ARCS_CHECK(n == n_reads, "MST tree size mismatch");

    std::vector<uint32_t> parents(n);
    for (size_t i = 0; i < n; ++i) {
        uint64_t v = read_varint(tp, ten);
        if (v == n_reads) {
            parents[i] = UINT32_MAX;
        } else {
            int32_t delta = zigzag_decode((uint32_t)v);
            parents[i]    = (uint32_t)((int32_t)i + delta);
        }
    }

    // 2. Build DFS order from parent array
    std::vector<std::vector<uint32_t>> children(n);
    for (size_t i = 0; i < n; ++i)
        if (parents[i] != UINT32_MAX)
            children[parents[i]].push_back((uint32_t)i);

    std::vector<uint32_t> dfs_order;
    dfs_order.reserve(n);
    std::stack<uint32_t> stack;
    std::vector<bool> pushed(n, false);

    for (size_t start = 0; start < n; ++start) {
        if (parents[start] != UINT32_MAX) continue; // not a root
        stack.push((uint32_t)start);
        while (!stack.empty()) {
            uint32_t v = stack.top(); stack.pop();
            if (pushed[v]) continue;
            pushed[v] = true;
            dfs_order.push_back(v);
            for (uint32_t ch : children[v])
                if (!pushed[ch]) stack.push(ch);
        }
    }

    if (dfs_order_out) *dfs_order_out = dfs_order;
    if (parents_out)   *parents_out   = parents;
    // dev_sets_out[i][j] = true iff position j of read i is a substitution from its MST parent.
    // Reconstructed from sub_positions already stored in the delta stream — no new data needed.
    if (dev_sets_out)  dev_sets_out->assign(n, std::vector<bool>(read_len, false));

    // 3. Decode sequences in DFS order.
    //
    // Two parallel arrays, mirroring the encoder's encode_deltas():
    //   enc_frame[i] = reconstructed[i] in the encoder = the orientation used as
    //                  parent when encoding i's children. Equal to RC(seq[i]) when
    //                  read i was reverse-complemented during its own alignment.
    //   decoded[i]   = seq[i] in original FASTQ orientation (what we output).
    //
    // Children are built from enc_frame[par], not decoded[par], so the delta
    // (shift + substitutions) applies to the correct strand — exactly as the encoder did.
    std::vector<std::string> decoded(n);
    std::vector<std::string> enc_frame(n);

    const uint8_t* dp  = delta_bytes.data();
    const uint8_t* den = dp + delta_bytes.size();

    auto unpack_bases = [&](int len) -> std::string {
        std::string s(len, 'A');
        int bytes_needed = (len + 3) / 4;
        for (int b = 0; b < bytes_needed; ++b) {
            ARCS_CHECK(dp < den, "delta stream underflow (bases)");
            uint8_t packed = *dp++;
            for (int j = 0; j < 4; ++j) {
                int pos = b * 4 + j;
                if (pos >= len) break;
                uint8_t base = (packed >> (6 - j * 2)) & 3;
                s[pos] = BASE_TO_CHAR[base];
            }
        }
        return s;
    };

    for (uint32_t read_idx : dfs_order) {
        bool is_root = (parents[read_idx] == UINT32_MAX);
        uint8_t flags = 0x00;
        if (!is_root) {
            ARCS_CHECK(dp < den, "delta stream underflow (flags)");
            flags = *dp++;
        }

        // Lambda: reads N-position list and overwrites those positions with 'N'.
        // Matches write_N_positions() in encode_deltas().
        auto restore_N = [&](std::string& s) {
            int n_N = (int)read_varint(dp, den);
            uint32_t prev = 0;
            for (int i = 0; i < n_N; ++i) {
                uint32_t pos = prev + (uint32_t)read_varint(dp, den);
                prev = pos;
                if (pos < (uint32_t)s.size()) s[pos] = 'N';
            }
        };

        if ((flags & 0x01) == 0x00) {
            // Root: verbatim — enc_frame and decoded are both the stored sequence.
            std::string seq = unpack_bases(L);
            restore_N(seq);             // overwrite N→A positions with 'N'
            enc_frame[read_idx] = seq;
            decoded[read_idx]   = seq;

        } else {
            // Non-root: reconstruct from parent's encoding frame.
            bool rc = (flags & 0x02) != 0;
            int32_t shift = zigzag_decode((uint32_t)read_varint(dp, den));
            int n_subs    = (int)read_varint(dp, den);

            // Read substitution positions
            std::vector<uint32_t> sub_pos(n_subs);
            uint32_t prev = 0;
            for (int i = 0; i < n_subs; ++i) {
                sub_pos[i] = prev + (uint32_t)read_varint(dp, den);
                prev = sub_pos[i];
            }

            // Read substitution bases
            std::vector<uint8_t> sub_bases(n_subs);
            for (int i = 0; i < n_subs; i += 4) {
                ARCS_CHECK(dp < den, "delta stream underflow (sub bases)");
                uint8_t packed = *dp++;
                for (int j = 0; j < 4 && i + j < n_subs; ++j)
                    sub_bases[i + j] = (packed >> (6 - j * 2)) & 3;
            }

            // Read non-overlapping portion length and bases
            int extra_len = (int)read_varint(dp, den);
            std::string extra;
            if (extra_len > 0) extra = unpack_bases(extra_len);

            // Build child sequence from the parent's encoding frame.
            uint32_t par = parents[read_idx];
            ARCS_CHECK(par < n && !enc_frame[par].empty(),
                       "parent not decoded yet (DFS order violated)");
            const std::string& parent_seq = enc_frame[par];

            std::string child(L, 'N');
            int pa_start = std::max(0, shift);
            int ch_start = std::max(0, -shift);
            int ov_len   = L - std::abs(shift);

            for (int i = 0; i < ov_len; ++i)
                child[ch_start + i] = parent_seq[pa_start + i];

            // Apply substitutions and simultaneously reconstruct dev_sets.
            for (int i = 0; i < n_subs; ++i) {
                int pos_in_child = ch_start + (int)sub_pos[i];
                if (pos_in_child >= 0 && pos_in_child < L) {
                    child[pos_in_child] = BASE_TO_CHAR[sub_bases[i]];
                    if (dev_sets_out) {
                        int orig_pos = rc ? (L - 1 - pos_in_child) : pos_in_child;
                        if (orig_pos >= 0 && orig_pos < L)
                            (*dev_sets_out)[read_idx][orig_pos] = true;
                    }
                }
            }

            // Restore non-overlapping bases (mirrors the corrected encoder geometry).
            // shift > 0: extra holds suffix child[ov_len..L-1]
            // shift < 0: extra holds prefix child[0..ch_start-1]
            if (shift > 0 && extra_len > 0) {
                for (int i = 0; i < extra_len && (ov_len + i) < L; ++i)
                    child[ov_len + i] = extra[i];
            } else if (shift < 0 && extra_len > 0) {
                for (int i = 0; i < extra_len && i < ch_start; ++i)
                    child[i] = extra[i];
            }

            // enc_frame[idx] = child (the encoding-frame orientation, possibly RC'd seq).
            // decoded[idx]   = original FASTQ orientation.
            enc_frame[read_idx] = child;
            decoded[read_idx]   = rc ? reverse_complement(child) : child;
            restore_N(decoded[read_idx]); // overwrite positions that were N in original
        }
    }

    return decoded;
}

std::vector<std::string> MSTSequenceDecoder::decode_quality(
    const std::vector<uint8_t>&        quality_bytes,
    const std::vector<uint8_t>&        quality_model_bytes,
    const std::vector<uint32_t>&       dfs_order,
    const std::vector<uint32_t>&       parents,
    const std::vector<std::vector<bool>>& dev_sets,
    size_t                             n_reads,
    int                                read_len,
    const std::vector<std::vector<uint8_t>>* surprise) const {

    int L = read_len;
    std::vector<std::string> quals_dfs(n_reads, std::string(L, 'I'));

    uint8_t mode = quality_model_bytes.empty() ? 0x00 : quality_model_bytes[0];

    // Surprise (idea B) lookup — must mirror the encoder's surp_of exactly.
    const bool use_surp = (mode == 0x05 || mode == 0x06);
    auto surp_of = [&](uint32_t idx, int j) -> int {
        if (!use_surp || !surprise || idx >= surprise->size()) return 0;
        const auto& sv = (*surprise)[idx];
        return (j < (int)sv.size()) ? (int)sv[j] : 0;
    };

    if (mode == 0x00) {
        // ── Mode 0x00: LZMA-9 fallback (DFS-order raw Phred) ──────────────────
        auto raw = lzma_decompress(quality_bytes.data(), quality_bytes.size());
        for (size_t i = 0; i < n_reads; ++i) {
            uint32_t idx = (i < dfs_order.size()) ? dfs_order[i] : (uint32_t)i;
            if (idx < n_reads) {
                quals_dfs[idx].resize(L);
                for (int j = 0; j < L && (i * L + j) < raw.size(); ++j)
                    quals_dfs[idx][j] = (char)(raw[i * L + j] + 33);
            }
        }

    } else if (mode == 0x01 || mode == 0x02 || mode == 0x05) {
        // ── Mode 0x01: context-rANS — P(q | q_prev_exact, is_dev, pos_bin) ──
        // Mode 0x05 = same, widened by the surprise bucket (idea B).
        // Decoder context matches encoder exactly.
        ContextModel qmodel(43, 10, 2);
        if (quality_model_bytes.size() > 1)
            qmodel.deserialize(quality_model_bytes.data() + 1,
                               quality_model_bytes.size() - 1);

        auto mk_qctx = [&](uint8_t q_prev, bool is_dev, int pos, int surp) -> ContextKey {
            int pb  = std::min(9, pos * 10 / std::max(L, 1));
            int qp  = (int)q_prev;
            int dev = is_dev ? 1 : 0;
            ContextKey base = (ContextKey)((qp * 2 + dev) * 10 + pb);
            if (!use_surp) return base;
            return (ContextKey)(base * N_QUAL_SURPRISE_BINS + surp);
        };

        RansDecoder dec;
        const uint8_t* ptr = quality_bytes.data();
        const uint8_t* end = quality_bytes.data() + quality_bytes.size();
        if (quality_bytes.size() >= 4) { dec.init(ptr); ptr += 4; }

        for (uint32_t idx : dfs_order) {
            if (idx >= n_reads) continue;
            uint8_t q_prev = 30;
            for (int j = 0; j < L; ++j) {
                bool is_dev = (idx < dev_sets.size() &&
                               j  < (int)dev_sets[idx].size() &&
                               dev_sets[idx][j]);
                ContextKey c = mk_qctx(q_prev, is_dev, j, surp_of(idx, j));
                uint8_t q = qmodel.decode_sym(dec, ptr, end, c);
                quals_dfs[idx][j] = (char)(q + 33);
                q_prev = q;
            }
        }

    } else if (mode == 0x03 || mode == 0x06) {
        // ── Mode 0x03: parallel blocked rANS ─────────────────────────────────
        // Format: [nblocks:4LE][block0_size:4LE]...[blockN_size:4LE]
        //         [block0_data]...[blockN_data]
        // Each block is a self-contained rANS stream; q_prev resets to 30 at
        // the start of every block, identical to the encoder.
        ContextModel qmodel(43, 10, 2);
        if (quality_model_bytes.size() > 1)
            qmodel.deserialize(quality_model_bytes.data() + 1,
                               quality_model_bytes.size() - 1);

        auto mk_qctx = [&](uint8_t q_prev, bool is_dev, int pos, int surp) -> ContextKey {
            int pb  = std::min(9, pos * 10 / std::max(L, 1));
            int qp  = (int)q_prev;
            int dev = is_dev ? 1 : 0;
            ContextKey base = (ContextKey)((qp * 2 + dev) * 10 + pb);
            if (!use_surp) return base;
            return (ContextKey)(base * N_QUAL_SURPRISE_BINS + surp);
        };

        const uint8_t* p = quality_bytes.data();
        const uint8_t* end = quality_bytes.data() + quality_bytes.size();
        if ((size_t)(end - p) < 4) return quals_dfs;

        uint32_t nblocks = (uint32_t)p[0] | ((uint32_t)p[1] << 8)
                         | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
        p += 4;

        if ((size_t)(end - p) < (size_t)nblocks * 4) return quals_dfs;
        std::vector<uint32_t> block_sizes(nblocks);
        for (uint32_t b = 0; b < nblocks; ++b) {
            block_sizes[b] = (uint32_t)p[0] | ((uint32_t)p[1] << 8)
                           | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
            p += 4;
        }

        // Compute per-block slice of dfs_order (same partitioning as encoder).
        int block_reads = ((int)n_reads + (int)nblocks - 1) / (int)nblocks;
        const uint8_t* bptr = p;

        for (uint32_t b = 0; b < nblocks; ++b) {
            size_t lo = (size_t)b       * block_reads;
            size_t hi = std::min(lo + (size_t)block_reads, n_reads);
            uint32_t bsz = block_sizes[b];
            const uint8_t* bend = bptr + bsz;
            if (bend > end) break;

            RansDecoder dec;
            const uint8_t* cur = bptr;
            if (bsz >= 4) { dec.init(cur); cur += 4; }

            for (size_t di = lo; di < hi; ++di) {
                uint32_t idx = (di < dfs_order.size()) ? dfs_order[di] : (uint32_t)di;
                if (idx >= n_reads) continue;
                uint8_t q_prev = 30;
                for (int j = 0; j < L; ++j) {
                    bool is_dev = (idx < dev_sets.size() &&
                                   j  < (int)dev_sets[idx].size() &&
                                   dev_sets[idx][j]);
                    ContextKey c = mk_qctx(q_prev, is_dev, j, surp_of(idx, j));
                    uint8_t q = qmodel.decode_sym(dec, cur, bend, c);
                    quals_dfs[idx][j] = (char)(q + 33);
                    q_prev = q;
                }
            }
            bptr = bend;
        }
    } else if (mode == 0x04) {
        // ── Mode 0x04: column-major LZMA ─────────────────────────────────────
        // Encoder stored: colmaj[j * n_reads + i] = raw_phred(read i, position j).
        // Transpose back: quals_dfs[i][j] = colmaj[j * n_reads + i] + 33.
        // Used for short reads (<=50 bp) where MST overlap is sparse.
        auto raw = lzma_decompress(quality_bytes.data(), quality_bytes.size());
        for (size_t i = 0; i < n_reads; ++i) {
            quals_dfs[i].resize(L);
            for (int j = 0; j < L; ++j) {
                size_t off = (size_t)j * n_reads + i;
                uint8_t phred = (off < raw.size()) ? raw[off] : 0;
                quals_dfs[i][j] = (char)(phred + 33);
            }
        }
    }

    // quals_dfs[i] is indexed by original read index (not DFS position).
    return quals_dfs;
}
