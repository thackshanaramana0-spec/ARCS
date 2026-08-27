#pragma once
#include <atomic>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

// ── Generic parallel-built flat k-mer occurrence table ───────────────────────
// Extracted from kmer_anchor_pg.cpp's FlatKmerIndex/build_index_parallel
// (Method A's proven index-build parallelization: ~15-18x speedup, zero
// quality cost, confirmed lossless) — generalized here to carry an arbitrary
// uint64_t value per (k-mer, value) pair, not just a fixed rid/view/off
// packing, so it can also serve Method B's fallback placement pass (which
// indexes PG POSITIONS, not per-read occurrences) without duplicating the
// two-phase counting-sort-scatter logic a second time.
//
// Two-phase build: pass 1, each thread independently emits its own local
// (k-mer, value) pairs for its slice of the input (via the caller-supplied
// EmitFn) — no shared state, no locks. Pass 2, every thread's local pairs are
// scattered into ONE shared flat table using a single atomic fetch_add per
// bucket to claim a slot — zero retries, zero true locks, because each
// insert's destination bucket is independent of every other insert's.
//
// km is built via repeated (v<<2)|base, so its LOW bits are dominated by the
// LAST few bases only — using them directly as a bucket index (raw `& mask`)
// causes very uneven bucket occupancy (many distinct k-mers sharing a common
// suffix motif collide into the same bucket and fight over its bucket_cap
// slots). Fix: a real avalanche mix (splitmix64 finalizer) before masking, so
// bucket assignment no longer correlates with k-mer suffix content — this was
// found and fixed the hard way during Method A's own index-build work
// (a first attempt without this regressed pg_len 34.0M -> 49.0M).
inline uint64_t flat_kmer_index_mix64(uint64_t x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

struct FlatKmerIndex {
    uint64_t table_mask = 0;
    int      bucket_cap = 0;
    std::vector<uint64_t> km;   // UINT64_MAX = empty slot
    std::vector<uint64_t> val;

    template <class F>
    void for_each(uint64_t query_km, F&& f) const {
        uint64_t b = flat_kmer_index_mix64(query_km) & table_mask;
        const uint64_t* kp = &km[(size_t)b * (size_t)bucket_cap];
        const uint64_t* vp = &val[(size_t)b * (size_t)bucket_cap];
        for (int s = 0; s < bucket_cap; ++s) {
            if (kp[s] == UINT64_MAX) break;   // compact fill guaranteed by atomic slot claim
            if (kp[s] == query_km) f(vp[s]);
        }
    }

    bool empty() const { return km.empty(); }
};

// EmitFn signature: void(size_t lo, size_t hi, std::vector<std::pair<uint64_t,uint64_t>>& out)
// — called once per thread with that thread's [lo,hi) slice of [0, n_items);
// must push_back every (kmer, value) pair its slice contributes.
template <class EmitFn>
FlatKmerIndex build_flat_kmer_index_parallel(size_t n_items, int bucket_cap,
                                              int n_threads, EmitFn&& emit) {
    int T = std::max(1, n_threads);
    if ((size_t)T > std::max((size_t)1, n_items)) T = std::max(1, (int)n_items);
    std::vector<std::vector<std::pair<uint64_t,uint64_t>>> local((size_t)T);
    {
        std::vector<std::thread> ths;
        ths.reserve((size_t)T);
        for (int t = 0; t < T; ++t) {
            ths.emplace_back([&, t]{
                size_t lo = n_items * (size_t)t / (size_t)T;
                size_t hi = n_items * (size_t)(t + 1) / (size_t)T;
                emit(lo, hi, local[(size_t)t]);
            });
        }
        for (auto& th : ths) th.join();
    }

    size_t total = 0;
    for (auto& v : local) total += v.size();

    // Real memory cost is table_size * bucket_cap slots (each bucket carries
    // bucket_cap inline entries), not table_size alone — size so that TOTAL
    // SLOT CAPACITY is ~2x the occurrence count (average bucket well under
    // half full), then cap total slot count outright as a hard memory guard.
    size_t target_buckets = (total * 2) / (size_t)std::max(1, bucket_cap) + 1;
    size_t table_size = 1024;
    while (table_size < target_buckets && table_size < (64ull << 20)) table_size <<= 1;
    while (table_size * (size_t)bucket_cap > (2ull << 30)) table_size >>= 1; // hard cap: 2G slots total

    FlatKmerIndex idx;
    idx.table_mask = table_size - 1;
    idx.bucket_cap = bucket_cap;
    idx.km.assign(table_size * (size_t)bucket_cap, UINT64_MAX);
    idx.val.assign(table_size * (size_t)bucket_cap, 0);

    std::vector<std::atomic<uint32_t>> bucket_count(table_size);
    for (auto& c : bucket_count) c.store(0, std::memory_order_relaxed);

    {
        std::vector<std::thread> ths;
        ths.reserve((size_t)T);
        for (int t = 0; t < T; ++t) {
            ths.emplace_back([&, t]{
                for (auto& kv : local[(size_t)t]) {
                    uint64_t b = flat_kmer_index_mix64(kv.first) & idx.table_mask;
                    uint32_t slot = bucket_count[b].fetch_add(1, std::memory_order_relaxed);
                    if (slot < (uint32_t)bucket_cap) {
                        size_t pos = (size_t)b * (size_t)bucket_cap + slot;
                        idx.km[pos]  = kv.first;
                        idx.val[pos] = kv.second;
                    }
                }
            });
        }
        for (auto& th : ths) th.join();
    }
    return idx;
}
