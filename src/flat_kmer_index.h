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
    std::vector<uint64_t> km;   // UINT64_MAX = empty slot (fixed-bucket layout only)
    std::vector<uint64_t> val;

    // ── CSR layout (ARCS_FAST_UPGRADE) ───────────────────────────────────────
    // Non-empty `offset` selects it: bucket b owns entries [offset[b],
    // offset[b+1]) of km/val, so every bucket is exactly the size it needs.
    //
    // The fixed-bucket layout above allocates table_size * bucket_cap slots of
    // 16 bytes whether or not they hold anything, and table_size is rounded up
    // to a power of two -- for ~17.5M pseudogenome k-mers that is 1,048,576
    // buckets x 64 slots = 67M slots = 1.07 GB to store 17.5M items, about 61
    // bytes apiece. CSR stores the same items in ~297 MB.
    //
    // It also fixes a correctness wart rather than only a size one: once a
    // fixed bucket fills, further k-mers hashing there are DROPPED, so the
    // index silently forgets occurrences and the fallback placement misses
    // matches it should have found. CSR has no overflow, so nothing is lost.
    std::vector<uint32_t> offset;

    template <class F>
    void for_each(uint64_t query_km, F&& f) const {
        const uint64_t b = flat_kmer_index_mix64(query_km) & table_mask;
        if (!offset.empty()) {
            const uint32_t e = offset[(size_t)b + 1];
            for (uint32_t i = offset[(size_t)b]; i < e; ++i)
                if (km[i] == query_km) f(val[i]);
            return;
        }
        const uint64_t* kp = &km[(size_t)b * (size_t)bucket_cap];
        const uint64_t* vp = &val[(size_t)b * (size_t)bucket_cap];
        for (int s = 0; s < bucket_cap; ++s) {
            if (kp[s] == UINT64_MAX) break;   // compact fill guaranteed by atomic slot claim
            if (kp[s] == query_km) f(vp[s]);
        }
    }

    bool empty() const { return km.empty(); }
};

// ── CSR build: exact-size buckets via counting sort ──────────────────────────
// Same EmitFn contract as build_flat_kmer_index_parallel below, same for_each
// afterwards, but the table is sized to what the data actually needs instead of
// table_size * bucket_cap fixed slots.
//
// Three passes over the staged pairs: count per bucket, prefix-sum into offsets,
// then scatter. Textbook counting sort. Buckets are sized exactly, so no item is
// ever dropped for want of a slot -- unlike the fixed-bucket build, where a full
// bucket silently discards the k-mers that hash to it.
//
// Emission stays parallel (it is the expensive part); the scatter runs serially
// in global order -- thread 0's pairs, then thread 1's, and so on -- so the
// finished index is identical for any thread count.
template <class EmitFn>
FlatKmerIndex build_flat_kmer_index_csr(size_t n_items, int n_threads, EmitFn&& emit) {
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
    for (const auto& v : local) total += v.size();

    FlatKmerIndex idx;
    if (total == 0) { idx.table_mask = 0; return idx; }

    // One bucket per ~4 items: enough to keep chains short, while the table
    // itself stays a small fraction of the entry arrays (4 bytes per bucket
    // against 12 per entry). Entry count is what dominates now, not slack.
    size_t table_size = 1024;
    const size_t want = total / 4 + 1;
    while (table_size < want && table_size < (64ull << 20)) table_size <<= 1;
    idx.table_mask = table_size - 1;

    // total must be addressable by the uint32 offsets; fall back rather than wrap.
    if (total > (size_t)UINT32_MAX) return idx;

    idx.offset.assign(table_size + 1, 0);
    for (const auto& v : local)
        for (const auto& kv : v)
            ++idx.offset[(size_t)(flat_kmer_index_mix64(kv.first) & idx.table_mask) + 1];
    for (size_t b = 0; b < table_size; ++b) idx.offset[b + 1] += idx.offset[b];

    idx.km.resize(total);
    idx.val.resize(total);
    std::vector<uint32_t> cursor(idx.offset.begin(), idx.offset.end() - 1);
    for (const auto& v : local) {
        for (const auto& kv : v) {
            const size_t b = (size_t)(flat_kmer_index_mix64(kv.first) & idx.table_mask);
            const uint32_t at = cursor[b]++;
            idx.km[at]  = kv.first;
            idx.val[at] = kv.second;
        }
    }
    return idx;
}

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

    // Scatter SERIALLY, in global emission order (thread 0's slice, then
    // thread 1's, ...). Each thread's slice is a fixed range of [0, n_items),
    // so that concatenation IS the single global order, which makes the
    // finished index identical no matter how many threads emitted it.
    //
    // This was a parallel scatter using bucket_count[b].fetch_add(). That made
    // the index depend on which thread reached a bucket first -- and not only
    // in slot ORDER: once a bucket reaches bucket_cap the remaining items are
    // DROPPED, so the race decided *which occurrences the index forgets*, and
    // therefore which fallback placement each read found. That is a semantic
    // difference, and it was one of the two reasons the same input did not
    // produce the same archive twice.
    //
    // Serialising costs almost nothing: the expensive phase is emit() computing
    // the k-mers (still parallel above); this loop only hashes and stores, over
    // an item count on the order of the pseudogenome length.
    std::vector<uint32_t> bucket_count(table_size, 0);
    for (int t = 0; t < T; ++t) {
        for (auto& kv : local[(size_t)t]) {
            uint64_t b = flat_kmer_index_mix64(kv.first) & idx.table_mask;
            uint32_t slot = bucket_count[b]++;
            if (slot < (uint32_t)bucket_cap) {
                size_t pos = (size_t)b * (size_t)bucket_cap + slot;
                idx.km[pos]  = kv.first;
                idx.val[pos] = kv.second;
            }
        }
    }
    return idx;
}
