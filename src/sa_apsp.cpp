#include "sa_apsp.h"
#include "libsais/libsais.h"
#include "libsais/libsais64.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <thread>

// ── Suffix array construction: prefix doubling (Karp-Miller-Rosenberg) with
// radix/counting sort per round — true O(n log n), not comparison-sort based
// (which would add another log n factor and isn't tractable at the scale
// suffix arrays here actually need to run at, hundreds of millions of chars).
void SuffixArray::build(const std::string& text) {
    const int n = (int)text.size();
    sa.assign((size_t)n, 0);
    lcp.assign((size_t)n, 0);
    if (n == 0) return;
    if (n == 1) { sa[0] = 0; return; }

    std::vector<int> rank((size_t)n), tmp((size_t)n);
    for (int i = 0; i < n; ++i) { sa[(size_t)i] = (uint32_t)i; rank[(size_t)i] = (unsigned char)text[(size_t)i]; }

    // Initial counting sort by single byte (alphabet size 256).
    {
        std::vector<int> cnt(257, 0);
        for (int i = 0; i < n; ++i) cnt[(size_t)rank[(size_t)i] + 1]++;
        for (int v = 1; v <= 256; ++v) cnt[(size_t)v] += cnt[(size_t)v - 1];
        std::vector<uint32_t> sorted(n);
        for (int i = 0; i < n; ++i) sorted[(size_t)cnt[(size_t)rank[(size_t)i]]++] = (uint32_t)i;
        sa = sorted;
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; ++i)
            tmp[sa[(size_t)i]] = tmp[sa[(size_t)i - 1]] + (rank[sa[(size_t)i]] != rank[sa[(size_t)i - 1]] ? 1 : 0);
        rank = tmp;
    }

    for (int k = 1; rank[sa[(size_t)n - 1]] < n - 1; k <<= 1) {
        // key2key(i) in [0, n]: 0 means "no second half" (i+k>=n), else rank+1
        // (shifted so it's non-negative — avoids ever casting -1 to size_t).
        auto key2key = [&](int i) -> int { return (i + k < n) ? (rank[(size_t)(i + k)] + 1) : 0; };

        // Standard stable counting sort, backward iteration + pre-decrement —
        // avoids the off-by-one ambiguity of a forward/post-increment variant.
        // Pass 1 (LSD): sort ALL indices by key2key.
        std::vector<uint32_t> by_key2((size_t)n);
        {
            std::vector<int> cnt((size_t)n + 1, 0);
            for (int i = 0; i < n; ++i) cnt[(size_t)key2key(i)]++;
            for (int v = 1; v < (int)cnt.size(); ++v) cnt[(size_t)v] += cnt[(size_t)v - 1];
            for (int i = n - 1; i >= 0; --i) by_key2[(size_t)(--cnt[(size_t)key2key(i)])] = (uint32_t)i;
        }
        // Pass 2 (MSD): stable-sort that order by rank (the more significant
        // key) — composes into a correct full sort by (rank, key2key).
        {
            std::vector<int> cnt((size_t)n + 1, 0);
            for (int i = 0; i < n; ++i) cnt[(size_t)rank[(size_t)i]]++;
            for (int v = 1; v < (int)cnt.size(); ++v) cnt[(size_t)v] += cnt[(size_t)v - 1];
            std::vector<uint32_t> new_sa((size_t)n);
            for (int idx = n - 1; idx >= 0; --idx) {
                uint32_t i = by_key2[(size_t)idx];
                new_sa[(size_t)(--cnt[(size_t)rank[(size_t)i]])] = i;
            }
            sa = new_sa;
        }

        tmp[sa[0]] = 0;
        for (int i = 1; i < n; ++i) {
            bool same = rank[sa[(size_t)i]] == rank[sa[(size_t)i - 1]] &&
                        key2key((int)sa[(size_t)i]) == key2key((int)sa[(size_t)i - 1]);
            tmp[sa[(size_t)i]] = tmp[sa[(size_t)i - 1]] + (same ? 0 : 1);
        }
        rank = tmp;
    }

    // Kasai's algorithm: LCP array in O(n), using the fully-refined rank[]
    // (now the inverse permutation of sa[]) computed above.
    int h = 0;
    for (int i = 0; i < n; ++i) {
        if (rank[(size_t)i] > 0) {
            int j = (int)sa[(size_t)rank[(size_t)i] - 1];
            while (i + h < n && j + h < n && text[(size_t)(i + h)] == text[(size_t)(j + h)]) ++h;
            lcp[(size_t)rank[(size_t)i]] = (uint32_t)h;
            if (h > 0) --h;
        } else {
            h = 0;
        }
    }
}

void SuffixArray::build_libsais(const std::string& text, int threads) {
    // Size computed via size_t from the start — never cast to a signed
    // 32-bit int before this point, which is exactly what let a >2^31-char
    // real input (GIAB HG002, ~3.8B chars) silently wrap to a negative
    // value and crash downstream in the previous version of this function.
    const size_t n64 = text.size();
    if (n64 > (size_t)UINT32_MAX) {
        // Our sa/lcp storage is deliberately uint32_t (see sa_apsp.h) — no
        // realistic dataset should ever reach this, but fail loudly rather
        // than silently truncate if one somehow does.
        throw std::runtime_error("SuffixArray::build_libsais: text exceeds uint32_t position range (need wider storage)");
    }
    sa.assign(n64, 0);
    lcp.assign(n64, 0);
    if (n64 == 0) return;
    if (n64 == 1) { sa[0] = 0; return; }

    if (n64 <= (size_t)INT32_MAX) {
        // Fits the 32-bit libsais API directly — and, crucially, int32_t and
        // uint32_t are the same width, so libsais can write its output STRAIGHT
        // into our own sa/lcp storage. This previously allocated a separate
        // int32 SA and LCP and then copied 4n+4n bytes across, which held two
        // full copies of each array live at once for no reason: 8 bytes per
        // char of text, ~2.1 GB on a 257 Mchar read set. Every value libsais
        // writes here is a text position or an LCP length, both non-negative
        // and < n <= INT32_MAX, so the reinterpretation is exact in both
        // directions (it is the same object representation, not a conversion).
        // PLCP is genuinely transient — libsais_lcp_omp consumes it to produce
        // LCP — so it keeps its own buffer, scoped to be freed immediately.
        const int32_t n = (int32_t)n64;
        int32_t* SAp  = reinterpret_cast<int32_t*>(sa.data());
        int32_t* LCPp = reinterpret_cast<int32_t*>(lcp.data());
        int rc1 = libsais_gsa_omp((const uint8_t*)text.data(), SAp, n, 0, nullptr, threads);
        if (rc1 != 0) throw std::runtime_error("libsais_gsa_omp failed");
        {
            std::vector<int32_t> PLCP((size_t)n);
            int rc2 = libsais_plcp_gsa_omp((const uint8_t*)text.data(), SAp, PLCP.data(), n, threads);
            if (rc2 != 0) throw std::runtime_error("libsais_plcp_gsa_omp failed");
            int rc3 = libsais_lcp_omp(PLCP.data(), SAp, LCPp, n, threads);
            if (rc3 != 0) throw std::runtime_error("libsais_lcp_omp failed");
        }
    } else {
        // Text exceeds INT32_MAX chars (~2.147B) — the real case that
        // crashed on GIAB HG002 (~3.8B chars). Use libsais64's int64_t-indexed
        // API for construction, then downcast into our uint32_t storage
        // (safe: n64 already confirmed <= UINT32_MAX above).
        const int64_t n = (int64_t)n64;
        std::vector<int64_t> SA((size_t)n), PLCP((size_t)n), LCP((size_t)n);
        int64_t rc1 = libsais64_gsa_omp((const uint8_t*)text.data(), SA.data(), n, 0, nullptr, threads);
        if (rc1 != 0) throw std::runtime_error("libsais64_gsa_omp failed");
        int64_t rc2 = libsais64_plcp_gsa_omp((const uint8_t*)text.data(), SA.data(), PLCP.data(), n, threads);
        if (rc2 != 0) throw std::runtime_error("libsais64_plcp_gsa_omp failed");
        int64_t rc3 = libsais64_lcp_omp(PLCP.data(), SA.data(), LCP.data(), n, threads);
        if (rc3 != 0) throw std::runtime_error("libsais64_lcp_omp failed");
        for (size_t i = 0; i < n64; ++i) sa[i]  = (uint32_t)SA[i];
        for (size_t i = 0; i < n64; ++i) lcp[i] = (uint32_t)LCP[i];
    }
}

namespace {
// SEP is byte value 0 — guaranteed distinct from every real base byte the pg
// pipeline ever puts through here (encode_base rejects non-ACGT already, and
// reads_both_views entries are always plain ACGT strings), so a suffix
// starting at (or comparison reaching) a separator can never spuriously
// match real sequence content: Kasai's algorithm does a real byte comparison,
// and SEP != any ACGT byte, so LCP naturally stops there — no boundary bug to
// guard against separately.
constexpr char SEP = '\0';
}

std::vector<std::vector<APSPCandidate>> build_apsp_candidates(
    const std::vector<std::string>& reads_both_views,
    uint32_t n_reads, int max_cands, uint32_t min_overlap, uint32_t search_cap) {

    const size_t m = reads_both_views.size(); // = 2*n_reads
    std::vector<uint32_t> seg_start(m), seg_len(m);
    std::string T;
    // pos_to_seg[pos] = which segment (index into seg_start/seg_len) owns
    // absolute text position `pos` — one entry per character of T, built in
    // the SAME linear, sequential-write pass that builds T itself below (no
    // extra asymptotic cost). Exists to replace locate()'s binary search:
    // real profiling (perf record -g -e cache-misses, see session notes) on
    // this exact function found locate()'s binary search over `seg_start`
    // (m ~ millions of entries) to be the single largest cache-miss
    // contributor in the ENTIRE compress pipeline — ~22-23% of every cache
    // miss in the whole program — because SA-order visits positions in
    // content order, not text-position order, so successive locate() calls
    // during a walk have no spatial locality to exploit (verified, not
    // assumed). A direct array read replaces ~log2(m) quasi-random
    // comparisons with one O(1) lookup. Real memory cost: 4 bytes per
    // character of T (e.g. ~1.2GB at this project's ~300M-char yeast test
    // scale) — a genuine, deliberate trade of memory for the confirmed
    // speed win, not free; at far larger inputs (hundreds of GB text) this
    // would add proportionally, which matters for that scale's own
    // already-tighter memory budget (see SuffixArray::build_libsais notes).
    std::vector<uint32_t> pos_to_seg;
    {
        size_t total = 0;
        for (auto& s : reads_both_views) total += s.size() + 1;
        // seg_start stores absolute positions in T as uint32_t (matching
        // SuffixArray's own storage decision, see sa_apsp.h) — guard against
        // the total exceeding that range up front, rather than silently
        // wrapping partway through the build below.
        if (total > (size_t)UINT32_MAX) {
            throw std::runtime_error("build_apsp_candidates: concatenated text exceeds uint32_t position range");
        }
        T.reserve(total);
        pos_to_seg.resize(total);
        for (size_t i = 0; i < m; ++i) {
            seg_start[i] = (uint32_t)T.size();
            seg_len[i]   = (uint32_t)reads_both_views[i].size();
            T += reads_both_views[i];
            T += SEP;
            // Fill this segment's full footprint (content + its own trailing
            // SEP byte) — matches locate()'s existing convention of
            // resolving the separator position to its OWNING segment too
            // (off == seg_len[e] there signals "this is the separator").
            std::fill(pos_to_seg.begin() + seg_start[i],
                      pos_to_seg.begin() + seg_start[i] + seg_len[i] + 1,
                      (uint32_t)i);
        }
    }

    const bool SA_TIMING = getenv("ARCS_VODBG_TIMING") != nullptr;
    auto _sa_t0 = std::chrono::steady_clock::now();

    // libsais (O(n) induced-sorting, OpenMP-parallel SA+PLCP+LCP) is now the
    // default — confirmed via test_libsais_apsp.cpp that its derived APSP
    // candidates are brute-force-correct despite differing from our own
    // prefix-doubling build()'s raw (sa,lcp) arrays at separator-tie
    // positions (see SuffixArray::build_libsais's comment). Our own
    // build() stays available as ARCS_VODBG_SA_SERIAL=1 for comparison/
    // fallback, not removed.
    SuffixArray SA;
    int sa_threads = (int)std::thread::hardware_concurrency();
    if (sa_threads < 1) sa_threads = 1;
    if (const char* s = getenv("ARCS_VODBG_SA_THREADS")) { int v = atoi(s); if (v >= 1) sa_threads = v; }
    if (getenv("ARCS_VODBG_SA_SERIAL")) {
        // Legacy prefix-doubling path — comparison/fallback only, never the
        // default. It internally uses signed 32-bit indices throughout and
        // was never widened (not worth it for a superseded fallback), so
        // guard it explicitly here rather than let it silently corrupt on
        // the exact class of large input that broke the libsais path before
        // it was fixed.
        if (T.size() > (size_t)INT32_MAX) {
            throw std::runtime_error("build_apsp_candidates: ARCS_VODBG_SA_SERIAL (legacy prefix-doubling) "
                                      "does not support texts beyond INT32_MAX chars; unset it to use libsais");
        }
        SA.build(T);
    } else {
        SA.build_libsais(T, sa_threads);
    }
    // Root fix for the crash found on real GIAB HG002 data (~3.8B chars):
    // this was previously `const int n = (int)T.size();` — a premature
    // signed-32-bit cast that silently wrapped negative before any size
    // check ever ran, producing a bogus huge size_t a few lines down and
    // crashing with "cannot create std::vector larger than max_size()".
    // T.size() itself is already guarded to fit uint32_t above; int64_t
    // here is just the safe common type for the loop arithmetic below.
    const int64_t n = (int64_t)T.size();
    if (SA_TIMING) {
        auto t1 = std::chrono::steady_clock::now();
        fprintf(stderr, "[VB-TIMING]   sa_construct(%s): %.2fs (n=%lld, threads=%d)\n",
                getenv("ARCS_VODBG_SA_SERIAL") ? "prefix-doubling+Kasai" : "libsais",
                std::chrono::duration<double>(t1 - _sa_t0).count(), (long long)n, sa_threads);
        _sa_t0 = t1;
    }

    // rank_at_seg[e] = position in SA order of the suffix that starts at
    // seg_start[e] (i.e. the whole of entry e's read).
    //
    // This used to be a full inverse suffix array, rank_of[pos] for every one
    // of the n text positions -- 4 bytes per char, ~1.0 GB on a 257 Mchar read
    // set. But it has exactly one reader (the `uint32_t r = rank_at_seg[e];`
    // below), and that reader only ever asks about a segment START. So 99.3%
    // of the entries were built and never read: m is ~1.7M against n ~257M.
    // Storing one entry per segment answers every query that is actually made,
    // for 1/150th of the memory and better locality.
    //
    // Built in the same single pass over SA: for each suffix position p, ask
    // pos_to_seg which entry owns p and keep it only if p is that entry's
    // start. One extra comparison per position, no extra pass.
    std::vector<uint32_t> rank_at_seg(m, 0);
    for (int64_t i = 0; i < n; ++i) {
        const uint32_t p = SA.sa[(size_t)i];
        const uint32_t e = pos_to_seg[p];
        if (seg_start[e] == p) rank_at_seg[e] = (uint32_t)i;
    }

    // pos -> (entry index, offset within that entry), or npos if pos lands on
    // a separator byte (not a valid read-content position). Was a binary
    // search over seg_start (confirmed via real perf profiling to be the
    // single largest cache-miss contributor in the whole compress pipeline —
    // see pos_to_seg's own comment above); now a direct O(1) array read.
    auto locate = [&](uint32_t pos) -> int64_t {
        size_t e = pos_to_seg[pos];
        uint32_t off = pos - seg_start[e];
        if (off >= seg_len[e]) return -1; // the separator byte itself
        return (int64_t)((e << 32) | off);
    };

    std::vector<std::vector<APSPCandidate>> out(m);

    // Per-entry APSP discovery is embarrassingly parallel: entry e only ever
    // reads the already-fully-built, shared, read-only SA/LCP/rank_at_seg/locate
    // structures and writes exclusively to out[e] — no other entry's slot is
    // ever touched. Unlike Method A's read-sharding mistake, splitting the
    // RANGE OF e here changes nothing about what any single e can see (SA
    // already covers every read globally) — it's pure parallel-map, zero
    // coordination needed, not a coverage trade-off.
    int apsp_threads = (int)std::thread::hardware_concurrency();
    if (apsp_threads < 1) apsp_threads = 1;
    if (const char* s = getenv("ARCS_VODBG_APSP_THREADS")) { int v = atoi(s); if (v >= 1) apsp_threads = v; }
    if ((size_t)apsp_threads > m) apsp_threads = std::max(1, (int)m);

    auto apsp_worker = [&](size_t lo, size_t hi) {
    for (size_t e = lo; e < hi; ++e) {
        uint32_t Lr = seg_len[e];
        if (Lr < min_overlap) continue;
        uint32_t r = rank_at_seg[e];   // SA rank of this entry's full-read suffix

        std::vector<APSPCandidate> cands;
        // `running` is the raw LCP between the two suffixes, which can be
        // spuriously inflated by exactly one separator byte matching another
        // separator byte (both are the same value) right after two reads of
        // otherwise-equal length — real content never causes this since SEP
        // differs from every ACGT byte, but read_e's separator and read_oe's
        // separator can coincide one position past where the real match
        // should stop. So the overlap length isn't "however far raw LCP
        // reaches" — it's *determined* by where oe ends (off + ov must equal
        // seg_len[oe] by definition of a suffix-prefix overlap); accept iff
        // the real common-prefix evidence (running, capped at e's own length
        // too) covers at least that far, and use the derived length, not the
        // possibly-inflated raw one.
        auto consider = [&](int64_t loc, uint32_t running) {
            if (loc < 0 || running < min_overlap) return;
            size_t oe  = (size_t)((uint64_t)loc >> 32);
            uint32_t off = (uint32_t)((uint64_t)loc & 0xFFFFFFFFu);
            if (oe == e) return; // self
            if (off >= seg_len[oe]) return;
            uint32_t ov = seg_len[oe] - off; // the only overlap length that makes this a real suffix-prefix match
            if (ov < min_overlap || ov > Lr || ov > running) return;
            cands.push_back({(uint32_t)(oe / 2), (uint8_t)(oe % 2), ov});
        };

        // Walk left (toward smaller SA rank), tracking running-min LCP.
        {
            uint32_t running = UINT32_MAX;
            uint32_t steps = 0;
            for (int64_t i = (int64_t)r - 1; i >= 0 && steps < search_cap; --i, ++steps) {
                running = std::min(running, SA.lcp[(size_t)(i + 1)]);
                if (running < min_overlap) break;
                consider(locate(SA.sa[(size_t)i]), running);
            }
        }
        // Walk right (toward larger SA rank).
        {
            uint32_t running = UINT32_MAX;
            uint32_t steps = 0;
            for (int64_t i = (int64_t)r + 1; i < n && steps < search_cap; ++i, ++steps) {
                running = std::min(running, SA.lcp[(size_t)i]);
                if (running < min_overlap) break;
                consider(locate(SA.sa[(size_t)i]), running);
            }
        }

        std::sort(cands.begin(), cands.end(), [](const APSPCandidate& a, const APSPCandidate& b) {
            return a.overlap > b.overlap;
        });
        if ((int)cands.size() > max_cands) cands.resize((size_t)max_cands);
        out[e] = std::move(cands);
    }
    }; // end apsp_worker

    if (apsp_threads <= 1) {
        apsp_worker(0, m);
    } else {
        std::vector<std::thread> ths;
        ths.reserve((size_t)apsp_threads);
        for (int t = 0; t < apsp_threads; ++t) {
            size_t lo = m * (size_t)t / (size_t)apsp_threads;
            size_t hi = m * (size_t)(t + 1) / (size_t)apsp_threads;
            ths.emplace_back(apsp_worker, lo, hi);
        }
        for (auto& th : ths) th.join();
    }
    if (SA_TIMING) {
        auto t1 = std::chrono::steady_clock::now();
        fprintf(stderr, "[VB-TIMING]   apsp_walk: %.2fs (threads=%d, m=%zu)\n",
                std::chrono::duration<double>(t1 - _sa_t0).count(), apsp_threads, m);
        _sa_t0 = t1;
    }

    return out;
}
