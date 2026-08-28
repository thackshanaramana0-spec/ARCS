#include "arcs_threads.h"
#include "repeat_elim.h"
#include "sa_apsp.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <vector>

namespace {

// pg text here is always pure ACGT (see chain_encoder.h: "pg (ACGT only, N ->
// 'A')"), so a simple 4-way lookup is exact — no need for the fuller
// reverse_complement() used elsewhere for arbitrary IUPAC/lowercase input.
inline char complement_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return c;
    }
}

struct BestMatch {
    uint32_t len   = 0;
    uint32_t src   = 0;
    uint8_t  is_rc = 0;
};

}  // namespace

bool repeat_elim_encode(const std::string& pg, uint32_t K,
                         std::string& literal_out,
                         std::vector<RepeatMatch>& matches_out) {
    return repeat_elim_encode_cap(pg, K, 64ULL * 1024 * 1024, literal_out, matches_out);
}

bool repeat_elim_encode_cap(const std::string& pg, uint32_t K,
                            uint64_t /*cap_entries*/,  // vestigial: the old sparse-hash
                            // index used this to bound table size; the SA/LPF approach
                            // below is exact at every position, no sampling density to
                            // cap. Kept in the signature so callers (encoder.cpp) and its
                            // ARCS_REPEAT_ELIM_CAP knob don't need touching.
                            std::string& literal_out,
                            std::vector<RepeatMatch>& matches_out) {
    literal_out.clear();
    matches_out.clear();
    const size_t n = pg.size();

    if (n > 0xFFFFFFFFULL) { literal_out = pg; return true; }  // wire offsets are u32
    if (n < K || K < 8)     { literal_out = pg; return true; }

    // T's size (2*n + 2) must itself fit our uint32_t-indexed SA/LCP storage
    // (see sa_apsp.h / SuffixArray::build_libsais's own UINT32_MAX guard) —
    // bail to a plain literal pass (never incorrect, just misses the
    // opportunity) rather than risk it. repeat_elim is an opt-in research
    // pass, not exercised at multi-billion-char pg scale today.
    if ((uint64_t)n * 2 + 2 > (uint64_t)UINT32_MAX) { literal_out = pg; return true; }

    // T = pg + SEP + reverse_complement(pg) + SEP — see repeat_elim.h for
    // why this single concatenation lets one suffix array answer both
    // "forward repeat" and "reverse-complement repeat" queries.
    std::string T;
    T.reserve(2 * n + 2);
    T.append(pg);
    T.push_back('\0');
    for (size_t t = 0; t < n; ++t) T.push_back(complement_base(pg[n - 1 - t]));
    T.push_back('\0');

    int threads = arcs_threads();
    if (threads < 1) threads = 1;
    if (const char* s = getenv("ARCS_REPEAT_ELIM_SA_THREADS")) { int v = atoi(s); if (v >= 1) threads = v; }

    SuffixArray SA;
    SA.build_libsais(T, threads);

    const size_t Tn = T.size();
    std::vector<uint32_t> rank_of(Tn);
    for (size_t i = 0; i < Tn; ++i) rank_of[SA.sa[i]] = (uint32_t)i;

    uint32_t search_cap = 2000;
    if (const char* s = getenv("ARCS_REPEAT_ELIM_SEARCHCAP")) { long v = atol(s); if (v >= 10) search_cap = (uint32_t)v; }

    // For every pg position i, the best available earlier match (forward or
    // RC) depends only on the static SA/LCP/rank_of tables built above —
    // never on which OTHER positions the greedy parse below ends up
    // consuming — so it's embarrassingly parallel across i, same principle
    // as build_apsp_candidates's own per-entry walk.
    std::vector<BestMatch> bm(n);
    auto worker = [&](size_t lo, size_t hi) {
        for (size_t i = lo; i < hi; ++i) {
            if (i + K > n) continue;
            uint32_t r = rank_of[i];
            uint32_t best_len = 0, best_src = 0;
            uint8_t  best_rc  = 0;

            auto consider_forward = [&](uint32_t p, uint32_t running) {
                if (p >= (uint32_t)i) return;  // only strictly-earlier, already-decoded content is usable
                uint32_t usable = std::min(running, (uint32_t)(i - p));
                if (usable > best_len) { best_len = usable; best_src = p; best_rc = 0; }
            };
            // A match to reverse_complement(pg[s_end-len : s_end)) shows up as a
            // plain forward LCP between T-suffix(i) and T-suffix at Q-offset
            // qoff, where s_end = n - qoff is FIXED per candidate (independent of
            // the eventual match length: RC(pg[s:s+len))[t] = complement(pg[s+len-1-t]),
            // and complement(pg[u]) = Q[n-1-u], so pg[i+t] = Q[n-s_end+t] with
            // s_end = s+len held constant as t, and hence len, varies) — so only
            // s_end <= i (the source range ends at/before the current position)
            // needs checking, exactly like the forward case.
            auto consider_rc = [&](uint32_t qoff, uint32_t running) {
                uint32_t s_end = (uint32_t)n - qoff;
                if (s_end > (uint32_t)i) return;
                uint32_t usable = std::min(running, s_end);
                if (usable > best_len) { best_len = usable; best_src = s_end - usable; best_rc = 1; }
            };

            uint32_t running = UINT32_MAX;
            uint32_t steps = 0;
            for (int64_t rr = (int64_t)r - 1; rr >= 0 && steps < search_cap; --rr, ++steps) {
                running = std::min(running, (uint32_t)SA.lcp[(size_t)(rr + 1)]);
                if (running < K) break;
                uint32_t p = SA.sa[(size_t)rr];
                if (p < (uint32_t)n) consider_forward(p, running);
                else if (p > (uint32_t)n && p <= (uint32_t)(2 * n)) consider_rc(p - (uint32_t)(n + 1), running);
            }
            running = UINT32_MAX;
            steps = 0;
            for (size_t rr = (size_t)r + 1; rr < Tn && steps < search_cap; ++rr, ++steps) {
                running = std::min(running, (uint32_t)SA.lcp[rr]);
                if (running < K) break;
                uint32_t p = SA.sa[rr];
                if (p < (uint32_t)n) consider_forward(p, running);
                else if (p > (uint32_t)n && p <= (uint32_t)(2 * n)) consider_rc(p - (uint32_t)(n + 1), running);
            }

            if (best_len >= K) bm[i] = BestMatch{best_len, best_src, best_rc};
        }
    };

    if (threads <= 1 || n < 100000) {
        worker(0, n);
    } else {
        std::vector<std::thread> ths;
        ths.reserve((size_t)threads);
        for (int t = 0; t < threads; ++t) {
            size_t lo = n * (size_t)t / (size_t)threads;
            size_t hi = n * (size_t)(t + 1) / (size_t)threads;
            ths.emplace_back(worker, lo, hi);
        }
        for (auto& th : ths) th.join();
    }

    // ── Measurement-only diagnostic (ARCS_REPEAT_ELIM_DP_DEBUG=1): does the
    // greedy left-to-right parse below leave real bytes on the table vs an
    // OPTIMAL parse over the SAME precomputed bm[] candidate set? Classic LZ
    // "optimal parsing" question (same family of local-vs-global-greedy gap
    // already fixed once for assembly's own growth strategy) — bm[] is
    // already computed for EVERY position, so this costs nothing extra to
    // check. Pure side-channel: never touches literal_out/matches_out, never
    // changes real encode behavior, only prints a comparison to stderr.
    if (getenv("ARCS_REPEAT_ELIM_DP_DEBUG")) {
        // Cost-per-unit calibrated from real measured numbers on this
        // project's own yeast_sub.fq run (not arbitrary): literal bytes cost
        // ~0.244 bytes each post-adaptive-coding (1.955 bpb observed on the
        // repeat-elim literal); each match costs ~4.7 compressed bytes
        // (off+len+rc side-stream total / match count, observed).
        const double LIT_COST   = 0.244;
        const double MATCH_COST = 4.7;
        std::vector<double> cost((size_t)n + 1);
        std::vector<uint8_t> take_match((size_t)n, 0);
        cost[n] = 0.0;
        for (size_t ii = n; ii-- > 0; ) {
            double lit_opt = LIT_COST + cost[ii + 1];
            double best = lit_opt;
            bool take = false;
            if (bm[ii].len >= K && ii + bm[ii].len <= n) {
                double match_opt = MATCH_COST + cost[ii + bm[ii].len];
                if (match_opt < best) { best = match_opt; take = true; }
            }
            cost[ii] = best;
            take_match[ii] = take ? 1 : 0;
        }
        size_t dp_literal_bytes = 0, dp_matches = 0;
        {
            size_t ii = 0;
            while (ii < n) {
                if (ii + K <= n && take_match[ii]) { ++dp_matches; ii += bm[ii].len; }
                else { ++dp_literal_bytes; ++ii; }
            }
        }
        size_t greedy_literal_bytes = 0, greedy_matches = 0;
        {
            size_t ii = 0, rs = 0;
            while (ii + K <= n) {
                if (bm[ii].len >= K) { greedy_literal_bytes += (ii - rs); ++greedy_matches; ii += bm[ii].len; rs = ii; }
                else ++ii;
            }
            greedy_literal_bytes += (n - rs);
        }
        fprintf(stderr, "[REPEAT-ELIM-DP] greedy: literal=%zu matches=%zu est_cost=%.0f | optimal: literal=%zu matches=%zu est_cost=%.0f (delta=%.0f bytes)\n",
                greedy_literal_bytes, greedy_matches, greedy_literal_bytes * LIT_COST + greedy_matches * MATCH_COST,
                dp_literal_bytes, dp_matches, cost[0], greedy_literal_bytes * LIT_COST + greedy_matches * MATCH_COST - cost[0]);
    }

    // Serial greedy left-to-right parse: same shape as the previous
    // implementation, just consuming exact precomputed matches instead of
    // hash-bucket-derived approximate ones.
    const char* data = pg.data();
    literal_out.reserve(n);
    size_t i = 0, run_start = 0;
    while (i + K <= n) {
        const BestMatch& m = bm[i];
        if (m.len >= K) {
            size_t run_len = i - run_start;
            literal_out.append(data + run_start, run_len);
            matches_out.push_back({(uint32_t)run_len, m.src, m.len, m.is_rc});
            i += m.len;
            run_start = i;
        } else {
            i++;
        }
    }
    literal_out.append(data + run_start, n - run_start);  // trailing literal run
    return true;
}

std::string repeat_elim_decode(const std::string& literal,
                                const std::vector<RepeatMatch>& matches,
                                uint32_t orig_len) {
    std::string out(orig_len, '\0');
    char*  dst = &out[0];
    const char* lit = literal.data();
    size_t litpos = 0;
    size_t pos    = 0;
    for (const auto& m : matches) {
        memcpy(dst + pos, lit + litpos, m.run_len);
        litpos += m.run_len;
        pos    += m.run_len;
        if (m.is_rc) {
            // dst[pos+t] = complement(dst[m.src + m.len-1-t]) — reverse scan, so no
            // aliasing hazard even though src and dest can be close (unlike the
            // forward memmove case, this always reads strictly right-to-left from
            // a region fully to the left of pos, per the encode-side invariant).
            for (uint32_t t = 0; t < m.len; ++t)
                dst[pos + t] = complement_base(dst[m.src + m.len - 1 - t]);
        } else {
            memmove(dst + pos, dst + m.src, m.len);  // m.src < pos always (encode invariant)
        }
        pos += m.len;
    }
    size_t tail = literal.size() - litpos;
    memcpy(dst + pos, lit + litpos, tail);
    return out;
}
