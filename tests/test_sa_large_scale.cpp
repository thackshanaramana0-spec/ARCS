// Targets the EXACT bug found on real GIAB HG002 data: SuffixArray::build_libsais
// (and build_apsp_candidates' own local size computation) used to compute text
// size via a premature signed-32-bit cast, silently overflowing negative for
// any text beyond INT32_MAX (~2.147B) chars, before ever reaching a size
// check — producing "cannot create std::vector larger than max_size()" deep
// in a downstream vector construction. This test builds a real text PAST that
// exact threshold and verifies the libsais64 dispatch path produces a valid,
// correctly-sorted suffix array with correct LCP values — not just "doesn't
// crash". Full O(n^2) brute-force verification is infeasible at this scale
// (unlike the project's smaller SA/APSP tests), so correctness is checked via:
// (1) SA is a valid permutation of [0,n) (bitmap occupancy check), (2) a large
// random sample of adjacent-rank pairs have their recorded LCP verified
// against the ACTUAL common-prefix length in the text (bounded comparison
// length, since full unbounded comparison at this scale could itself be slow
// in pathological cases — 2000 chars is far more than any real LCP value
// expected in essentially-random synthetic content), and (3) that sampled
// suffix pairs are correctly ordered (SA[i-1] <= SA[i] lexicographically, same
// bounded check).
#include "sa_apsp.h"
#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <random>
#include <string>
#include <thread>
#include <vector>

int main() {
    const size_t TARGET = (size_t)INT32_MAX + (200ull * 1024 * 1024); // ~2.147B + 200MB margin
    printf("[large_scale] building text of %zu chars (target > INT32_MAX=%d)...\n", TARGET, INT32_MAX);

    std::string text;
    text.resize(TARGET);
    {
        std::mt19937_64 rng(42);
        const char alpha[4] = {'A','C','G','T'};
        // Fill in chunks for speed (avoid per-byte RNG call overhead at this scale).
        // Reserve the final byte for the required GSA null terminator (below) —
        // fill only [0, TARGET-1).
        size_t i = 0;
        while (i < TARGET - 1) {
            uint64_t r = rng();
            for (int b = 0; b < 8 && i < TARGET - 1; ++b, r >>= 2) {
                text[i++] = alpha[r & 3];
            }
        }
        // libsais64_gsa (and libsais_gsa) require T[n-1] == 0 — this is a real
        // precondition of the GSA API, already satisfied by the actual pipeline
        // (build_apsp_candidates always appends a trailing SEP='\0' to its last
        // segment too); the test text must honor the same contract.
        text[TARGET - 1] = '\0';
    }
    printf("[large_scale] text built (n=%zu, exceeds INT32_MAX: %s)\n",
           text.size(), text.size() > (size_t)INT32_MAX ? "yes" : "NO -- TEST INVALID");
    if (text.size() <= (size_t)INT32_MAX) { printf("[large_scale] FAILED: test text not actually large enough\n"); return 1; }

    SuffixArray SA;
    int threads = (int)std::thread::hardware_concurrency();
    if (threads < 1) threads = 1;
    printf("[large_scale] calling build_libsais with %d threads...\n", threads);
    try {
        SA.build_libsais(text, threads);
    } catch (const std::exception& e) {
        printf("[large_scale] FAILED: build_libsais threw: %s\n", e.what());
        return 1;
    }
    printf("[large_scale] build_libsais returned successfully (sa.size()=%zu, lcp.size()=%zu)\n",
           SA.sa.size(), SA.lcp.size());

    bool ok = true;
    const size_t n = text.size();
    if (SA.sa.size() != n || SA.lcp.size() != n) {
        printf("[large_scale] FAILED: sa/lcp size mismatch (expected %zu)\n", n);
        ok = false;
    }

    // (1) SA is a valid permutation of [0,n) — bitmap occupancy check (~n/8 bytes,
    // cheap even at this scale: ~268MB for n~2.35B).
    if (ok) {
        std::vector<bool> seen(n, false);
        size_t bad = 0;
        for (size_t i = 0; i < n; ++i) {
            uint32_t v = SA.sa[i];
            if ((size_t)v >= n || seen[v]) { ++bad; if (bad <= 5) printf("[large_scale] bad SA entry at i=%zu: sa=%u\n", i, v); }
            else seen[v] = true;
        }
        if (bad > 0) { printf("[large_scale] FAILED: SA is not a valid permutation (%zu bad entries)\n", bad); ok = false; }
        else printf("[large_scale] SA permutation check: OK\n");
    }

    // (2) + (3) Sampled LCP correctness + SA ordering, bounded comparison length.
    if (ok) {
        std::mt19937_64 rng2(7);
        std::uniform_int_distribution<size_t> dist(1, n - 1);
        const int SAMPLES = 20000;
        const size_t MAXCMP = 2000;
        long lcp_bad = 0, order_bad = 0;
        for (int s = 0; s < SAMPLES; ++s) {
            size_t i = dist(rng2);
            uint32_t a = SA.sa[i - 1], b = SA.sa[i];
            size_t max_len = std::min({MAXCMP, n - (size_t)a, n - (size_t)b});
            size_t real_lcp = 0;
            while (real_lcp < max_len && text[(size_t)a + real_lcp] == text[(size_t)b + real_lcp]) ++real_lcp;
            uint32_t recorded = SA.lcp[i];
            // Only check when the true LCP is within our bounded comparison window
            // (if real_lcp hit MAXCMP without finding a mismatch, the true value
            // could be longer than we measured — skip those rare cases rather than
            // report a false failure).
            if (real_lcp < max_len) {
                if (recorded != real_lcp) {
                    ++lcp_bad;
                    if (lcp_bad <= 5) printf("[large_scale] LCP mismatch at rank %zu: recorded=%u real=%zu\n", i, recorded, real_lcp);
                }
                // Ordering: suffix at a should be <= suffix at b lexicographically
                // up to the compared window (real_lcp chars equal, then either one
                // ran out (shorter is smaller) or a real differing byte at real_lcp).
                if (real_lcp < max_len) {
                    bool a_shorter_or_equal_here = ((size_t)a + real_lcp >= n) || ((size_t)b + real_lcp < n && text[(size_t)a + real_lcp] <= text[(size_t)b + real_lcp]);
                    if (!a_shorter_or_equal_here) { ++order_bad; if (order_bad <= 5) printf("[large_scale] order violation at rank %zu\n", i); }
                }
            }
        }
        printf("[large_scale] sampled %d pairs: lcp_bad=%ld order_bad=%ld\n", SAMPLES, lcp_bad, order_bad);
        if (lcp_bad > 0 || order_bad > 0) { ok = false; }
    }

    printf(ok ? "\n[large_scale] ALL LARGE-SCALE SA CHECKS PASSED\n" : "\n[large_scale] LARGE-SCALE SA CHECK(S) FAILED\n");
    return ok ? 0 : 1;
}
