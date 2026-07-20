// Direct coder microbenchmark: times ONLY dna_encode on a large pg, isolating
// the per-base model cost from the rest of the pipeline. ARCSDNA_ORDERS (env,
// read once by the coder's load_cfg) selects the model set — run one config per
// process. Reports encode wall time + output size. Reads/chain_order/pg_pos are
// empty (pre-seeding is optional; we measure the coder's own hot loop).
#include "dna_coder.h"
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <chrono>

int main(int argc, char** argv) {
    size_t N = (argc > 1) ? (size_t)atol(argv[1]) : 5'000'000;  // pg length in bases
    // Build a pg with realistic local structure: order-1 Markov-ish + periodic
    // repeats, so high-order models have something to latch onto (not white noise,
    // which would make all models useless and hide their cost differences).
    std::string pg; pg.reserve(N);
    static const char B[4] = {'A','C','G','T'};
    uint64_t x = 0x9E3779B97F4A7C15ULL;
    std::string motif;
    for (int i = 0; i < 200; ++i) { x = x*6364136223846793005ULL + 1; motif += B[(x>>60)&3]; }
    for (size_t i = 0; i < N; ++i) {
        x = x*6364136223846793005ULL + 1442695040888963407ULL;
        // 70% follow a repeating motif (long-range structure), 30% random
        if ((x & 0x3FF) < 720) pg += motif[i % motif.size()];
        else                   pg += B[(x >> 61) & 3];
    }
    std::vector<std::string> reads;             // empty: measure coder alone
    std::vector<uint32_t> chain_order, pg_pos;

    // one timed encode (plus a warm-up discarded)
    auto out0 = dna_encode(pg, reads, chain_order, pg_pos);
    double best = 1e9; size_t osz = out0.size();
    for (int r = 0; r < 3; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        auto out = dna_encode(pg, reads, chain_order, pg_pos);
        auto t1 = std::chrono::steady_clock::now();
        double s = std::chrono::duration<double>(t1 - t0).count();
        if (s < best) best = s;
        osz = out.size();
    }
    const char* ord = getenv("ARCSDNA_ORDERS");
    printf("orders=%-28s pg=%zu  encode=%.3fs (best of 3)  out=%zu B  %.3f bpb\n",
           ord ? ord : "(default)", N, best, osz, 8.0 * osz / N);
    return 0;
}
