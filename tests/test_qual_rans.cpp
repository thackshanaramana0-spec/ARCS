#include "../src/rans_model.h"
#include <cstdio>
#include <vector>
#include <random>

// Isolate the exact encode/decode pattern used by encode_quality_rans:
// per-symbol context, encode in reverse, serialize+deserialize, decode forward.
int main() {
    std::mt19937 rng(42);
    int N = 50000;            // symbols
    int alpha = 43;
    int n_ctx = 550;

    // Random symbols and contexts
    std::vector<uint8_t>    syms(N);
    std::vector<ContextKey> ctxs(N);
    for (int i = 0; i < N; ++i) {
        ctxs[i] = rng() % n_ctx;
        syms[i] = rng() % alpha;
    }

    // Train model
    ContextModel m(alpha, 10, 2);
    for (int i = 0; i < N; ++i) m.observe(ctxs[i], syms[i]);
    m.finalize();

    // Encode in reverse
    RansEncoder enc;
    std::vector<uint8_t> out;
    for (int i = N - 1; i >= 0; --i)
        m.encode_sym(out, enc, ctxs[i], syms[i]);
    enc.flush(out);
    std::reverse(out.begin(), out.end());
    printf("encoded %d syms -> %zu bytes\n", N, out.size());

    // Serialize + deserialize
    auto model_bytes = m.serialize();
    ContextModel d(alpha, 10, 2);
    d.deserialize(model_bytes.data(), model_bytes.size());

    // Decode forward
    RansDecoder dec;
    const uint8_t* p = out.data();
    const uint8_t* e = out.data() + out.size();
    dec.init(p); p += 4;

    int mismatches = 0, first_bad = -1;
    for (int i = 0; i < N; ++i) {
        uint8_t q = d.decode_sym(dec, p, e, ctxs[i]);
        if (q != syms[i]) {
            ++mismatches;
            if (first_bad < 0) first_bad = i;
        }
    }
    printf("mismatches: %d  first_bad: %d\n", mismatches, first_bad);

    // Also test WITHOUT serialize round-trip (same model object)
    RansDecoder dec2;
    const uint8_t* p2 = out.data();
    dec2.init(p2); p2 += 4;
    int mm2 = 0, fb2 = -1;
    for (int i = 0; i < N; ++i) {
        uint8_t q = m.decode_sym(dec2, p2, e, ctxs[i]);
        if (q != syms[i]) { ++mm2; if (fb2<0) fb2=i; }
    }
    printf("no-serialize mismatches: %d  first_bad: %d\n", mm2, fb2);

    // Localization: proven free-function block path (single context) at scale
    {
        uint32_t counts[43];
        for (int i=0;i<43;++i) counts[i]=100+i;
        FreqTable ft; ft.build(counts,43);
        std::vector<uint8_t> s2(N);
        for (int i=0;i<N;++i) s2[i]=rng()%43;
        auto e2 = rans_encode_block(s2.data(), N, ft);
        auto d2 = rans_decode_block(e2.data(), e2.size(), ft, N);
        int bm=0; for(int i=0;i<N;++i) if(d2[i]!=s2[i]) ++bm;
        printf("block-path (single ctx) mismatches: %d\n", bm);
    }

    if (mismatches == 0 && mm2 == 0) { printf("PASS\n"); return 0; }
    printf("FAIL\n"); return 1;
}
