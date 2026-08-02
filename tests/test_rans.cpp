#include "../third_party/rans_byte.h"
#include "../src/rans_model.h"
#include <cassert>
#include <cstdio>
#include <vector>

int main() {
    // Test FreqTable encode/decode
    uint32_t counts[4] = {100, 50, 30, 20}; // A, C, G, T
    FreqTable ft;
    ft.build(counts, 4);

    // All frequencies must sum to SCALE
    uint32_t sum = 0;
    for (int i = 0; i < 4; ++i) sum += ft.syms[i].freq;
    assert(sum == SCALE);

    // Test encode + decode round-trip
    std::vector<uint8_t> symbols = {0, 1, 2, 3, 0, 0, 1, 2, 3, 1};
    auto encoded = rans_encode_block(symbols.data(), symbols.size(), ft);
    assert(!encoded.empty());
    auto decoded = rans_decode_block(encoded.data(), encoded.size(), ft, symbols.size());
    assert(decoded.size() == symbols.size());
    for (size_t i = 0; i < symbols.size(); ++i)
        assert(decoded[i] == symbols[i]);

    // Test ContextModel
    ContextModel model(5, 10, 4); // 5-symbol alphabet (ACGTN), 10 pos bins
    // Observe some data
    for (int ctx = 0; ctx < 40; ++ctx)
        for (uint8_t sym = 0; sym < 5; ++sym)
            for (int cnt = 0; cnt < 10; ++cnt)
                model.observe((ContextKey)ctx, sym);
    model.finalize();

    // Encode sequence with context model
    std::vector<uint8_t> seq = {0, 1, 2, 3, 0, 1, 2, 3};
    auto enc_fn = [](size_t i, const std::vector<uint8_t>&) -> ContextKey {
        return (ContextKey)(i % 40);
    };
    auto encoded2 = model.encode_sequence(seq, enc_fn);
    assert(!encoded2.empty());

    auto decoded2 = model.decode_sequence(encoded2.data(), encoded2.size(),
                                          seq.size(), enc_fn);
    assert(decoded2.size() == seq.size());
    for (size_t i = 0; i < seq.size(); ++i)
        assert(decoded2[i] == seq[i]);

    printf("PASS: test_rans (encoded=%zu bytes for %zu symbols)\n",
           encoded.size(), symbols.size());
    return 0;
}
