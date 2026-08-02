#include "bsc_codec.h"
#include "libbsc/libbsc.h"
#include <cstring>
#include <stdexcept>
#include <string>

static bool bsc_initialised = false;
static void ensure_bsc_init() {
    if (!bsc_initialised) {
        int rc = bsc_init(LIBBSC_FEATURE_FASTMODE | LIBBSC_FEATURE_MULTITHREADING);
        if (rc != LIBBSC_NO_ERROR)
            throw std::runtime_error("bsc_init failed: " + std::to_string(rc));
        bsc_initialised = true;
    }
}

std::vector<uint8_t> bsc_compress_buf(const std::vector<uint8_t>& data) {
    ensure_bsc_init();
    if (data.empty()) return {};
    int n = (int)data.size();
    std::vector<uint8_t> out(n + LIBBSC_HEADER_SIZE);
    int rc = bsc_compress(
        data.data(), out.data(), n,
        LIBBSC_DEFAULT_LZPHASHSIZE, LIBBSC_DEFAULT_LZPMINLEN,
        LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_STATIC,
        LIBBSC_FEATURE_FASTMODE | LIBBSC_FEATURE_MULTITHREADING);
    if (rc < 0)
        throw std::runtime_error("bsc_compress failed: " + std::to_string(rc));
    out.resize((size_t)rc);
    return out;
}

std::vector<uint8_t> bsc_decompress_buf(const uint8_t* data, size_t size) {
    ensure_bsc_init();
    if (size == 0) return {};
    int block_size = 0, data_size = 0;
    int rc = bsc_block_info(data, (int)size, &block_size, &data_size,
                            LIBBSC_FEATURE_FASTMODE | LIBBSC_FEATURE_MULTITHREADING);
    if (rc != LIBBSC_NO_ERROR)
        throw std::runtime_error("bsc_block_info failed: " + std::to_string(rc));
    std::vector<uint8_t> out((size_t)data_size);
    rc = bsc_decompress(data, block_size, out.data(), data_size,
                        LIBBSC_FEATURE_FASTMODE | LIBBSC_FEATURE_MULTITHREADING);
    if (rc != LIBBSC_NO_ERROR)
        throw std::runtime_error("bsc_decompress failed: " + std::to_string(rc));
    return out;
}
