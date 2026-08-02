#pragma once
#include <cstdint>
#include <vector>
#include <stdexcept>

// Buffer-level BSC compress/decompress for MST delta/tree blobs.
// Uses libbsc (BWT + QLFC) which outperforms LZMA on structured delta streams.

std::vector<uint8_t> bsc_compress_buf(const std::vector<uint8_t>& data);
std::vector<uint8_t> bsc_decompress_buf(const uint8_t* data, size_t size);
