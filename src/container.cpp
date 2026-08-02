#include "container.h"
#include <cstring>
#include <zlib.h>
#include <stdexcept>

#ifdef ARCS_HAVE_LZMA
#include <lzma.h>
#endif

#ifdef ARCS_HAVE_ZSTD
#include <zstd.h>
#endif

static constexpr size_t ZSTD_LARGE_THRESHOLD = 65536; // route blobs >=64KB to zstd

// ── Byte order helpers ────────────────────────────────────────────────────────
static void write_be16(uint8_t* p, uint16_t v) {
    p[0] = (v >> 8); p[1] = v;
}
static void write_be32(uint8_t* p, uint32_t v) {
    p[0]=(v>>24); p[1]=(v>>16); p[2]=(v>>8); p[3]=v;
}
static void write_be64(uint8_t* p, uint64_t v) {
    write_be32(p,   (uint32_t)(v >> 32));
    write_be32(p+4, (uint32_t)(v));
}
static uint16_t read_be16(const uint8_t* p) {
    return ((uint16_t)p[0]<<8)|p[1];
}
static uint32_t read_be32(const uint8_t* p) {
    return ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|((uint32_t)p[2]<<8)|p[3];
}
static uint64_t read_be64(const uint8_t* p) {
    return ((uint64_t)read_be32(p) << 32) | read_be32(p+4);
}

// ── LZMA/zlib compression ─────────────────────────────────────────────────────
std::vector<uint8_t> zlib_compress(const std::vector<uint8_t>& data, int level) {
    if (data.empty()) return {};
    uLongf bound = compressBound((uLong)data.size());
    std::vector<uint8_t> out(bound + 4); // 4 bytes for original size header

    // Store original size first (for decompression)
    write_be32(out.data(), (uint32_t)data.size());
    uLongf out_len = bound;
    int ret = compress2(out.data() + 4, &out_len,
                        data.data(), (uLong)data.size(), level);
    ARCS_CHECK(ret == Z_OK, "zlib compress failed");
    out.resize(4 + out_len);
    return out;
}

std::vector<uint8_t> zlib_decompress(const uint8_t* data, size_t len) {
    if (len < 4) return {};
    uint32_t orig_size = read_be32(data);
    std::vector<uint8_t> out(orig_size);
    uLongf out_len = orig_size;
    int ret = uncompress(out.data(), &out_len, data + 4, (uLong)(len - 4));
    ARCS_CHECK(ret == Z_OK, "zlib decompress failed");
    out.resize(out_len);
    return out;
}

std::vector<uint8_t> lzma_compress(const std::vector<uint8_t>& data, int preset) {
#ifdef ARCS_HAVE_LZMA
    if (data.empty()) return {};
    // LZMA encoding: prepend original size (8 bytes)
    std::vector<uint8_t> out(8 + data.size() * 2 + 1024); // worst-case
    lzma_stream strm = LZMA_STREAM_INIT;
    lzma_ret ret = lzma_easy_encoder(&strm, preset, LZMA_CHECK_CRC64);
    ARCS_CHECK(ret == LZMA_OK, "lzma_easy_encoder failed");

    // Store original size
    write_be64(out.data(), data.size());

    strm.next_in   = data.data();
    strm.avail_in  = data.size();
    strm.next_out  = out.data() + 8;
    strm.avail_out = out.size() - 8;

    ret = lzma_code(&strm, LZMA_FINISH);
    lzma_end(&strm);
    ARCS_CHECK(ret == LZMA_STREAM_END, "lzma_code failed");
    out.resize(8 + (out.size() - 8 - strm.avail_out));
    return out;
#else
    // No liblzma: fall back to zlib
    return zlib_compress(data, 9);
#endif
}

std::vector<uint8_t> lzma_decompress(const uint8_t* data, size_t len) {
#ifdef ARCS_HAVE_LZMA
    if (len < 8) return {};
    uint64_t orig_size = read_be64(data);
    std::vector<uint8_t> out(orig_size);

    lzma_stream strm = LZMA_STREAM_INIT;
    lzma_ret ret = lzma_auto_decoder(&strm, UINT64_MAX, 0);
    ARCS_CHECK(ret == LZMA_OK, "lzma_auto_decoder failed");

    strm.next_in   = data + 8;
    strm.avail_in  = len - 8;
    strm.next_out  = out.data();
    strm.avail_out = orig_size;

    ret = lzma_code(&strm, LZMA_FINISH);
    lzma_end(&strm);
    ARCS_CHECK(ret == LZMA_STREAM_END, "lzma_code decomp failed");
    out.resize(orig_size - strm.avail_out);
    return out;
#else
    // No liblzma: fall back to zlib
    return zlib_decompress(data, len);
#endif
}

std::vector<uint8_t> zstd_compress(const std::vector<uint8_t>& data, int level) {
#ifdef ARCS_HAVE_ZSTD
    if (data.empty()) return {};
    size_t bound = ZSTD_compressBound(data.size());
    std::vector<uint8_t> out(8 + bound);
    write_be64(out.data(), (uint64_t)data.size());
    size_t r = ZSTD_compress(out.data() + 8, bound, data.data(), data.size(), level);
    ARCS_CHECK(!ZSTD_isError(r),
               std::string("zstd compress: ") + ZSTD_getErrorName(r));
    out.resize(8 + r);
    return out;
#else
    return lzma_compress(data, 9);
#endif
}

std::vector<uint8_t> zstd_decompress(const uint8_t* data, size_t len) {
#ifdef ARCS_HAVE_ZSTD
    if (len < 8) return {};
    uint64_t orig_size = read_be64(data);
    std::vector<uint8_t> out(orig_size);
    size_t r = ZSTD_decompress(out.data(), orig_size, data + 8, len - 8);
    ARCS_CHECK(!ZSTD_isError(r),
               std::string("zstd decompress: ") + ZSTD_getErrorName(r));
    out.resize(r);
    return out;
#else
    return lzma_decompress(data, len);
#endif
}

// Detect zstd magic at offset 8 (after our 8-byte orig_size header).
// zstd frame magic: 0x28 0xB5 0x2F 0xFD
static bool is_zstd_blob(const uint8_t* data, size_t len) {
    return len >= 12 &&
           data[8] == 0x28 && data[9] == 0xB5 &&
           data[10] == 0x2F && data[11] == 0xFD;
}

std::vector<uint8_t> arcs_compress(const std::vector<uint8_t>& data, int lzma_level) {
#if defined(ARCS_HAVE_ZSTD) && defined(ARCS_HAVE_LZMA)
    // zstd is opt-in (ARCS_COMPRESS_ZSTD=1): use for large blobs only.
    // Default stays LZMA-9 — chain-pg blobs are small enough that LZMA wins on
    // both ratio and speed. zstd may help on very large blob datasets (>1MB blobs).
    if (lzma_level == 9 &&
        data.size() >= ZSTD_LARGE_THRESHOLD &&
        getenv("ARCS_COMPRESS_ZSTD") != nullptr) {
        return zstd_compress(data, 19);
    }
#endif
    return lzma_compress(data, lzma_level);
}

std::vector<uint8_t> arcs_decompress(const uint8_t* data, size_t len) {
    if (is_zstd_blob(data, len))
        return zstd_decompress(data, len);
    return lzma_decompress(data, len);
}

// ── ARCSWriter ────────────────────────────────────────────────────────────────
ARCSWriter::ARCSWriter(const std::string& path) : path_(path) {
    // "w+b" (not "wb") so end_stream_blob can re-read the written region to compute
    // the streaming-blob CRC. Truncates+creates like "wb", but also allows reads.
    fp_ = fopen(path.c_str(), "w+b");
    ARCS_CHECK(fp_ != nullptr, "Cannot create output: " + path);

    // Reserve space for header + blob table (max 12 blobs)
    size_t reserved = HEADER_SIZE + (size_t)BlobType::MAX_BLOB_ID * BLOB_ENTRY_SIZE;
    std::vector<uint8_t> zeros(reserved, 0);
    fwrite(zeros.data(), 1, reserved, fp_);
    data_offset_ = reserved;
}

ARCSWriter::~ARCSWriter() {
    if (fp_) fclose(fp_);
}

size_t ARCSWriter::add_blob(BlobType type, const std::vector<uint8_t>& data) {
    ARCS_CHECK(fp_ != nullptr, "Writer not open");

    fseek(fp_, (long)data_offset_, SEEK_SET);
    fwrite(data.data(), 1, data.size(), fp_);

    uint32_t crc = crc32c(data.data(), data.size());
    blobs_.push_back({type, data_offset_, data.size(), crc});
    data_offset_ += data.size();

    return blobs_.size() - 1;
}

void ARCSWriter::begin_stream_blob() {
    ARCS_CHECK(fp_ != nullptr, "Writer not open");
    stream_start_ = data_offset_;
    fseek(fp_, (long)data_offset_, SEEK_SET);
}

void ARCSWriter::stream_write(const uint8_t* p, size_t n) {
    ARCS_CHECK(fp_ != nullptr, "Writer not open");
    fseek(fp_, (long)data_offset_, SEEK_SET);
    fwrite(p, 1, n, fp_);
    data_offset_ += n;
}

void ARCSWriter::stream_patch(size_t file_offset, const uint8_t* p, size_t n) {
    ARCS_CHECK(fp_ != nullptr, "Writer not open");
    fseek(fp_, (long)file_offset, SEEK_SET);
    fwrite(p, 1, n, fp_);
    fseek(fp_, (long)data_offset_, SEEK_SET);   // restore append position
}

void ARCSWriter::end_stream_blob(BlobType type) {
    ARCS_CHECK(fp_ != nullptr, "Writer not open");
    fflush(fp_);
    size_t total = data_offset_ - stream_start_;
    // CRC over the written region by re-reading in a 64 KB buffer (bounded RAM).
    uint32_t crc = 0xFFFFFFFFu;
    std::vector<uint8_t> buf(1 << 16);
    fseek(fp_, (long)stream_start_, SEEK_SET);
    size_t left = total;
    while (left > 0) {
        size_t want = left < buf.size() ? left : buf.size();
        size_t got = fread(buf.data(), 1, want, fp_);
        ARCS_CHECK(got == want, "stream blob CRC re-read short");
        crc = crc32c(buf.data(), got, crc);
        left -= got;
    }
    fseek(fp_, (long)data_offset_, SEEK_SET);
    blobs_.push_back({type, stream_start_, total, crc});
}

void ARCSWriter::finalize(const ARCSHeader& hdr) {
    ARCS_CHECK(fp_ != nullptr, "Writer not open");

    // Seek back to start and write header + blob table
    fseek(fp_, 0, SEEK_SET);

    // Header
    uint8_t header[HEADER_SIZE] = {0};
    memcpy(header, ARCS_MAGIC.data(), 4);
    header[4] = hdr.version;
    header[5] = hdr.flags;
    write_be16(header + 6, (uint16_t)blobs_.size());
    write_be64(header + 8,  hdr.n_reads);
    write_be64(header + 16, hdr.n_mapped);
    write_be64(header + 24, hdr.genome_len);
    // Float: store as 4 bytes IEEE 754
    memcpy(header + 32, &hdr.akc_score, 4);
    header[36] = hdr.k;
    header[37] = hdr.regime;

    fwrite(header, 1, HEADER_SIZE, fp_);

    // Blob table
    for (const auto& b : blobs_) {
        uint8_t entry[BLOB_ENTRY_SIZE] = {0};
        write_be32(entry,      (uint32_t)b.type);
        write_be64(entry + 4,  b.offset);
        write_be64(entry + 12, b.size);
        write_be32(entry + 20, b.crc32);
        fwrite(entry, 1, BLOB_ENTRY_SIZE, fp_);
    }

    fflush(fp_);
    fclose(fp_);
    fp_ = nullptr;
}

// ── ARCSReader ────────────────────────────────────────────────────────────────
ARCSReader::ARCSReader(const std::string& path) {
    fp_ = fopen(path.c_str(), "rb");
    ARCS_CHECK(fp_ != nullptr, "Cannot open: " + path);
    parse_header();
}

ARCSReader::~ARCSReader() {
    if (fp_) fclose(fp_);
}

void ARCSReader::parse_header() {
    uint8_t header[HEADER_SIZE];
    size_t r = fread(header, 1, HEADER_SIZE, fp_);
    ARCS_CHECK(r == HEADER_SIZE, "Truncated header");

    // Check magic
    ARCS_CHECK(memcmp(header, ARCS_MAGIC.data(), 4) == 0,
               "Not an ARCS file (bad magic)");

    hdr_.version  = header[4];
    hdr_.flags    = header[5];
    uint16_t nb   = read_be16(header + 6);
    hdr_.n_blobs  = nb;
    hdr_.n_reads  = read_be64(header + 8);
    hdr_.n_mapped = read_be64(header + 16);
    hdr_.genome_len = read_be64(header + 24);
    memcpy(&hdr_.akc_score, header + 32, 4);
    hdr_.k      = header[36];
    hdr_.regime = header[37];

    // Read blob table
    blobs_.resize(nb);
    for (uint16_t i = 0; i < nb; ++i) {
        uint8_t entry[BLOB_ENTRY_SIZE];
        r = fread(entry, 1, BLOB_ENTRY_SIZE, fp_);
        ARCS_CHECK(r == BLOB_ENTRY_SIZE, "Truncated blob table");
        blobs_[i].type   = (BlobType)read_be32(entry);
        blobs_[i].offset = read_be64(entry + 4);
        blobs_[i].size   = read_be64(entry + 12);
        blobs_[i].crc32  = read_be32(entry + 20);
    }
}

bool ARCSReader::has_blob(BlobType type) const {
    for (const auto& b : blobs_)
        if (b.type == type) return true;
    return false;
}

std::vector<uint8_t> ARCSReader::read_blob(BlobType type) const {
    for (const auto& b : blobs_) {
        if (b.type != type) continue;

        fseek(fp_, (long)b.offset, SEEK_SET);
        std::vector<uint8_t> data(b.size);
        size_t r = fread(data.data(), 1, b.size, fp_);
        ARCS_CHECK(r == b.size, "Truncated blob data");

        uint32_t crc = crc32c(data.data(), data.size());
        ARCS_CHECK(crc == b.crc32, "Blob CRC32 mismatch — data corruption");

        return data;
    }
    return {};
}

std::vector<std::vector<uint8_t>> ARCSReader::read_all_blobs(BlobType type) const {
    std::vector<std::vector<uint8_t>> result;
    for (const auto& b : blobs_) {
        if (b.type != type) continue;

        fseek(fp_, (long)b.offset, SEEK_SET);
        std::vector<uint8_t> data(b.size);
        size_t r = fread(data.data(), 1, b.size, fp_);
        ARCS_CHECK(r == b.size, "Truncated blob data (read_all_blobs)");

        uint32_t crc = crc32c(data.data(), data.size());
        ARCS_CHECK(crc == b.crc32, "Blob CRC32 mismatch (read_all_blobs)");

        result.push_back(std::move(data));
    }
    return result;
}
