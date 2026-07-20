#pragma once
#include "common.h"
#include <vector>
#include <array>
#include <string>
#include <fstream>

// ── ARCS Container Format v1 ──────────────────────────────────────────────────
// Binary format for fully lossless FASTQ archival.
// All multi-byte integers are big-endian.
//
// Header (fixed 64 bytes):
//   [0..3]   magic: 'A','R','C','S'
//   [4]      version: 1
//   [5]      flags:
//               bit 0: paired_end
//               bit 1: pe_delta_encoding used
//               bit 2: gutensor stored
//               bit 3: names stored (first-token only)
//               bit 4: amplicon path (BSC)
//   [6..7]   n_blobs: number of blobs (uint16_t)
//   [8..15]  n_reads: total read count (uint64_t)
//   [16..23] n_mapped: mapped read count (uint64_t)
//   [24..31] genome_len: pseudogenome length (uint64_t)
//   [32..35] akc_score_float: AKC score (float32)
//   [36]     k: k-mer size used
//   [37]     regime: DataRegime enum
//   [38..63] reserved (0)
//
// Blob table (follows header): n_blobs × 24 bytes each
//   [0..3]   blob_id: BlobType enum (uint32_t)
//   [4..11]  offset: byte offset from file start (uint64_t)
//   [12..19] size: compressed size in bytes (uint64_t)
//   [20..23] crc32: CRC32 of compressed blob data (uint32_t)
//
// Blobs (at their respective offsets):
//   Raw compressed bytes (LZMA or rANS or BSC-encoded)

enum class BlobType : uint32_t {
    GENOME         = 0,  // Pseudogenome G* (LZMA-9 compressed)
    NAMES          = 1,  // First-token read names (LZMA-9 compressed)
    SE_POSITIONS   = 2,  // SE position deltas (varint, then rANS)
    PE_R2_DELTAS   = 3,  // PE R2 position deltas (varint, then rANS)
    MISMATCHES     = 4,  // Mismatch (offset, base) pairs, context rANS
    QUALITY_MODEL  = 5,  // Serialized ContextModel
    QUALITY_DATA   = 6,  // Quality symbols, context rANS
    UNMAPPED       = 7,  // Unmapped reads verbatim FASTQ (LZMA-9)
    GUTENSOR       = 8,  // GU-TENSOR codon token stream
    STRAND_FLAGS   = 9,  // Strand orientation bits (1 per read)
    PE_STATS       = 10, // PEPositionStats (insert mean, std)
    // MST tree encoding blobs (replaces GENOME+SE_POSITIONS+MISMATCHES in MST mode)
    MST_TREE       = 11, // Parent pointer array (compressed)
    MST_DELTAS     = 12, // Per-read delta stream (shift + substitutions)
    MST_RC_FLAGS   = 13, // RC flags for MST-encoded reads
    // V7: amplicon dedup support
    COUNT_DATA     = 14, // uint32_le count array for PCR-dedup unique sequences (LZMA-9)
    // V7*: byte-exact amplicon quality recovery
    QUALITY_PERM   = 15, // uint32_le perm[p]=orig_read_idx at cluster-sorted pos p (LZMA-9)
    // Chunked-mode outer archive: single blob containing a size-prefixed manifest of
    // inner ARCS archives, one per chunk.  Format:
    //   [uint32_be n_chunks]
    //   repeated n_chunks times: [uint64_be chunk_size][chunk_size bytes of inner .arcs]
    CHUNK_DATA     = 17,
    CHAIN_STARTS   = 18, // bitset: bit i=1 if sorted read i is a chain-start (verbatim)
    CHAIN_READ_PERM= 19, // uint32_le[n]: sorted_pos -> original_fastq_index (LZMA-9)
    MST_TREE_DFS   = 20, // DFS-position parent back-deltas: varint(dfs_pos - parent_dfs_pos)
    CHAIN_VSEQ     = 21, // chain verbatim sequences + N-pos (LZMA-9, separate from delta fields)
    CHAIN_PG_SEQ   = 22, // chain-pg pseudogenome: [1B codec(0x01=LZMA,0x03=GeCo3)][payload]
    CHAIN_PG_POS   = 23, // chain-pg pg positions: LZMA-9 of uint32_be[n_reads] (chain order)
    CHAIN_PG_AUX   = 24, // chain-pg RC+mm+N+qmm data: LZMA-9 of columnar binary (chain order)
    PLUS_LINES     = 25, // arbitrary '+' line suffixes (assembly order, newline-joined, LZMA-9)
    MAX_BLOB_ID    = 16  // reserved blob-table slots (chunked outer archive uses only 1)
};

static constexpr size_t HEADER_SIZE    = 64;
static constexpr size_t BLOB_ENTRY_SIZE = 24;
static constexpr std::array<char,4> ARCS_MAGIC = {'A','R','C','S'};

struct BlobEntry {
    BlobType type   = BlobType::MAX_BLOB_ID;
    uint64_t offset = 0;
    uint64_t size   = 0;
    uint32_t crc32  = 0;
};

struct ARCSHeader {
    uint8_t  version    = 1;
    uint8_t  flags      = 0;
    uint16_t n_blobs    = 0;
    uint64_t n_reads    = 0;
    uint64_t n_mapped   = 0;
    uint64_t genome_len = 0;
    float    akc_score  = 0.0f;
    uint8_t  k          = 0;
    uint8_t  regime     = 0;
    uint16_t read_len   = 0;  // base read length (stored in reserved bytes [38..39])
};

// ── Container writer ──────────────────────────────────────────────────────────
class ARCSWriter {
public:
    explicit ARCSWriter(const std::string& path);
    ~ARCSWriter();

    // Add a blob. data is raw compressed bytes.
    // Returns blob index.
    size_t add_blob(BlobType type, const std::vector<uint8_t>& data);

    // ── Streaming blob (bounded RAM) ───────────────────────────────────────────
    // Write a blob incrementally so the whole payload never resides in RAM (used by
    // chunked mode, where the manifest can reach GB scale on whole-genome files).
    // begin → repeated stream_write → (optional stream_patch) → end_stream_blob.
    // The blob CRC is computed by re-reading the written region at end (sequential,
    // 64 KB buffer), so any in-place patch is included.
    void   begin_stream_blob();
    void   stream_write(const uint8_t* p, size_t n);
    size_t stream_offset() const { return data_offset_; }   // current file write offset
    void   stream_patch(size_t file_offset, const uint8_t* p, size_t n);  // overwrite in place
    void   end_stream_blob(BlobType type);

    // Finalize: write header + blob table at start of file.
    void finalize(const ARCSHeader& hdr);

private:
    std::string          path_;
    FILE*                fp_    = nullptr;
    std::vector<BlobEntry> blobs_;
    size_t               data_offset_ = 0; // current write position
    size_t               stream_start_ = 0; // start offset of an open streaming blob
};

// ── Container reader ──────────────────────────────────────────────────────────
class ARCSReader {
public:
    explicit ARCSReader(const std::string& path);
    ~ARCSReader();

    const ARCSHeader& header() const { return hdr_; }

    // Check if a blob exists.
    bool has_blob(BlobType type) const;

    // Read blob data. Returns raw compressed bytes (first match).
    std::vector<uint8_t> read_blob(BlobType type) const;

    // Read ALL blobs of a given type in table order.
    // Used by chunked-mode decoder to iterate CHUNK_DATA blobs.
    std::vector<std::vector<uint8_t>> read_all_blobs(BlobType type) const;

private:
    FILE*                  fp_ = nullptr;
    ARCSHeader             hdr_;
    std::vector<BlobEntry> blobs_;

    void parse_header();
};

// ── LZMA compression wrappers ─────────────────────────────────────────────────
// These compress/decompress blobs before storage.
// If LZMA not available (no liblzma), falls back to zlib-9.
std::vector<uint8_t> lzma_compress(const std::vector<uint8_t>& data, int preset = 9);
std::vector<uint8_t> lzma_decompress(const uint8_t* data, size_t len);

// zlib fallback (always available)
std::vector<uint8_t> zlib_compress(const std::vector<uint8_t>& data, int level = 9);
std::vector<uint8_t> zlib_decompress(const uint8_t* data, size_t len);
