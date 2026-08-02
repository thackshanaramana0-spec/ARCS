#pragma once
#include "common.h"
#include <zlib.h>
#include <functional>
#include <memory>
#include <cstring>

// ── Streaming FASTQ reader ────────────────────────────────────────────────────
// Supports: plain FASTQ, gzipped FASTQ (.gz)
// Streaming: reads one record at a time — never loads full file into RAM
// SE and PE modes: for PE, call next() on two readers in lockstep

class FASTQReader {
public:
    explicit FASTQReader(const std::string& path);
    ~FASTQReader();

    // Read next record. Returns false on EOF or error.
    bool next(Read& r);

    // Read up to n records. Returns count read (< n on EOF).
    size_t read_batch(std::vector<Read>& batch, size_t n);

    // Seek to approximate fraction [0,1] of file (for stratified sampling).
    // Only valid for gzipped files if index is built; for plain files, uses fseek.
    // Returns actual fraction seeked to (may differ from requested).
    double seek_fraction(double frac);

    // Number of reads read so far.
    size_t reads_read() const { return reads_read_; }

    bool is_open() const { return open_; }
    const std::string& path() const { return path_; }

private:
    std::string path_;
    bool        gzipped_ = false;
    bool        open_    = false;
    gzFile      gz_      = nullptr;
    FILE*       fp_      = nullptr;
    size_t      reads_read_ = 0;

    // Line buffer
    static constexpr size_t BUF_SIZE = 1 << 20; // 1 MB
    std::vector<char> line_buf_;

    bool read_line(std::string& out);
    bool fill_record(Read& r);
};

// ── Streaming FASTQ writer ────────────────────────────────────────────────────
class FASTQWriter {
public:
    explicit FASTQWriter(const std::string& path);
    ~FASTQWriter();

    void write(const Read& r);
    // final_newline=false omits the trailing line terminator (for a file whose last
    // line has no newline) — used only for the very last record.
    void write(const Read& r, bool final_newline);
    void flush();
    // Emit CRLF ("\r\n") line endings instead of LF, to reproduce Windows-style
    // input byte-exactly. Off by default (standard LF FASTQ).
    void set_crlf(bool on) { crlf_ = on; }

private:
    FILE*   fp_   = nullptr;
    bool    open_ = false;
    bool    crlf_ = false;
    // Large user-space output buffer: records are assembled here via memcpy and
    // flushed to disk in big chunks, avoiding per-record fprintf format parsing and
    // small stdio flushes (the dominant cost of the decoder's write phase).
    std::vector<char> obuf_;
    static constexpr size_t OBUF_FLUSH = 8u << 20;   // flush at 8 MB
    void obuf_flush();
};

// Returns true if the file uses CRLF ("\r\n") line endings (peeks the first line).
bool file_uses_crlf(const std::string& path);
// Returns true if the file's last byte is a newline (usual case).
bool file_ends_with_newline(const std::string& path);

// ── Paired-end reader ─────────────────────────────────────────────────────────
// Reads R1 and R2 in lockstep, yields pairs.
struct ReadPair {
    Read r1, r2;
};

class PEFASTQReader {
public:
    PEFASTQReader(const std::string& path1, const std::string& path2);

    // Read next pair. Returns false on EOF.
    bool next(ReadPair& pair);

    size_t pairs_read() const { return pairs_read_; }

private:
    FASTQReader r1_, r2_;
    size_t pairs_read_ = 0;
};

// ── Stratified sampler ────────────────────────────────────────────────────────
// Reads n_total reads from n_strata equally-spaced file positions.
// Used for AKC pilot and quality model training.
std::vector<Read> sample_stratified(
    const std::string& path,
    size_t n_total,
    size_t n_strata = 3
);
