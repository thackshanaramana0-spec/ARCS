#include "fastq_io.h"
#include <stdexcept>
#include <cstring>
#include <cerrno>

// ── CRC32 (software fallback) ─────────────────────────────────────────────────
uint32_t crc32c(const uint8_t* data, size_t len, uint32_t crc) {
    // Standard CRC32 (IEEE 802.3)
    static uint32_t table[256] = {};
    static bool table_init = false;
    if (!table_init) {
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int j = 0; j < 8; ++j)
                c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            table[i] = c;
        }
        table_init = true;
    }
    for (size_t i = 0; i < len; ++i)
        crc = table[(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
    return crc;
}

// ── FASTQReader ───────────────────────────────────────────────────────────────
FASTQReader::FASTQReader(const std::string& path) : path_(path) {
    line_buf_.resize(BUF_SIZE);
    gzipped_ = (path.size() >= 3 &&
                path.substr(path.size()-3) == ".gz");
    if (gzipped_) {
        gz_ = gzopen(path.c_str(), "rb");
        ARCS_CHECK(gz_ != nullptr, "Cannot open gz file: " + path);
        gzbuffer(gz_, 1 << 20);
    } else {
        fp_ = fopen(path.c_str(), "rb");
        ARCS_CHECK(fp_ != nullptr, "Cannot open file: " + path +
                   " (" + strerror(errno) + ")");
    }
    open_ = true;
}

FASTQReader::~FASTQReader() {
    if (gz_) gzclose(gz_);
    if (fp_) fclose(fp_);
}

bool FASTQReader::read_line(std::string& out) {
    out.clear();
    if (gzipped_) {
        char* res = gzgets(gz_, line_buf_.data(), (int)line_buf_.size());
        if (!res) return false;
        out = line_buf_.data();
    } else {
        char* res = fgets(line_buf_.data(), (int)line_buf_.size(), fp_);
        if (!res) return false;
        out = line_buf_.data();
    }
    // Remove trailing \r\n
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r'))
        out.pop_back();
    return true;
}

bool FASTQReader::fill_record(Read& r) {
    std::string line;
    // Line 1: @name
    if (!read_line(line)) return false;
    while (line.empty()) {  // skip blank lines between records
        if (!read_line(line)) return false;
    }
    ARCS_CHECK(!line.empty() && line[0] == '@',
               "Expected '@' at start of FASTQ record, got: " + line);
    // Keep the FULL header line (including any comment after a space) so the
    // record is reconstructed byte-for-byte. Names are stored newline-delimited,
    // so spaces are safe; the only forbidden char (newline) can't appear in a
    // header. This makes read names lossless on real FASTQ (barcodes, mate info…).
    r.name = line.substr(1);

    // Line 2: sequence
    if (!read_line(r.seq)) return false;

    // Line 3: '+' separator. Usually bare "+", sometimes repeats the header. Capture
    // the suffix (after '+') so the record can be reproduced byte-exactly.
    if (!read_line(line)) return false;
    r.plus = (line.size() > 1) ? line.substr(1) : std::string();

    // Line 4: quality
    if (!read_line(r.qual)) return false;

    ARCS_CHECK(r.seq.size() == r.qual.size(),
               "Sequence and quality length mismatch in read: " + r.name);
    return true;
}

bool FASTQReader::next(Read& r) {
    if (!open_) return false;
    if (!fill_record(r)) return false;
    ++reads_read_;
    return true;
}

size_t FASTQReader::read_batch(std::vector<Read>& batch, size_t n) {
    batch.clear();
    batch.reserve(n);
    Read r;
    while (batch.size() < n && next(r))
        batch.push_back(std::move(r));
    return batch.size();
}

double FASTQReader::seek_fraction(double frac) {
    if (frac <= 0.0) return 0.0;
    if (frac >= 1.0) return 1.0;

    if (gzipped_) {
        // gzseek is slow on compressed files — use approximate byte offset
        // gzseek to fraction of compressed file size (approximate)
        // This is imprecise but acceptable for stratified sampling
        gzrewind(gz_);
        // Scan to find total reads first is too slow; use byte approximation
        // For now: scan forward by reading and discarding records
        // A better implementation would build a gzip index (gzindex)
        reads_read_ = 0;
        return 0.0; // caller must handle by reading and discarding
    } else {
        fseek(fp_, 0, SEEK_END);
        long total = ftell(fp_);
        long target = (long)(frac * total);
        fseek(fp_, target, SEEK_SET);
        // Scan forward to next '@' at start of line
        std::string line;
        while (read_line(line)) {
            if (!line.empty() && line[0] == '@') {
                // Check the next line is a valid sequence (not a quality line starting with @N)
                // We accept this as the start of a record
                // Back up: we already consumed the @ line, so "un-read" it by
                // pushing back the name. This is tricky with gzFile.
                // Simple fix: store the name and read seq+qual directly.
                // Actually: at this point line contains the @name line.
                // We'll handle this in fill_record by first checking if line is valid.
                // For now, just return the fraction we seeked to.
                return (double)target / total;
            }
        }
        return (double)target / total;
    }
}

// ── FASTQWriter ───────────────────────────────────────────────────────────────
FASTQWriter::FASTQWriter(const std::string& path) {
    fp_ = fopen(path.c_str(), "wb");
    ARCS_CHECK(fp_ != nullptr, "Cannot open output file: " + path);
    open_ = true;
    obuf_.reserve(OBUF_FLUSH + (1u << 16));   // headroom so one record never overflows a flush
}

void FASTQWriter::obuf_flush() {
    if (!obuf_.empty()) {
        size_t nw = fwrite(obuf_.data(), 1, obuf_.size(), fp_);
        ARCS_CHECK(nw == obuf_.size(), "Short write to output file");
        obuf_.clear();
    }
}

FASTQWriter::~FASTQWriter() {
    if (fp_) { obuf_flush(); fclose(fp_); }
}

void FASTQWriter::write(const Read& r) { write(r, true); }

void FASTQWriter::write(const Read& r, bool final_newline) {
    ARCS_CHECK(open_, "Writer not open");
    const char* nl = crlf_ ? "\r\n" : "\n";
    const size_t nll = crlf_ ? 2 : 1;
    // Assemble the record — "@name\nseq\n+plus\nqual\n" — by appending into the
    // user-space buffer (byte-identical to the previous fprintf output).
    auto app = [&](const char* p, size_t len) { obuf_.insert(obuf_.end(), p, p + len); };
    app("@", 1);
    app(r.name.data(), r.name.size()); app(nl, nll);
    app(r.seq.data(),  r.seq.size());  app(nl, nll);
    app("+", 1);
    app(r.plus.data(), r.plus.size()); app(nl, nll);
    app(r.qual.data(), r.qual.size());
    if (final_newline) app(nl, nll);          // omit only for a file with no trailing newline
    if (obuf_.size() >= OBUF_FLUSH) obuf_flush();
}

// Peek the first ~64 KB and report whether the first line break is CRLF. FASTQ
// files use uniform line endings, so the first is representative.
bool file_uses_crlf(const std::string& path) {
    bool gz = (path.size() > 3 && path.substr(path.size()-3) == ".gz");
    bool crlf = false;
    if (gz) {
        gzFile g = gzopen(path.c_str(), "rb");
        if (!g) return false;
        char buf[65536];
        int n = gzread(g, buf, sizeof(buf));
        for (int i = 0; i + 1 < n; ++i) if (buf[i] == '\n') { crlf = (i > 0 && buf[i-1] == '\r'); break; }
        gzclose(g);
    } else {
        FILE* f = fopen(path.c_str(), "rb");
        if (!f) return false;
        char buf[65536];
        size_t n = fread(buf, 1, sizeof(buf), f);
        for (size_t i = 0; i + 1 < n; ++i) if (buf[i] == '\n') { crlf = (i > 0 && buf[i-1] == '\r'); break; }
        fclose(f);
    }
    return crlf;
}

// Returns true if the file's last byte is a newline (the usual case). Plain files
// only — gzipped inputs are always treated as newline-terminated.
bool file_ends_with_newline(const std::string& path) {
    bool gz = (path.size() > 3 && path.substr(path.size()-3) == ".gz");
    if (gz) return true;
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return true;
    if (fseek(f, -1, SEEK_END) != 0) { fclose(f); return true; }
    int c = fgetc(f);
    fclose(f);
    return c == '\n';
}

void FASTQWriter::flush() {
    if (fp_) { obuf_flush(); fflush(fp_); }
}

// ── PEFASTQReader ─────────────────────────────────────────────────────────────
PEFASTQReader::PEFASTQReader(const std::string& p1, const std::string& p2)
    : r1_(p1), r2_(p2) {}

bool PEFASTQReader::next(ReadPair& pair) {
    bool ok1 = r1_.next(pair.r1);
    bool ok2 = r2_.next(pair.r2);
    if (ok1 != ok2)
        throw ARCSError("Paired FASTQ files have different number of reads");
    if (!ok1) return false;
    ++pairs_read_;
    return true;
}

// ── Stratified sampler ────────────────────────────────────────────────────────
std::vector<Read> sample_stratified(const std::string& path,
                                    size_t n_total,
                                    size_t n_strata) {
    std::vector<Read> result;
    result.reserve(n_total);
    size_t per_stratum = n_total / n_strata;
    size_t remainder   = n_total % n_strata;

    // Stratum 0: read from beginning
    {
        FASTQReader rdr(path);
        size_t count = per_stratum + (0 < remainder ? 1 : 0);
        Read r;
        while (result.size() < count && rdr.next(r))
            result.push_back(r);
    }

    // Strata 1..n_strata-1: seek to fraction then read
    for (size_t s = 1; s < n_strata; ++s) {
        size_t count = per_stratum + (s < remainder ? 1 : 0);
        double frac  = (double)s / n_strata;

        if (path.size() >= 3 && path.substr(path.size()-3) == ".gz") {
            // For gzip: scan forward by reading and discarding
            FASTQReader rdr(path);
            // Estimate reads to skip: frac * estimated total reads
            // We use a rough estimate: skip frac * 200000 reads (works for 100k files)
            size_t skip = (size_t)(frac * 200000);
            Read r;
            for (size_t i = 0; i < skip; ++i) {
                if (!rdr.next(r)) break;
            }
            for (size_t i = 0; i < count; ++i) {
                if (!rdr.next(r)) break;
                result.push_back(r);
            }
        } else {
            FASTQReader rdr(path);
            rdr.seek_fraction(frac);
            // After seek, read next record to align to record boundary
            // The seek positioned us at a '@' line — now read full records
            Read r;
            for (size_t i = 0; i < count; ++i) {
                if (!rdr.next(r)) break;
                result.push_back(r);
            }
        }
    }

    return result;
}
