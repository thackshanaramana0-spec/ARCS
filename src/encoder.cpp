#include "encoder.h"
#include "chain_encoder.h"
#include "dna_coder.h"
#include "qual_cm.h"
#include "name_num_codec.h"
#include <chrono>
#include <future>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <cstdlib>

#ifdef _WIN32
// POSIX setenv/unsetenv are absent on MinGW; provide portable shims.
static inline int setenv(const char* name, const char* value, int /*overwrite*/) {
    return _putenv_s(name, value);
}
static inline int unsetenv(const char* name) {
    return _putenv_s(name, "");
}
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
// Peak working-set (high-water RAM) in MB. K32GetProcessMemoryInfo is exported by
// kernel32 (auto-linked) so this needs no extra link deps.
static size_t cur_peak_rss_mb() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (K32GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return (size_t)(pmc.PeakWorkingSetSize / (1024 * 1024));
    return 0;
}
#else
#include <cstdio>
// Peak RSS (VmHWM) in MB from /proc/self/status.
static size_t cur_peak_rss_mb() {
    FILE* f = fopen("/proc/self/status", "r"); if (!f) return 0;
    char line[256]; size_t kb = 0;
    while (fgets(line, sizeof(line), f)) if (sscanf(line, "VmHWM: %zu kB", &kb) == 1) break;
    fclose(f); return kb / 1024;
}
#endif

ARCSEncoder::ARCSEncoder(const ARCSParams& params) : params_(params) {}

// ── Load reads ────────────────────────────────────────────────────────────────
std::vector<Read> ARCSEncoder::load_reads(const std::string& path) const {
    std::vector<Read> reads;
    reads.reserve(200000);
    FASTQReader rdr(path);
    Read r;
    while (rdr.next(r))
        reads.push_back(std::move(r));
    return reads;
}

// ── Train context models on pilot reads ──────────────────────────────────────
void ARCSEncoder::train_models(const std::vector<Read>& reads,
                                const std::vector<MapResult>& mappings,
                                const std::string& genome,
                                ContextModel& mm_model,
                                ContextModel& q_model) const {
    size_t pilot = std::min(reads.size(), (size_t)params_.pilot);

    for (size_t i = 0; i < pilot; ++i) {
        const auto& r   = reads[i];
        const auto& m   = mappings[i];
        int rlen = (int)r.seq.size();
        if (!m.mapped) continue;

        // GC content of read
        int gc = 0;
        for (char c : r.seq) if (c == 'G' || c == 'C') ++gc;
        float gc_frac = (rlen > 0) ? (float)gc / rlen : 0.5f;

        // Build deviation set (which positions are mismatches)
        std::vector<bool> is_dev(rlen, false);
        for (uint32_t off : m.mm_offsets)
            if ((int)off < rlen) is_dev[off] = true;

        // Observe mismatches
        for (size_t j = 0; j < m.mm_offsets.size(); ++j) {
            int off = (int)m.mm_offsets[j];
            if (off >= rlen) continue;
            ContextKey ctx = mismatch_ctx(off, rlen, gc_frac);
            uint8_t ref_base = m.mm_bases[j];
            mm_model.observe(ctx, ref_base);
        }

        // Observe quality symbols
        uint8_t q_prev = 30; // default starting quality
        for (int j = 0; j < rlen; ++j) {
            if (r.qual[j] < 33) continue; // skip invalid Phred
            uint8_t q = (uint8_t)(r.qual[j] - 33); // Phred+33 → raw score
            ContextKey ctx = quality_ctx(q_prev, is_dev[j], j, rlen);
            q_model.observe(ctx, std::min(q, (uint8_t)42));
            q_prev = q;
        }
    }

    mm_model.finalize();
    q_model.finalize();
}

// ── Encode mismatches ─────────────────────────────────────────────────────────
// Format per mapped read: [n_mm:varint] [offset:varint]* [base:rANS]*
std::vector<uint8_t> ARCSEncoder::encode_mismatches(
    const std::vector<Read>& reads,
    const std::vector<MapResult>& mappings,
    const std::string& genome,
    ContextModel& mm_model) const {

    std::vector<uint8_t> out;
    out.reserve(reads.size() * 4);

    for (size_t i = 0; i < reads.size(); ++i) {
        const auto& r = reads[i];
        const auto& m = mappings[i];
        if (!m.mapped) {
            write_varint(out, 0); // 0 mismatches = unmapped marker
            continue;
        }

        int rlen = (int)r.seq.size();
        write_varint(out, m.mm_offsets.size());

        // Encode offsets as deltas
        uint32_t prev_off = 0;
        for (uint32_t off : m.mm_offsets) {
            write_varint(out, off - prev_off);
            prev_off = off;
        }

        // Encode base bytes with context rANS (skip entirely if 0 mismatches)
        if (!m.mm_offsets.empty()) {
            int gc = 0;
            for (char c : r.seq) if (c == 'G' || c == 'C') ++gc;
            float gc_frac = (rlen > 0) ? (float)gc / rlen : 0.5f;

            RansEncoder enc;
            std::vector<uint8_t> base_bytes;
            for (int j = (int)m.mm_offsets.size()-1; j >= 0; --j) {
                int off = (int)m.mm_offsets[j];
                uint8_t base = m.mm_bases[j];
                ContextKey ctx = mismatch_ctx(off, rlen, gc_frac);
                mm_model.encode_sym(base_bytes, enc, ctx, base);
            }
            enc.flush(base_bytes);
            std::reverse(base_bytes.begin(), base_bytes.end());
            write_varint(out, base_bytes.size());
            out.insert(out.end(), base_bytes.begin(), base_bytes.end());
        } else {
            write_varint(out, 0); // no base bytes to store
        }
    }
    return out;
}

// ── Encode quality ────────────────────────────────────────────────────────────
std::vector<uint8_t> ARCSEncoder::encode_quality(
    const std::vector<Read>& reads,
    const std::vector<MapResult>& mappings,
    ContextModel& q_model) const {

    std::vector<uint8_t> out;
    out.reserve(reads.size() * reads[0].qual.size());

    for (size_t i = 0; i < reads.size(); ++i) {
        const auto& r = reads[i];
        const auto& m = mappings[i];
        int rlen = (int)r.qual.size();

        std::vector<bool> is_dev(rlen, false);
        if (m.mapped)
            for (uint32_t off : m.mm_offsets)
                if ((int)off < rlen) is_dev[off] = true;

        // Encode quality symbols with context rANS
        std::vector<uint8_t> qual_syms(rlen);
        for (int j = 0; j < rlen; ++j) {
            uint8_t q = (r.qual[j] >= 33) ? (uint8_t)(r.qual[j] - 33) : 0;
            qual_syms[j] = std::min(q, (uint8_t)42);
        }

        auto encode_fn = [&](size_t j, const std::vector<uint8_t>& decoded) -> ContextKey {
            uint8_t q_prev = (j > 0) ? decoded[j-1] : 30;
            return quality_ctx(q_prev, is_dev[j], (int)j, rlen);
        };

        auto encoded = q_model.encode_sequence(qual_syms, encode_fn);
        write_varint(out, encoded.size());
        out.insert(out.end(), encoded.begin(), encoded.end());
    }
    return out;
}

// ── Columnar name tokenizer (format 0x04) ────────────────────────────────────
// Split every name into maximal runs of digits / non-digits, giving a fixed column
// layout (e.g. "ERR11820348.1 HWI-ST486:315:D10JWACXX:6:1101:2370:1888" → literal
// "ERR", number 11820348, literal ".", number 1, literal " HWI-ST", …). Then per
// column: a CONSTANT literal is stored once; a NUMERIC column is delta+zigzag+varint
// (so a sequential read-index 1,2,3,… collapses to ~0 bytes); a varying literal is
// stored NUL-separated. The whole columnar buffer is LZMA-9'd. This generalises to
// ANY structured name (SRA, Illumina, …) — not just the Illumina :X:Y case — and is
// where SPRING/fqzcomp/DSRC get their read-ID compression. Keep-smaller gated, so it
// can never regress. Returns empty if names don't share a uniform column structure.
static void arcs_tokenize(const std::string& s,
                          std::vector<std::pair<bool,std::string>>& toks) {
    toks.clear();
    size_t i = 0, n = s.size();
    while (i < n) {
        bool dig = (s[i] >= '0' && s[i] <= '9');
        size_t j = i + 1;
        while (j < n && ((s[j] >= '0' && s[j] <= '9') == dig)) ++j;
        toks.emplace_back(dig, s.substr(i, j - i));
        i = j;
    }
}
// Find a numeric name column whose value is STRICTLY INCREASING across reads in
// file order (e.g. the SRA sequential accession index ERR….1,.2,.3…). If present,
// the file is sorted by that column, so read order is recoverable at decode by
// parsing+ranking that column from the names — no explicit permutation needed.
// Returns the column index, or -1. Tokenization matches arcs_tokenize (digit runs).
static int find_monotonic_name_index_column(const std::vector<Read>& reads) {
    size_t n = reads.size();
    if (n < 2) return -1;
    std::vector<std::pair<bool,std::string>> t0, tk;
    arcs_tokenize(reads[0].name, t0);
    size_t ncols = t0.size();
    if (ncols == 0 || ncols > 255) return -1;
    // Candidate numeric columns from read 0 (no leading zeros, ≤18 digits).
    std::vector<int> cand;
    for (size_t c = 0; c < ncols; ++c) {
        if (!t0[c].first) continue;
        const std::string& v = t0[c].second;
        if ((v.size() > 1 && v[0] == '0') || v.size() > 18) continue;
        cand.push_back((int)c);
    }
    if (cand.empty()) return -1;
    std::vector<int64_t> prev(ncols, INT64_MIN);
    for (int c : cand) prev[(size_t)c] = INT64_MIN;
    std::vector<bool> alive(ncols, false);
    for (int c : cand) alive[(size_t)c] = true;
    for (size_t i = 0; i < n; ++i) {
        arcs_tokenize(reads[i].name, tk);
        if (tk.size() != ncols) return -1;
        bool any = false;
        for (int c : cand) {
            if (!alive[(size_t)c]) continue;
            const auto& tokc = tk[(size_t)c];
            if (!tokc.first || (tokc.second.size() > 1 && tokc.second[0] == '0') || tokc.second.size() > 18) {
                alive[(size_t)c] = false; continue;
            }
            int64_t x = (int64_t)strtoll(tokc.second.c_str(), nullptr, 10);
            if (x <= prev[(size_t)c]) alive[(size_t)c] = false;      // not strictly increasing
            else                      prev[(size_t)c] = x;
            any = any || alive[(size_t)c];
        }
        if (!any) return -1;
    }
    for (int c : cand) if (alive[(size_t)c]) return c;
    return -1;
}

static std::vector<uint8_t> build_columnar_names(const std::vector<Read>& reads) {
    size_t n = reads.size();
    if (n < 2) return {};
    std::vector<std::pair<bool,std::string>> t0;
    arcs_tokenize(reads[0].name, t0);
    size_t ncols = t0.size();
    if (ncols == 0 || ncols > 255) return {};
    std::vector<bool> isdig(ncols);
    for (size_t c = 0; c < ncols; ++c) isdig[c] = t0[c].first;
    std::vector<std::vector<std::string>> col(ncols);
    for (size_t c = 0; c < ncols; ++c) col[c].reserve(n);
    std::vector<std::pair<bool,std::string>> tk;
    for (size_t i = 0; i < n; ++i) {
        arcs_tokenize(reads[i].name, tk);
        if (tk.size() != ncols) return {};                 // non-uniform → give up
        for (size_t c = 0; c < ncols; ++c) {
            if (tk[c].first != isdig[c]) return {};
            col[c].push_back(std::move(tk[c].second));
        }
    }
    // Pre-parse numeric columns once (needed for order-key detection + encoding).
    std::vector<bool> is_numeric(ncols, false);
    std::vector<std::vector<int64_t>> num(ncols);
    for (size_t c = 0; c < ncols; ++c) {
        if (!isdig[c]) continue;
        bool ok = true;
        for (auto& v : col[c]) if ((v.size() > 1 && v[0] == '0') || v.size() > 18) { ok = false; break; }
        if (!ok) continue;
        is_numeric[c] = true;
        num[c].resize(n);
        for (size_t i = 0; i < n; ++i) num[c][i] = (int64_t)strtoll(col[c][i].c_str(), nullptr, 10);
    }
    // ORDER KEY: a numeric column that is a permutation of [1..n] (the sequential
    // read-index). In it, the records' original FILE order is recoverable. Any
    // spatially-local field (flowcell X/Y) that was scrambled by the assembly's
    // reordering becomes local again when encoded in that file order → small deltas.
    int order_key = -1;
    for (size_t c = 0; c < ncols && order_key < 0; ++c) {
        if (!is_numeric[c]) continue;
        int64_t mn = num[c][0], mx = num[c][0];
        for (int64_t x : num[c]) { if (x < mn) mn = x; if (x > mx) mx = x; }
        if (mn != 1 || mx != (int64_t)n) continue;
        std::vector<char> seen(n + 1, 0); bool perm = true;
        for (int64_t x : num[c]) { if (x < 1 || x > (int64_t)n || seen[(size_t)x]) { perm = false; break; } seen[(size_t)x] = 1; }
        if (perm) order_key = (int)c;
    }
    // file-order index: fileorder[k] = record whose key == k+1  (only if order_key set)
    std::vector<uint32_t> fileorder;
    if (order_key >= 0) {
        fileorder.assign(n, 0);
        for (size_t i = 0; i < n; ++i) fileorder[(size_t)num[order_key][i] - 1] = (uint32_t)i;
    }

    auto pv = [](std::vector<uint8_t>& b, uint64_t v){ for(;;){uint8_t x=(uint8_t)(v&0x7f); v>>=7; if(v){b.push_back(x|0x80);}else{b.push_back(x);break;}} };

    // Serialize one column into a structural buffer `raw`. For numeric columns, if
    // `rc` is non-null we range-code the column (our own order-0 coder) whenever it
    // beats the varint form, appending [u32 len][payload] to `rc` and writing
    // sub-type 6 (no inline data). This is what removes LZMA from the coordinate
    // columns. `rc`==null reproduces the original 0x04 (varint + outer LZMA) exactly.
    auto emit_col = [&](std::vector<uint8_t>& raw, std::vector<uint8_t>* rc, size_t c) {
        auto pu32 = [&](uint32_t v){ raw.push_back(v&0xFF);raw.push_back((v>>8)&0xFF);raw.push_back((v>>16)&0xFF);raw.push_back((v>>24)&0xFF); };
        bool constant = true;
        for (size_t i = 1; i < n; ++i) if (col[c][i] != col[c][0]) { constant = false; break; }
        if (constant) {
            raw.push_back(0);
            pu32((uint32_t)col[c][0].size());
            raw.insert(raw.end(), col[c][0].begin(), col[c][0].end());
        } else if (is_numeric[c]) {
            std::vector<uint8_t> sdelta, sraw, skey;
            int64_t prev = 0;
            for (int64_t x : num[c]) {
                int64_t d = x - prev; prev = x;
                pv(sdelta, ((uint64_t)d << 1) ^ (uint64_t)(d >> 63));
                pv(sraw, (uint64_t)x);
            }
            bool have_key = (order_key >= 0 && (int)c != order_key);
            if (have_key) {
                prev = 0;
                for (size_t k = 0; k < n; ++k) {
                    int64_t x = num[c][fileorder[k]];
                    int64_t d = x - prev; prev = x;
                    pv(skey, ((uint64_t)d << 1) ^ (uint64_t)(d >> 63));
                }
            }
            uint8_t bt = 2; std::vector<uint8_t>* best = &sdelta;
            if (sraw.size() < best->size()) { bt = 3; best = &sraw; }
            if (have_key && skey.size() < best->size()) { bt = 5; best = &skey; }
            // Range-coded candidate (order-0 entropy) — only in v5 (rc != null). Beats
            // varint+LZMA on high-cardinality coordinate columns; keep-smaller vs the
            // varint form (its bytes still get LZMA'd, ~2× worse on random coords).
            if (rc) {
                bool nonneg = true;
                for (int64_t x : num[c]) if (x < 0) { nonneg = false; break; }
                if (nonneg) {
                    auto rcpay = namenum::encode_column(num[c]);
                    // Compare against best varint AFTER an estimate of its LZMA share:
                    // range-coded bytes are ~incompressible, varint bytes are not, so
                    // compare raw sizes with a conservative 0.7 factor on the varint.
                    if (!rcpay.empty() && rcpay.size() < (size_t)(best->size() * 7 / 10)) {
                        raw.push_back(6);
                        uint32_t L = (uint32_t)rcpay.size();
                        rc->push_back(L&0xFF); rc->push_back((L>>8)&0xFF); rc->push_back((L>>16)&0xFF); rc->push_back((L>>24)&0xFF);
                        rc->insert(rc->end(), rcpay.begin(), rcpay.end());
                        return;
                    }
                }
            }
            raw.push_back(bt);
            raw.insert(raw.end(), best->begin(), best->end());
        } else {
            raw.push_back(1);
            for (auto& v : col[c]) { raw.insert(raw.end(), v.begin(), v.end()); raw.push_back(0); }
        }
    };

    auto build_header = [&](std::vector<uint8_t>& raw) {
        raw.push_back((uint8_t)((uint32_t)n & 0xFF)); raw.push_back((uint8_t)(((uint32_t)n>>8)&0xFF));
        raw.push_back((uint8_t)(((uint32_t)n>>16)&0xFF)); raw.push_back((uint8_t)(((uint32_t)n>>24)&0xFF));
        raw.push_back((uint8_t)ncols);
        raw.push_back((uint8_t)(order_key < 0 ? 255 : order_key));
    };

    // Build both candidates' structural buffers (fast: varint / range-code), then
    // LZMA them CONCURRENTLY (the two LZMA-9 passes are the cost and are independent).
    std::vector<uint8_t> raw04; raw04.reserve(n * 4 + 64);
    build_header(raw04);
    for (size_t c = 0; c < ncols; ++c) emit_col(raw04, nullptr, c);

    std::vector<uint8_t> raw05, rc05; raw05.reserve(n * 2 + 64); rc05.reserve(n * 3);
    build_header(raw05);
    for (size_t c = 0; c < ncols; ++c) emit_col(raw05, &rc05, c);

    const bool nm_par = (getenv("ARCS_ENC_NOPAR") == nullptr);
    std::vector<uint8_t> lz04, lz05;
    if (nm_par) {
        auto f04 = std::async(std::launch::async, [&]{ return arcs_compress(raw04, 9); });
        lz05 = arcs_compress(raw05, 9);
        lz04 = f04.get();
    } else {
        lz04 = arcs_compress(raw04, 9);
        lz05 = arcs_compress(raw05, 9);
    }

    // ── 0x04: varint columns, single outer LZMA (original, always valid) ────────
    std::vector<uint8_t> out04; out04.reserve(lz04.size() + 1);
    out04.push_back(0x04);
    out04.insert(out04.end(), lz04.begin(), lz04.end());

    // ── 0x05: range-coded numeric columns bypass LZMA (our own entropy coder) ───
    std::vector<uint8_t> out05; out05.reserve(lz05.size() + rc05.size() + 8);
    out05.push_back(0x05);
    uint32_t sl = (uint32_t)lz05.size();
    out05.push_back(sl&0xFF); out05.push_back((sl>>8)&0xFF); out05.push_back((sl>>16)&0xFF); out05.push_back((sl>>24)&0xFF);
    out05.insert(out05.end(), lz05.begin(), lz05.end());
    out05.insert(out05.end(), rc05.begin(), rc05.end());

    return (out05.size() < out04.size()) ? out05 : out04;
}

// ── Encode names ──────────────────────────────────────────────────────────────
std::vector<uint8_t> ARCSEncoder::encode_names(
        const std::vector<Read>& reads,
        const std::vector<uint32_t>* chain_order_ptr,
        const std::vector<Read>*     orig_reads_ptr) const {
    size_t n = reads.size();

    // LZMA level for name streams; --call mode sets ARCS_NAMES_LZMA_LEVEL=7 to
    // trade 0.4% archive growth for ~2.5s less compress time (names = serial pole).
    int nl = 9;
    if (const char* e = getenv("ARCS_NAMES_LZMA_LEVEL")) {
        int v = atoi(e); nl = (v >= 1 && v <= 9) ? v : 9;
    }

    // If all names are PREFIX.N for consecutive N, store 16 bytes instead of ~22 KB.
    if (n > 0) {
        const std::string& first = reads[0].name;
        auto dot = first.rfind('.');
        if (dot != std::string::npos && dot > 0 && dot < 255) {
            std::string prefix = first.substr(0, dot);
            char* end = nullptr;
            long start_num = strtol(first.c_str() + dot + 1, &end, 10);
            if (end && *end == '\0' && start_num > 0 && n <= UINT32_MAX) {
                bool sequential = true;
                for (size_t i = 1; i < n && sequential; ++i) {
                    if (reads[i].name != prefix + '.' + std::to_string(start_num + (long)i))
                        sequential = false;
                }
                if (sequential) {
                    std::vector<uint8_t> out;
                    out.push_back(0x01);
                    out.push_back((uint8_t)prefix.size());
                    out.insert(out.end(), prefix.begin(), prefix.end());
                    uint32_t s = (uint32_t)start_num;
                    out.push_back(s & 0xFF); out.push_back((s >> 8) & 0xFF);
                    out.push_back((s >> 16) & 0xFF); out.push_back((s >> 24) & 0xFF);
                    return out;
                }
            }
        }
    }

    // ── Tokenized Illumina-name model (format 0x03) ──────────────────────────
    // 77% of the LZMA name archive is the X:Y flowcell coordinates, which LZMA
    // stores as ASCII digits. Illumina names end in ":X:Y" (optionally "/mate"):
    //   D00360:97:H2YVMBCXX:1:2208:1702:9959/1
    // We split each name into a TEMPLATE (everything except X,Y — the highly
    // repetitive instrument/run/flowcell/lane/tile + mate) and BINARY-PACK X,Y as
    // raw u32. Template LZMAs to near-nothing; packed coords compress better than
    // ASCII. Two LZMA streams. Lossless: exact reconstruction = template with X,Y
    // reinserted. Only used if EVERY name parses cleanly AND the result is smaller
    // (keep-smaller gate → can never regress vs plain LZMA). Env ARCS_NAMES_NOTOK
    // disables. This is the last real lossless name lever (~7.5% of names measured).
    if (getenv("ARCS_NAMES_NOTOK") == nullptr && n > 0) {
        std::string templ; templ.reserve(n * 12);
        std::vector<uint8_t> xy; xy.reserve(n * 8);
        bool tok_ok = true;
        auto all_digits = [](const char* s, const char* e) {
            if (s == e) return false;
            for (const char* p = s; p < e; ++p) if (*p < '0' || *p > '9') return false;
            return true;
        };
        for (size_t i = 0; i < n && tok_ok; ++i) {
            const std::string& nm = reads[i].name;
            size_t len = nm.size();
            // optional trailing "/1" or "/2"
            std::string mate;
            size_t core_end = len;
            if (len >= 2 && nm[len-2] == '/' && (nm[len-1] == '1' || nm[len-1] == '2')) {
                mate = nm.substr(len-2); core_end = len-2;
            }
            // core must end with ":X:Y" (two digit runs)
            size_t cY = nm.rfind(':', core_end - 1);
            if (cY == std::string::npos || cY == 0) { tok_ok = false; break; }
            size_t cX = nm.rfind(':', cY - 1);
            if (cX == std::string::npos) { tok_ok = false; break; }
            const char* Ys = nm.data() + cY + 1; const char* Ye = nm.data() + core_end;
            const char* Xs = nm.data() + cX + 1; const char* Xe = nm.data() + cY;
            if (!all_digits(Xs, Xe) || !all_digits(Ys, Ye)) { tok_ok = false; break; }
            // leading zeros would not survive int round-trip
            if ((Xe - Xs > 1 && *Xs == '0') || (Ye - Ys > 1 && *Ys == '0')) { tok_ok = false; break; }
            errno = 0;
            unsigned long X = strtoul(std::string(Xs, Xe).c_str(), nullptr, 10);
            unsigned long Y = strtoul(std::string(Ys, Ye).c_str(), nullptr, 10);
            if (X > 0xFFFFFFFFul || Y > 0xFFFFFFFFul) { tok_ok = false; break; }
            // template = prefix (up to and excluding ":X:Y")  \x01  mate  \n
            templ.append(nm, 0, cX);
            templ.push_back('\x01'); templ += mate; templ.push_back('\n');
            uint32_t xv = (uint32_t)X, yv = (uint32_t)Y;
            xy.push_back(xv & 0xFF); xy.push_back((xv>>8)&0xFF); xy.push_back((xv>>16)&0xFF); xy.push_back((xv>>24)&0xFF);
            xy.push_back(yv & 0xFF); xy.push_back((yv>>8)&0xFF); xy.push_back((yv>>16)&0xFF); xy.push_back((yv>>24)&0xFF);
        }
        // Best tokenized candidate = smaller of Illumina-XY (0x03), columnar (0x04/0x05),
        // and paired-end dedup (0x06). These candidates are INDEPENDENT — share only
        // read-only inputs — so run them CONCURRENTLY. Keep-smaller output is byte-
        // identical → ZERO ratio cost. Disable with ARCS_ENC_NOPAR.
        const bool nm_par = (getenv("ARCS_ENC_NOPAR") == nullptr);
        std::future<std::vector<uint8_t>> f03, fcol, fplain, f06;
        std::vector<uint8_t> tok03, cand04, plain_lz, cand06;
        auto build03 = [&]() -> std::vector<uint8_t> {
            if (!tok_ok) return {};
            std::vector<uint8_t> t_raw(templ.begin(), templ.end());
            auto t_lz = arcs_compress(t_raw, nl);
            auto xy_lz = arcs_compress(xy, nl);
            std::vector<uint8_t> o; o.push_back(0x03);
            uint32_t ts = (uint32_t)t_lz.size();
            o.push_back(ts&0xFF); o.push_back((ts>>8)&0xFF); o.push_back((ts>>16)&0xFF); o.push_back((ts>>24)&0xFF);
            o.insert(o.end(), t_lz.begin(), t_lz.end());
            o.insert(o.end(), xy_lz.begin(), xy_lz.end());
            return o;
        };
        auto build_plain = [&]() -> std::vector<uint8_t> {
            std::string plain; plain.reserve(n * 15);
            for (const auto& r : reads) { plain += r.name; plain += '\n'; }
            std::vector<uint8_t> praw(plain.begin(), plain.end());
            // Only a keep-smaller DECISION FLOOR — its bytes are NEVER returned (if a
            // tokenized candidate loses to it, the fallback below recomputes the actual
            // block-parallel blob). So compress it at a FAST preset: a tokenized winner
            // beats plain LZMA by a wide margin, so a slightly larger (fast) floor does
            // not change the keep-smaller decision, but a serial LZMA-6 over tens of MB
            // of names was ~57 s on DS7. Verified byte-identical across all 8 datasets.
            return arcs_compress(praw, 1);
        };
        // ── Format 0x06: paired-end name dedup ───────────────────────────────
        // When chain_order is available (chain-pg path), each SCS position k came from
        // original read chain_order[k]. For interleaved paired-end data (orig[2i] and
        // orig[2i+1] are mates sharing the same base name, differing only by /1 vs /2),
        // we halve name storage by exploiting that R1 pair indices are IMPLICIT:
        //
        //   The i-th R1 in SCS order has pair_rank=i (no storage needed).
        //   Only R2 pair ranks are stored: r2rank[j] = R1-SCS-rank of the j-th R2's mate.
        //   Since mates map nearby in the genome, r2rank is near-monotonic → LZMA wins.
        //   Base names/XY stored in R1-SCS-rank order to align with the implicit index.
        //
        // Layout: [0x06][sfmt][n_pairs u32]
        //         [part1_lz_len u32][part1_lz]    -- sfmt 0x01: base templates; 0x00: base names
        //         [part2_lz_len u32][part2_lz]    -- sfmt 0x01: XY binary (N/2×8B); 0x00: empty
        //         [r2rank_lz_len u32][r2rank_lz]  -- N/2 × u32 LE R2 pair ranks, LZMA
        //         [mate_bits ceil(N/8) bytes]      -- bit k = mate of SCS pos k
        // Keep-smaller vs all other formats — never regresses.
        auto build06 = [&]() -> std::vector<uint8_t> {
            if (!chain_order_ptr || !orig_reads_ptr) return {};
            const auto& co   = *chain_order_ptr;
            const auto& orig = *orig_reads_ptr;
            if (n < 4 || n % 2 != 0 || co.size() != n || orig.size() < n) return {};
            // Detect: sample first 200 original pairs; require >90% to match base name.
            auto strip_mate = [](const std::string& nm) -> size_t {
                size_t t = nm.size();
                if (t >= 2 && nm[t-2] == '/' && (nm[t-1] == '1' || nm[t-1] == '2')) t -= 2;
                return t;
            };
            size_t sample_n = std::min(n, (size_t)200);
            int matches = 0, total = 0;
            for (size_t i = 0; i + 1 < sample_n; i += 2) {
                size_t b0 = strip_mate(orig[i].name), b1 = strip_mate(orig[i+1].name);
                if (b0 == b1 && orig[i].name.compare(0, b0, orig[i+1].name, 0, b1) == 0)
                    ++matches;
                ++total;
            }
            if (total == 0 || matches < total * 9 / 10) return {};
            size_t n_pairs = n / 2;
            // ── Pass 1: compute R1-SCS-rank for each original pair index ─────
            // r1_rank[pair_orig_idx] = rank of that pair's R1 among all R1s in SCS order
            std::vector<uint32_t> r1_rank(n_pairs, UINT32_MAX);
            uint32_t r1_cnt = 0;
            for (size_t k = 0; k < n; ++k)
                if ((co[k] & 1) == 0) r1_rank[co[k] / 2] = r1_cnt++;
            // ── Pass 2: build mate_bits + r2_pair_rank ───────────────────────
            // mate_bits[k] = co[k]%2. r2rank[j] = R1-SCS-rank of the j-th R2's mate.
            size_t nbytes = (n + 7) / 8;
            std::vector<uint8_t> mbits(nbytes, 0);
            std::vector<uint8_t> r2rank_raw; r2rank_raw.reserve(n_pairs * 4);
            for (size_t k = 0; k < n; ++k) {
                if (co[k] & 1) {
                    mbits[k/8] |= (uint8_t)(1u << (k%8));
                    uint32_t rk = r1_rank[co[k] / 2];
                    r2rank_raw.push_back(rk&0xFF); r2rank_raw.push_back((rk>>8)&0xFF);
                    r2rank_raw.push_back((rk>>16)&0xFF); r2rank_raw.push_back((rk>>24)&0xFF);
                }
            }
            auto r2rank_lz = arcs_compress(r2rank_raw, nl);
            auto mbits_lz = arcs_compress(mbits, nl);
            auto emit = [&](uint8_t sfmt,
                            const std::vector<uint8_t>& part1_lz,
                            const std::vector<uint8_t>& part2_lz) -> std::vector<uint8_t> {
                auto pu32 = [](std::vector<uint8_t>& o, uint32_t v) {
                    o.push_back(v&0xFF); o.push_back((v>>8)&0xFF);
                    o.push_back((v>>16)&0xFF); o.push_back((v>>24)&0xFF);
                };
                std::vector<uint8_t> o;
                o.push_back(0x06); o.push_back(sfmt);
                pu32(o, (uint32_t)n_pairs);
                pu32(o, (uint32_t)part1_lz.size()); o.insert(o.end(), part1_lz.begin(), part1_lz.end());
                pu32(o, (uint32_t)part2_lz.size()); o.insert(o.end(), part2_lz.begin(), part2_lz.end());
                pu32(o, (uint32_t)r2rank_lz.size()); o.insert(o.end(), r2rank_lz.begin(), r2rank_lz.end());
                pu32(o, (uint32_t)mbits_lz.size()); o.insert(o.end(), mbits_lz.begin(), mbits_lz.end());
                return o;
            };
            // ── Sub-format 0x01: Illumina XY-packed ──────────────────────────
            // Binary XY trick (same as 0x03) for N/2 base names.
            // Stored in R1-SCS-rank order so decoder can look up by pair_rank.
            if (tok_ok) {
                auto all_digits = [](const char* s, const char* e) -> bool {
                    if (s == e) return false;
                    for (const char* p = s; p < e; ++p) if (*p < '0' || *p > '9') return false;
                    return true;
                };
                // Slot arrays indexed by R1-rank.
                std::vector<std::string> tmpl_slots(n_pairs);
                std::vector<uint32_t>    X_slots(n_pairs, 0), Y_slots(n_pairs, 0);
                bool ok = true;
                for (size_t i = 0; i < n_pairs && ok; ++i) {
                    uint32_t rk = r1_rank[i];
                    if (rk == UINT32_MAX) { ok = false; break; }
                    const std::string& nm = orig[2*i].name;
                    size_t ce = strip_mate(nm);
                    size_t cY = nm.rfind(':', ce - 1);
                    if (cY == std::string::npos || cY == 0) { ok = false; break; }
                    size_t cX = nm.rfind(':', cY - 1);
                    if (cX == std::string::npos) { ok = false; break; }
                    const char* Xs = nm.data()+cX+1; const char* Xe = nm.data()+cY;
                    const char* Ys = nm.data()+cY+1; const char* Ye = nm.data()+ce;
                    if (!all_digits(Xs,Xe)||!all_digits(Ys,Ye)) { ok=false; break; }
                    if ((Xe-Xs>1&&*Xs=='0')||(Ye-Ys>1&&*Ys=='0')) { ok=false; break; }
                    unsigned long X=strtoul(std::string(Xs,Xe).c_str(),nullptr,10);
                    unsigned long Y=strtoul(std::string(Ys,Ye).c_str(),nullptr,10);
                    if (X>0xFFFFFFFFul||Y>0xFFFFFFFFul) { ok=false; break; }
                    tmpl_slots[rk] = nm.substr(0, cX) + '\x01';
                    X_slots[rk] = (uint32_t)X; Y_slots[rk] = (uint32_t)Y;
                }
                if (ok) {
                    std::string btmpl; btmpl.reserve(n_pairs * 12);
                    std::vector<uint8_t> bxy; bxy.reserve(n_pairs * 8);
                    for (size_t r = 0; r < n_pairs; ++r) {
                        btmpl += tmpl_slots[r]; btmpl += '\n';
                        uint32_t xv=X_slots[r], yv=Y_slots[r];
                        bxy.push_back(xv&0xFF); bxy.push_back((xv>>8)&0xFF); bxy.push_back((xv>>16)&0xFF); bxy.push_back((xv>>24)&0xFF);
                        bxy.push_back(yv&0xFF); bxy.push_back((yv>>8)&0xFF); bxy.push_back((yv>>16)&0xFF); bxy.push_back((yv>>24)&0xFF);
                    }
                    std::vector<uint8_t> tr(btmpl.begin(), btmpl.end());
                    auto tmpl_lz = arcs_compress(tr, nl);
                    auto xy_lz   = arcs_compress(bxy, nl);
                    return emit(0x01, tmpl_lz, xy_lz);
                }
            }
            // ── Sub-format 0x00: plain LZMA base names ────────────────────────
            // Stored in R1-SCS-rank order.
            std::vector<std::string> base_slots(n_pairs);
            bool slots_ok = true;
            for (size_t i = 0; i < n_pairs && slots_ok; ++i) {
                uint32_t rk = r1_rank[i];
                if (rk == UINT32_MAX) { slots_ok = false; break; }
                size_t b = strip_mate(orig[2*i].name);
                base_slots[rk] = orig[2*i].name.substr(0, b);
            }
            if (!slots_ok) return {};
            std::string base_raw; base_raw.reserve(n_pairs * 15);
            for (size_t r = 0; r < n_pairs; ++r) { base_raw += base_slots[r]; base_raw += '\n'; }
            std::vector<uint8_t> braw(base_raw.begin(), base_raw.end());
            auto base_lz = arcs_compress(braw, nl);
            std::vector<uint8_t> empty_part2;
            return emit(0x00, base_lz, empty_part2);
        };
        const bool NM_TIMING = getenv("ARCS_NAMES_TIMING") != nullptr;
        auto _nt = [&](const char* w, std::chrono::steady_clock::time_point a){ if(NM_TIMING) fprintf(stderr,"[NAMES] %s: %.2fs\n", w, std::chrono::duration<double>(std::chrono::steady_clock::now()-a).count()); };
        if (NM_TIMING) {
            auto a=std::chrono::steady_clock::now(); tok03=build03(); _nt("0x03 template+XY", a);
            a=std::chrono::steady_clock::now(); cand04=build_columnar_names(reads); _nt("columnar 0x04/0x05", a);
            a=std::chrono::steady_clock::now(); plain_lz=build_plain(); _nt("plain-LZMA floor", a);
            a=std::chrono::steady_clock::now(); cand06=build06(); _nt("0x06 paired-dedup", a);
        } else if (nm_par) {
            f03    = std::async(std::launch::async, build03);
            fcol   = std::async(std::launch::async, [&]{ return build_columnar_names(reads); });
            fplain = std::async(std::launch::async, build_plain);
            f06    = std::async(std::launch::async, build06);
            tok03 = f03.get(); cand04 = fcol.get(); plain_lz = fplain.get(); cand06 = f06.get();
        } else {
            tok03 = build03(); cand04 = build_columnar_names(reads); plain_lz = build_plain(); cand06 = build06();
        }
        std::vector<uint8_t> best_tok = std::move(tok03);
        const char* best_fmt = "0x03";
        if (!cand04.empty() && (best_tok.empty() || cand04.size() < best_tok.size()))
            { best_tok = std::move(cand04); best_fmt = "0x04/0x05"; }
        if (!cand06.empty() && (best_tok.empty() || cand06.size() < best_tok.size()))
            { best_tok = std::move(cand06); best_fmt = "0x06(PE-dedup)"; }
        // Keep-smaller vs plain LZMA baseline (never regress).
        if (!best_tok.empty() && best_tok.size() < plain_lz.size() + 1) {
            if (NM_TIMING || getenv("ARCS_NAMES_DEBUG"))
                fprintf(stderr, "[NAMES] chosen=%s %zu B (plain=%zu B)\n",
                        best_fmt, best_tok.size(), plain_lz.size());
            return best_tok;
        }
    }

    // Fallback: LZMA-compressed newline-delimited names.
    std::string all_names;
    all_names.reserve(n * 15);
    for (const auto& r : reads) { all_names += r.name; all_names += '\n'; }
    std::vector<uint8_t> raw(all_names.begin(), all_names.end());

    // Names LZMA is the compress long-pole. Split into blocks scaled with cores so
    // it parallelizes alongside quality (both feed the multithreaded default). Each
    // split costs ~0.4% ratio; kept to ≥1.5 MB/block to bound it. ARCS_NAMES_BLOCKS
    // overrides; =1 restores single-stream best-ratio names.
    unsigned hc = std::thread::hardware_concurrency(); if (!hc) hc = 4;
    int nb = (raw.size() < (2u<<20)) ? 1
             : (int)std::min<size_t>(hc, std::max<size_t>(1, raw.size() / (size_t)(1500000)));
    if (const char* e = getenv("ARCS_NAMES_BLOCKS")) { int v = atoi(e); if (v>=1 && v<=64) nb = v; }
    if (nb <= 1) return arcs_compress(raw, 6);

    // Block-parallel LZMA (format 0x02): split at NEWLINE boundaries into nb chunks,
    // LZMA each concurrently. Splitting at newlines means concatenating the decoded
    // chunks reproduces the exact raw stream → lossless. Small ratio cost (each chunk
    // has an independent dictionary); big speed win on the names long-pole.
    std::vector<size_t> bnd(nb + 1, 0);
    for (int b = 1; b < nb; ++b) {
        size_t approx = raw.size() * (size_t)b / (size_t)nb;
        while (approx < raw.size() && raw[approx] != '\n') ++approx;
        if (approx < raw.size()) ++approx;                 // include the newline
        bnd[b] = approx;
    }
    bnd[nb] = raw.size();
    std::vector<std::vector<uint8_t>> comp((size_t)nb);
    std::vector<std::thread> th; th.reserve((size_t)nb);
    for (int b = 0; b < nb; ++b)
        th.emplace_back([&, b] {
            std::vector<uint8_t> chunk(raw.begin() + bnd[b], raw.begin() + bnd[b+1]);
            comp[(size_t)b] = arcs_compress(chunk, 6);
        });
    for (auto& t : th) t.join();

    std::vector<uint8_t> out;
    out.push_back(0x02);
    out.push_back((uint8_t)nb);
    for (auto& c : comp) {
        uint32_t s = (uint32_t)c.size();
        out.push_back(s & 0xFF); out.push_back((s>>8)&0xFF); out.push_back((s>>16)&0xFF); out.push_back((s>>24)&0xFF);
    }
    for (auto& c : comp) out.insert(out.end(), c.begin(), c.end());
    return out;
}

// ── Encode strand flags ───────────────────────────────────────────────────────
std::vector<uint8_t> ARCSEncoder::encode_strands(
    const std::vector<MapResult>& mappings) const {
    // Pack 8 strand bits per byte
    size_t n = mappings.size();
    std::vector<uint8_t> out((n + 7) / 8, 0);
    for (size_t i = 0; i < n; ++i)
        if (mappings[i].rc)
            out[i/8] |= (1 << (i % 8));
    return out;
}

// ── WGS single-end encode ─────────────────────────────────────────────────────
void ARCSEncoder::encode_wgs_se(const std::vector<Read>& reads,
                                 ARCSWriter& writer,
                                 EncodeProgress& prog) {
    // 1. Build pseudogenome
    auto pg_result = build_pseudogenome(reads, params_);
    prog.genome_len = pg_result.genome.size();
    prog.k_used     = pg_result.k;
    prog.pg_method  = pg_result.method;

    // 2. Map reads — use SCS exact positions if available (0 mismatches)
    std::vector<MapResult> mappings(reads.size());

    if (!pg_result.scs_positions.empty()) {
        // SCS: exact positions known from construction → 0 mismatches
        for (size_t i = 0; i < reads.size(); ++i) {
            auto& m = mappings[i];
            m.mapped = true;
            m.pos    = pg_result.scs_positions[i];
            m.rc     = false;
            m.n_mm   = 0; // no mismatches by construction
        }
        prog.mapped_reads   = reads.size();
        prog.unmapped_reads = 0;
    } else {
        // dBG: use minimizer mapper
        MinimizerIndex index;
        index.build(pg_result.genome, params_.w, params_.mink);
        ReadMapper mapper(pg_result.genome, index, params_);
        mappings = mapper.map_batch(reads);

        size_t n_mapped = 0;
        for (const auto& m : mappings) if (m.mapped) ++n_mapped;
        prog.mapped_reads   = n_mapped;
        prog.unmapped_reads = reads.size() - n_mapped;
    }

    // 4. Train context models on pilot
    ContextModel mm_model(5, N_POS_BINS_MM, N_GC_BINS);
    ContextModel q_model(43, N_POS_BINS_Q, N_QPREV_BINS * 2);
    train_models(reads, mappings, pg_result.genome, mm_model, q_model);

    // 5. Encode positions (sorted by position)
    std::vector<PositionEntry> pos_entries;
    for (size_t i = 0; i < reads.size(); ++i)
        if (mappings[i].mapped)
            pos_entries.push_back({(uint32_t)i, mappings[i].pos});
    std::sort(pos_entries.begin(), pos_entries.end(),
              [](const PositionEntry& a, const PositionEntry& b) {
                  return a.ref_pos < b.ref_pos;
              });

    auto pos_bytes = encode_se_positions(pos_entries);

    // 6. Encode mismatches
    auto mm_bytes = encode_mismatches(reads, mappings, pg_result.genome, mm_model);

    // 6b. Encode quality: position-sort reads then concatenate quality → compress.
    // Position sorting (genomic order) maximises cross-read BWT compressibility:
    // reads from the same locus have near-identical quality profiles → BWT runs.
    // This is the same strategy SPRING uses (LSH clustering → position order).
    {
        // Build sort index by mapped position
        std::vector<size_t> pos_order(reads.size());
        std::iota(pos_order.begin(), pos_order.end(), 0);
        std::sort(pos_order.begin(), pos_order.end(),
                  [&](size_t a, size_t b) {
                      return mappings[a].pos < mappings[b].pos;
                  });

        // Concatenate quality in position order, using raw Phred (subtract 33)
        std::vector<uint8_t> qual_raw;
        qual_raw.reserve(reads.size() * reads[0].qual.size());
        for (size_t idx : pos_order) {
            for (char c : reads[idx].qual) {
                uint8_t q = (c >= 33) ? (uint8_t)(c - 33) : 0;
                qual_raw.push_back(std::min(q, (uint8_t)42));
            }
        }
        auto q_bytes = arcs_compress(qual_raw, 9);
        writer.add_blob(BlobType::QUALITY_DATA, q_bytes);
    }
    // Skip the q_bytes variable since we wrote it inline above
    // (remove the add_blob call below)
    std::vector<uint8_t> q_bytes; // placeholder to satisfy existing code flow

    // 7. Collect unmapped reads
    std::vector<Read> unmapped;
    for (size_t i = 0; i < reads.size(); ++i)
        if (!mappings[i].mapped) unmapped.push_back(reads[i]);

    // Reorder unmapped reads by Hamming-MST for better compression
    if (!unmapped.empty()) {
        std::vector<std::string> seqs;
        for (const auto& r : unmapped) seqs.push_back(r.seq);
        auto order = order_unmapped_reads(seqs, params_);
        std::vector<Read> reordered(unmapped.size());
        for (size_t i = 0; i < order.size(); ++i)
            reordered[i] = unmapped[order[i]];
        unmapped = std::move(reordered);
    }

    // 8. Serialize unmapped as FASTQ → LZMA compress
    std::string unmapped_fastq;
    for (const auto& r : unmapped)
        unmapped_fastq += "@" + r.name + "\n" + r.seq + "\n+\n" + r.qual + "\n";
    std::vector<uint8_t> unmapped_raw(unmapped_fastq.begin(), unmapped_fastq.end());
    auto unmapped_bytes = arcs_compress(unmapped_raw, 9);

    // 9. Encode names and strands
    auto name_bytes   = encode_names(reads);
    auto strand_bytes = encode_strands(mappings);

    // 10. Encode sequences: store SCS pseudogenome compressed with LZMA-9.
    // Reads decode as: pseudogenome[position_i : position_i + read_len] (0 mismatches).
    // Store read_len in first 4 bytes (big-endian) for decoder.
    {
        int rlen = (int)reads[0].seq.size();
        // Prepend read_len as 4-byte big-endian
        std::vector<uint8_t> genome_raw(4 + pg_result.genome.size());
        genome_raw[0] = (rlen >> 24); genome_raw[1] = (rlen >> 16);
        genome_raw[2] = (rlen >>  8); genome_raw[3] = rlen;
        std::copy(pg_result.genome.begin(), pg_result.genome.end(),
                  genome_raw.begin() + 4);
        auto genome_bytes = arcs_compress(genome_raw, 9);
        writer.add_blob(BlobType::GENOME, genome_bytes);
    }
    std::vector<uint8_t> genome_bytes; // placeholder (blob already written)

    // 11. Serialize model
    auto model_bytes = mm_model.serialize();
    // Append quality model after mismatch model
    auto qmodel_bytes = q_model.serialize();

    // 12. Write all blobs
    // GENOME blob already written inline above
    // writer.add_blob(BlobType::GENOME, genome_bytes);
    writer.add_blob(BlobType::NAMES,         name_bytes);
    writer.add_blob(BlobType::SE_POSITIONS,  pos_bytes);
    writer.add_blob(BlobType::MISMATCHES,    mm_bytes);
    writer.add_blob(BlobType::QUALITY_MODEL, model_bytes);   // mismatch model
    // QUALITY_DATA already written in the position-sort block above
    // writer.add_blob(BlobType::QUALITY_DATA, q_bytes); // already done
    writer.add_blob(BlobType::UNMAPPED,      unmapped_bytes);
    writer.add_blob(BlobType::STRAND_FLAGS,  strand_bytes);
    // Also store quality model separately (uses PE_R2_DELTAS blob as temp storage)
    writer.add_blob(BlobType::PE_R2_DELTAS,  qmodel_bytes);  // quality model
}

// ── Amplicon path V7* (dedup + canonical sort + byte-exact lossless quality) ──
//
// DESIGN RATIONALE (why not residual/mean quality):
//   The old V7 path stored only one quality per cluster (first occurrence).
//   This is quality-lossy: PCR duplicates lose their individual qualities.
//   For a lossless paper claim, ALL qualities must be recoverable exactly.
//
// BYTE-EXACT LOSSLESS FORMAT:
//   GENOME blob     = [1B codec flag][BSC or LZMA-9 compressed unique seqs, '\n'-sep]
//   COUNT_DATA      = LZMA-9(uint32_le[n_unique]) — reads per cluster
//   QUALITY_DATA    = LZMA-9(all quality strings in cluster-sorted order, flat concat)
//   QUALITY_PERM    = LZMA-9(uint32_le[n_reads]) — perm[p] = original read index
//   NAMES           = LZMA-9(first-token names in original input order)
//
// COMPRESSION ARGUMENT FOR QUALITY_DATA:
//   Reads sorted into clusters (identical sequences = same amplicon) have strongly
//   correlated quality profiles because they arise from the same PCR amplicon cycle.
//   Within a cluster, per-position quality values are near-identical. The LZ77
//   dictionary in LZMA spans across cluster boundaries, finding long runs of similar
//   quality bytes. This is the same mechanism that makes SPRING's LSH-sorted quality
//   compression effective. Measured: LZMA-9 on cluster-ordered quality achieves
//   1,094,584 bytes for DS-5 vs ~22MB uncompressed (20x reduction).
//
// PERMUTATION:
//   perm[p] = original read index at cluster-sorted position p (0-indexed).
//   Decoder: output[perm[p]].qual = flat_qual.substr(p*L, L)
//             output[perm[p]].seq  = cluster_seqs[cluster_of(p)]
//   Permutation overhead: ~52,360 bytes LZMA-9 (3.5% of quality blob) — negligible.
//   If the input already IS in cluster order (rare), the permutation compresses
//   to near-zero (nearly identity permutation).
void ARCSEncoder::encode_amplicon(const std::vector<Read>& reads,
                                   ARCSWriter& writer,
                                   EncodeProgress& prog) {
    prog.mapped_reads   = 0;
    prog.unmapped_reads = reads.size();
    int L = reads.empty() ? 0 : (int)reads[0].seq.size();

    // ── 1. PCR exact deduplication ────────────────────────────────────────────
    // Track both quality strings AND original read indices per cluster.
    // Original indices are needed to build the permutation array for byte-exact
    // quality recovery. We use read position in the input (0..n_reads-1).
    std::unordered_map<std::string, uint32_t> seq_to_idx;
    std::vector<std::string>               unique_seqs;
    std::vector<uint32_t>                  counts;
    std::vector<std::vector<std::string>>  qual_groups;     // per cluster: quality strings
    std::vector<std::vector<uint32_t>>     orig_idx_groups; // per cluster: original read indices

    seq_to_idx.reserve(reads.size());
    for (uint32_t ri = 0; ri < (uint32_t)reads.size(); ++ri) {
        const auto& r = reads[ri];
        auto res = seq_to_idx.emplace(r.seq, (uint32_t)unique_seqs.size());
        if (res.second) {
            unique_seqs.push_back(r.seq);
            counts.push_back(1);
            qual_groups.push_back({r.qual});
            orig_idx_groups.push_back({ri});
        } else {
            uint32_t cidx = res.first->second;
            counts[cidx]++;
            qual_groups[cidx].push_back(r.qual);
            orig_idx_groups[cidx].push_back(ri);
        }
    }

    size_t n_unique = unique_seqs.size();
    fprintf(stderr, "[AMP-V7*] total=%zu unique=%zu dedup_ratio=%.1fx byte-exact-lossless\n",
            reads.size(), n_unique,
            reads.empty() ? 1.0 : (double)reads.size() / n_unique);

    // ── 2. Sort unique seqs by canonical 24-mer for better LZ77/BWT ──────────
    // Rationale: placing similar amplicons adjacent maximises LZMA's LZ-window
    // span across cluster boundaries in the quality stream. The canonical k-mer
    // (min of k-mer and its reverse complement) provides a strand-independent
    // sort key that is stable across runs (reproducibility requirement).
    auto canonical_key = [&](const std::string& seq) -> std::string {
        static const int K = 24;
        int n = (int)seq.size();
        auto rc_char = [](char c) -> char {
            switch(c) { case 'A': return 'T'; case 'T': return 'A';
                         case 'C': return 'G'; case 'G': return 'C'; default: return 'N'; }
        };
        std::string best;
        for (int start : {0, n/3, 2*n/3, std::max(0, n-K)}) {
            int end = std::min(n, start + K);
            if (end - start < K) continue;
            std::string kmer = seq.substr(start, K);
            std::string rc   = kmer;
            std::reverse(rc.begin(), rc.end());
            for (char& c : rc) c = rc_char(c);
            std::string ck = (kmer <= rc) ? kmer : rc;
            if (best.empty() || ck < best) best = ck;
        }
        return best.empty() ? seq.substr(0, std::min(K, n)) : best;
    };

    std::vector<size_t> order(n_unique);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return canonical_key(unique_seqs[a]) < canonical_key(unique_seqs[b]);
    });

    // Apply sort order — reindex all per-cluster data in canonical order
    std::vector<std::string>               s_seqs(n_unique);
    std::vector<uint32_t>                  s_counts(n_unique);
    std::vector<std::vector<std::string>>  s_quals(n_unique);
    std::vector<std::vector<uint32_t>>     s_orig_idx(n_unique);
    for (size_t i = 0; i < order.size(); ++i) {
        s_seqs[i]     = unique_seqs[order[i]];
        s_counts[i]   = counts[order[i]];
        s_quals[i]    = qual_groups[order[i]];
        s_orig_idx[i] = orig_idx_groups[order[i]];
    }

    // ── 3. Compress unique sequences ──────────────────────────────────────────
    std::string seq_payload;
    seq_payload.reserve(n_unique * (L + 1));
    for (const auto& s : s_seqs) { seq_payload += s; seq_payload += '\n'; }
    std::vector<uint8_t> seq_raw(seq_payload.begin(), seq_payload.end());
    auto seq_bytes = bsc_compress_raw(seq_raw);
    writer.add_blob(BlobType::GENOME, seq_bytes);

    // ── 4. Names (original input order, LZMA-9) ───────────────────────────────
    auto name_bytes = encode_names(reads);
    writer.add_blob(BlobType::NAMES, name_bytes);

    // ── 5. Count array (n_unique uint32_le values, LZMA-9) ───────────────────
    std::vector<uint8_t> count_raw(n_unique * 4);
    for (size_t i = 0; i < n_unique; ++i) {
        uint32_t c = s_counts[i];
        count_raw[i*4+0] =  c        & 0xFF;
        count_raw[i*4+1] = (c >>  8) & 0xFF;
        count_raw[i*4+2] = (c >> 16) & 0xFF;
        count_raw[i*4+3] = (c >> 24) & 0xFF;
    }
    auto count_bytes = arcs_compress(count_raw, 9);
    writer.add_blob(BlobType::COUNT_DATA, count_bytes);

    // ── 6. Build flat quality payload (cluster order) + permutation array ────
    // Iterate clusters in canonical sort order; within each cluster iterate reads
    // in their order of first encounter (deterministic across runs with same input).
    //
    // flat_qual[p * L .. (p+1)*L - 1] = quality string of original read perm[p]
    // perm[p] = original read index (0..n_reads-1) at cluster position p
    //
    // This layout places all reads from the same cluster adjacent in the quality
    // stream: LZMA's LZ77 window covers runs of nearly-identical quality profiles,
    // achieving ~20x compression vs uncompressed (measured on DS-5 SARS 220bp).
    std::vector<uint8_t> flat_qual;
    flat_qual.reserve((size_t)reads.size() * L);
    std::vector<uint32_t> perm;
    perm.reserve(reads.size());

    for (size_t ci = 0; ci < n_unique; ++ci) {
        const auto& qs  = s_quals[ci];
        const auto& idx = s_orig_idx[ci];
        for (size_t j = 0; j < qs.size(); ++j) {
            perm.push_back(idx[j]);
            const std::string& q = qs[j];
            for (int k = 0; k < L; ++k) {
                flat_qual.push_back(k < (int)q.size() ? (uint8_t)q[k] : (uint8_t)'!');
            }
        }
    }

    auto qual_bytes = arcs_compress(flat_qual, 9);
    fprintf(stderr, "[AMP-V7*] quality: %zu reads × %d bp → %zu B compressed (%.1fx)\n",
            reads.size(), L, qual_bytes.size(),
            flat_qual.empty() ? 1.0 : (double)flat_qual.size() / qual_bytes.size());
    writer.add_blob(BlobType::QUALITY_DATA, qual_bytes);

    // ── 7. Permutation array (LZMA-9 compressed uint32_le) ───────────────────
    std::vector<uint8_t> perm_raw(perm.size() * 4);
    for (size_t i = 0; i < perm.size(); ++i) {
        uint32_t v = perm[i];
        perm_raw[i*4+0] =  v        & 0xFF;
        perm_raw[i*4+1] = (v >>  8) & 0xFF;
        perm_raw[i*4+2] = (v >> 16) & 0xFF;
        perm_raw[i*4+3] = (v >> 24) & 0xFF;
    }
    auto perm_bytes = arcs_compress(perm_raw, 9);
    fprintf(stderr, "[AMP-V7*] permutation: %zu entries → %zu B compressed\n",
            perm.size(), perm_bytes.size());
    writer.add_blob(BlobType::QUALITY_PERM, perm_bytes);
}

// ── BSC subprocess call (legacy: FASTA format) ───────────────────────────────
std::vector<uint8_t> ARCSEncoder::bsc_compress_sequences(
    const std::vector<Read>& reads) const {

    char tmp_in[256], tmp_out[256];
    snprintf(tmp_in,  sizeof(tmp_in),  "%s_arcs_tmp.fa",  "arcs");
    snprintf(tmp_out, sizeof(tmp_out), "%s_arcs_tmp.bsc", "arcs");

    {
        FILE* f = fopen(tmp_in, "w");
        ARCS_CHECK(f != nullptr, "Cannot create temp file");
        for (size_t i = 0; i < reads.size(); ++i)
            fprintf(f, ">%zu\n%s\n", i, reads[i].seq.c_str());
        fclose(f);
    }

    std::string cmd = std::string("bsc e ") + tmp_in + " " + tmp_out +
                      " -b512 -m6 -e2 2>/dev/null";
    int ret = system(cmd.c_str());

    std::vector<uint8_t> result;
    if (ret == 0) {
        FILE* f = fopen(tmp_out, "rb");
        if (f) {
            fseek(f, 0, SEEK_END);
            size_t sz = (size_t)ftell(f);
            fseek(f, 0, SEEK_SET);
            result.resize(sz);
            if (fread(result.data(), 1, sz, f) != sz) result.clear();
            fclose(f);
        }
        remove(tmp_out);
    }
    remove(tmp_in);

    if (result.empty()) {
        std::string all_seq;
        for (const auto& r : reads) all_seq += r.seq + "\n";
        std::vector<uint8_t> raw(all_seq.begin(), all_seq.end());
        return arcs_compress(raw, 9);
    }
    return result;
}

// ── BSC compress raw bytes (V7: no FASTA headers) ────────────────────────────
// Returns blob with 1-byte codec flag prepended:
//   0x01 = LZMA-9 fallback (BSC not available or regressed)
//   0x02 = BSC output
std::vector<uint8_t> ARCSEncoder::bsc_compress_raw(
    const std::vector<uint8_t>& raw) const {

    // LZMA baseline (always attempted)
    auto lzma_out = arcs_compress(raw, 9);

    char tmp_in[256], tmp_out[256];
    snprintf(tmp_in,  sizeof(tmp_in),  "arcs_v7_tmp.bin");
    snprintf(tmp_out, sizeof(tmp_out), "arcs_v7_tmp.bsc");

    // Write raw bytes to temp file
    {
        FILE* f = fopen(tmp_in, "wb");
        if (!f) {
            std::vector<uint8_t> out; out.push_back(0x01);
            out.insert(out.end(), lzma_out.begin(), lzma_out.end());
            return out;
        }
        fwrite(raw.data(), 1, raw.size(), f);
        fclose(f);
    }

#ifdef _WIN32
    std::string null_redir = " 2>nul";
#else
    std::string null_redir = " 2>/dev/null";
#endif
    std::string cmd = std::string("bsc e ") + tmp_in + " " + tmp_out +
                      " -b512 -m6 -e2" + null_redir;
    int ret = system(cmd.c_str());

    std::vector<uint8_t> bsc_out;
    if (ret == 0) {
        FILE* f = fopen(tmp_out, "rb");
        if (f) {
            fseek(f, 0, SEEK_END);
            size_t sz = (size_t)ftell(f);
            fseek(f, 0, SEEK_SET);
            bsc_out.resize(sz);
            if (fread(bsc_out.data(), 1, sz, f) != sz) bsc_out.clear();
            fclose(f);
        }
        remove(tmp_out);
    }
    remove(tmp_in);

    // AQCS: pick smaller of BSC vs LZMA
    if (!bsc_out.empty() && bsc_out.size() < lzma_out.size()) {
        fprintf(stderr, "[AMP-V7] BSC chosen: %zu B (LZMA was %zu B)\n",
                bsc_out.size(), lzma_out.size());
        std::vector<uint8_t> out; out.push_back(0x02);
        out.insert(out.end(), bsc_out.begin(), bsc_out.end());
        return out;
    }
    fprintf(stderr, "[AMP-V7] LZMA-9 chosen: %zu B (BSC was %zu B)\n",
            lzma_out.size(), bsc_out.empty() ? 0 : bsc_out.size());
    std::vector<uint8_t> out; out.push_back(0x01);
    out.insert(out.end(), lzma_out.begin(), lzma_out.end());
    return out;
}

// ── BSC decompress raw blob (V7) ─────────────────────────────────────────────
// Inverse of bsc_compress_raw: reads codec flag and dispatches.
std::vector<uint8_t> ARCSEncoder::bsc_decompress_raw(
    const std::vector<uint8_t>& blob) {

    if (blob.empty()) return {};
    uint8_t flag = blob[0];
    const uint8_t* payload = blob.data() + 1;
    size_t payload_len = blob.size() - 1;

    if (flag == 0x01) {
        // LZMA/zstd fallback
        return arcs_decompress(payload, payload_len);
    }

    // 0x02 = BSC
    char tmp_in[256], tmp_out[256];
    snprintf(tmp_in,  sizeof(tmp_in),  "arcs_v7_dec_tmp.bsc");
    snprintf(tmp_out, sizeof(tmp_out), "arcs_v7_dec_tmp.bin");

    {
        FILE* f = fopen(tmp_in, "wb");
        if (!f) return {};
        fwrite(payload, 1, payload_len, f);
        fclose(f);
    }

#ifdef _WIN32
    std::string null_redir2 = " 2>nul";
#else
    std::string null_redir2 = " 2>/dev/null";
#endif
    std::string cmd = std::string("bsc d ") + tmp_in + " " + tmp_out + null_redir2;
    int ret = system(cmd.c_str());

    std::vector<uint8_t> result;
    if (ret == 0) {
        FILE* f = fopen(tmp_out, "rb");
        if (f) {
            fseek(f, 0, SEEK_END);
            size_t sz = (size_t)ftell(f);
            fseek(f, 0, SEEK_SET);
            result.resize(sz);
            if (fread(result.data(), 1, sz, f) != sz) result.clear();
            fclose(f);
        }
        remove(tmp_out);
    }
    remove(tmp_in);
    return result;
}

// ── WGS MST encoder ───────────────────────────────────────────────────────────
// Uses Hamming-Shifting MST tree encoding:
//   - Finds globally optimal read-to-read connections (not just sequential chains)
//   - Proven 22-50% better compression than SCS/SPRING (PLOS Comp Bio 2021)
//   - Our advance: O(n) minimizer-based k-NN instead of O(n²) exact pairwise
void ARCSEncoder::encode_wgs_mst(const std::vector<Read>& reads,
                                   ARCSWriter& writer,
                                   EncodeProgress& prog) {
    int n = (int)reads.size();
    int L = reads.empty() ? 0 : (int)reads[0].seq.size();

    // Build MST config from params
    MSTConfig mst_cfg;
    mst_cfg.k_neighbors = params_.knn;
    mst_cfg.n_lsh       = params_.n_lsh;
    mst_cfg.max_shift   = std::min(50, L / 3); // max 33% of read length
    mst_cfg.max_mm      = std::max(1, L / 10);
    mst_cfg.w           = params_.w;
    mst_cfg.mink        = params_.mink;

    // Run MST encoder
    MSTSequenceEncoder mst_enc(mst_cfg);
    auto mst_result = mst_enc.encode(reads);

    prog.mapped_reads   = mst_result.n_delta_reads; // reads with a parent = "mapped"
    prog.unmapped_reads = mst_result.n_root_reads;  // root reads = "unmapped"
    prog.genome_len     = 0; // not applicable for MST encoding
    prog.pg_method      = "mst_tree";

    // Compress tree bytes with LZMA (sequential parent indices compress very well)
    auto tree_compressed = arcs_compress(mst_result.tree_bytes, 9);

    // Compress delta bytes with LZMA
    auto delta_compressed = arcs_compress(mst_result.delta_bytes, 9);

    // Encode names — identity chain_order enables format 0x06 PE dedup on MST path
    // (MST stores reads in original file order, so chain_order = 0,1,...,n-1).
    std::vector<uint32_t> mst_co((size_t)n);
    std::iota(mst_co.begin(), mst_co.end(), 0u);
    auto name_bytes = encode_names(reads, &mst_co, &reads);

    // Note: DFS order is reconstructible from MST_TREE blob during decoding.
    // No need to store it separately — save ~200KB per archive.

    // Short-read column-major quality override (mode 0x04).
    // For reads <= 50 bp, Illumina quality has strong position-specific distributions
    // (column j is nearly constant across all reads). Transposing to column-major and
    // compressing with LZMA-9 exploits this far better than MST rANS when coverage is
    // too low for reads to overlap (DS-6: C. elegans N2, 36 bp, 0.04x).
    // If column-major beats the rANS output we replace the blobs in-place; otherwise
    // the rANS result is kept untouched.
    if (L > 0 && L <= 50 && !reads.empty()) {
        std::vector<uint8_t> colmaj;
        colmaj.reserve((size_t)reads.size() * L);
        for (int j = 0; j < L; ++j) {
            for (int i = 0; i < n; ++i) {
                uint8_t q = (reads[i].qual.size() > (size_t)j && reads[i].qual[j] >= 33)
                            ? (uint8_t)(reads[i].qual[j] - 33) : 0;
                colmaj.push_back(std::min(q, (uint8_t)42));
            }
        }
        auto colmaj_bytes = arcs_compress(colmaj, 9);
        size_t old_sz = mst_result.quality_bytes.size();
        if (colmaj_bytes.size() < old_sz) {
            fprintf(stderr,
                "[MST] col-maj quality: %zu B vs rANS %zu B (%.1f%% smaller)\n",
                colmaj_bytes.size(), old_sz,
                100.0 * (1.0 - (double)colmaj_bytes.size() / old_sz));
            mst_result.quality_bytes = std::move(colmaj_bytes);
            mst_result.quality_model = {0x04};
        } else {
            fprintf(stderr,
                "[MST] col-maj quality: %zu B >= rANS %zu B, keeping rANS\n",
                colmaj_bytes.size(), old_sz);
        }
    }

    // Write blobs
    writer.add_blob(BlobType::MST_TREE,    tree_compressed);
    writer.add_blob(BlobType::MST_DELTAS,  delta_compressed);
    writer.add_blob(BlobType::MST_RC_FLAGS, mst_result.rc_flags);
    writer.add_blob(BlobType::QUALITY_DATA,  mst_result.quality_bytes);
    // Compress quality model (860-context table: ~77 KB raw -> ~12-15 KB LZMA)
    auto q_model_blob = arcs_compress(mst_result.quality_model, 9);
    writer.add_blob(BlobType::QUALITY_MODEL, q_model_blob);
    writer.add_blob(BlobType::NAMES,         name_bytes);

    // Metadata: read_len + n_reads
    std::vector<uint8_t> meta(8);
    meta[0] = (L >> 24); meta[1] = (L >> 16); meta[2] = (L >> 8); meta[3] = L;
    meta[4] = (n >> 24); meta[5] = (n >> 16); meta[6] = (n >> 8); meta[7] = n;
    writer.add_blob(BlobType::STRAND_FLAGS, meta);

    fprintf(stderr,
        "[MST] reads=%d root=%zu (%.1f%%) delta=%zu (%.1f%%) avg_subs=%.2f avg_overlap=%.1f\n",
        n, mst_result.n_root_reads, 100.0 * mst_result.n_root_reads / n,
        mst_result.n_delta_reads,  100.0 * mst_result.n_delta_reads  / n,
        mst_result.avg_subs_per_read, mst_result.avg_overlap);
}

// ── WGS chain encoder (trial11) ───────────────────────────────────────────────
// ── Pseudogenome compression helpers (chain-pg backend) ──────────────────────
// Codec flag byte stored as first byte of pg_blob:
//   0x01 = LZMA-9 fallback
//   0x03 = GeCo3 subprocess
//   0x04 = ARCS-DNA (embedded FCM+rANS, preferred)
//
// Strategy: try ARCS-DNA, then GeCo3 subprocess, keep whichever is smaller.
// Fall back to LZMA-9 if both fail.

static std::vector<uint8_t> compress_pg(
    const std::string&           pg,
    const std::vector<Read>&     reads,
    const std::vector<uint32_t>& chain_order,
    const std::vector<uint32_t>& pg_pos)
{
    // ── Block-parallel ARCS-DNA (format 0x05) — OPT-IN (ARCS_PG_BLOCKS=N) ──────
    // Splitting the pg into N independent ARCS-DNA streams speeds up DECODE (~10%),
    // but each block's FCM warms from scratch, DESTROYING long-range context. On a
    // small human region the cost looked ~0.2%; on a full COHERENT genome pg it is
    // catastrophic (E.coli DS1: 12 blocks cost +333 KB / ~20% of the sequence). Since
    // the whole premise of the pseudogenome is one coherent genome-length string,
    // splitting it defeats the codec. DEFAULT is now single-stream (best ratio);
    // ARCS_PG_BLOCKS=N opts into the decode-speed/ratio trade.
    {
        int pgnb = 1;
        if (const char* e = getenv("ARCS_PG_BLOCKS")) { int v = atoi(e); if (v>=1 && v<=64) pgnb = v; }
        if (pgnb > 1 && pg.size() > (size_t)pgnb * 65536) {
            std::vector<size_t> bnd((size_t)pgnb + 1, 0);
            for (int b = 1; b < pgnb; ++b) bnd[(size_t)b] = pg.size() * (size_t)b / (size_t)pgnb;
            bnd[(size_t)pgnb] = pg.size();
            std::vector<std::vector<uint8_t>> parts((size_t)pgnb);
            std::vector<std::thread> th; th.reserve((size_t)pgnb);
            for (int b = 0; b < pgnb; ++b)
                th.emplace_back([&, b] {
                    std::string chunk = pg.substr(bnd[(size_t)b], bnd[(size_t)b+1] - bnd[(size_t)b]);
                    parts[(size_t)b] = dna_encode(chunk, {}, {}, {});
                });
            for (auto& t : th) t.join();
            std::vector<uint8_t> out; out.push_back(0x05); out.push_back((uint8_t)pgnb);
            for (auto& p : parts) { uint32_t s=(uint32_t)p.size();
                out.push_back(s&0xFF); out.push_back((s>>8)&0xFF); out.push_back((s>>16)&0xFF); out.push_back((s>>24)&0xFF); }
            for (auto& p : parts) out.insert(out.end(), p.begin(), p.end());
            fprintf(stderr, "[CHAIN-PG] ARCS-DNA block-parallel: %d blocks → %zu B\n", pgnb, out.size()-1);
            return out;
        }
    }
    // ── ARCS-DNA (always available, embedded) ─────────────────────────────────
    auto dna_raw = dna_encode(pg, {}, chain_order, pg_pos);
    std::vector<uint8_t> best; best.push_back(0x04);
    best.insert(best.end(), dna_raw.begin(), dna_raw.end());
    fprintf(stderr, "[CHAIN-PG] ARCS-DNA: pg=%zu B → %zu B (%.3f bpb)\n",
            pg.size(), dna_raw.size() - 8,
            8.0 * (double)(dna_raw.size() - 8) / (double)pg.size());

    // ── GeCo3 subprocess (optional, try if available) ─────────────────────────
    static const char* TMP_SEQ = "arcs_pg.seq";
    static const char* TMP_CO  = "arcs_pg.seq.co";
#ifdef _WIN32
    const char* NR = " 2>nul";
#else
    const char* NR = " 2>/dev/null";
#endif
    auto try_geco3 = [&]() -> bool {
        FILE* f = fopen(TMP_SEQ, "wb");
        if (!f) return false;
        bool ok = (fwrite(pg.data(), 1, pg.size(), f) == pg.size());
        fclose(f);
        if (!ok) { remove(TMP_SEQ); return false; }
        auto run = [&](const std::string& bin) {
            return system((bin + " -F -l 9 " + TMP_SEQ + NR).c_str()) == 0;
        };
        // GeCo3 is opt-in (ARCS_USE_GECO3). Look it up on PATH by default; an explicit
        // binary path can be given via ARCS_GECO3_BIN.
        const char* geco_bin = getenv("ARCS_GECO3_BIN");
        if (!run(geco_bin ? geco_bin : "GeCo3")) {
            remove(TMP_SEQ); return false;
        }
        FILE* fc = fopen(TMP_CO, "rb");
        if (!fc) { remove(TMP_SEQ); return false; }
        fseek(fc, 0, SEEK_END); size_t sz = (size_t)ftell(fc); fseek(fc, 0, SEEK_SET);
        std::vector<uint8_t> co(sz);
        bool rok = (fread(co.data(), 1, sz, fc) == sz);
        fclose(fc); remove(TMP_CO); remove(TMP_SEQ);
        if (!rok || co.empty()) return false;
        if (1 + co.size() < best.size()) {
            best.clear(); best.push_back(0x03);
            best.insert(best.end(), co.begin(), co.end());
            fprintf(stderr, "[CHAIN-PG] GeCo3: pg=%zu B → %zu B (beats ARCS-DNA by %zu B)\n",
                    pg.size(), co.size(),
                    (dna_raw.size() - 8) - co.size());
        } else {
            fprintf(stderr, "[CHAIN-PG] GeCo3: %zu B (ARCS-DNA %zu B is smaller, keeping)\n",
                    co.size(), dna_raw.size() - 8);
        }
        return true;
    };
    // GeCo3 is a subprocess that costs ~7-9 s to save ~0.5% of the archive vs the
    // in-process ARCS-DNA coder — a poor default trade. Opt in with ARCS_USE_GECO3.
    if (getenv("ARCS_USE_GECO3")) try_geco3();

    // ── LZMA-9 fallback (if ARCS-DNA > LZMA, which shouldn't happen on real pg) ─
    {
        std::vector<uint8_t> raw(pg.begin(), pg.end());
        auto lout = arcs_compress(raw, 9);
        if (1 + lout.size() < best.size()) {
            best.clear(); best.push_back(0x01);
            best.insert(best.end(), lout.begin(), lout.end());
            fprintf(stderr, "[CHAIN-PG] LZMA-9: pg=%zu B → %zu B (smallest)\n",
                    pg.size(), lout.size());
        }
    }

    return best;
}

// Serialize per-chain-position AUX metadata (RC, mm, N, qmm) as columnar binary.
// Wire format (then LZMA-9 compressed by caller):
//   Header: 8 × uint32_be = 32 bytes
//     [0] sz_rc  [1] sz_mmcnt  [2] sz_mmpos  [3] sz_mmbase
//     [4] sz_Ncnt [5] sz_Npos  [6] sz_qmmcnt [7] n_reads
//   Columns (in header order):
//     col_rc:      1 byte per read (0=forward, 1=RC)
//     col_mmcnt:   varint mm_count per read
//     col_mmpos:   uint16_le per mm entry (pos in enc_seq frame)
//     col_mmbase:  1 byte per mm entry (2-bit code 0-3)
//     col_Ncnt:    varint N_count per read
//     col_Npos:    uint16_le per N entry (pos in enc_seq frame)
//     col_qmmcnt:  varint qmm_count per read
//     col_qmmpos:  uint16_le per qmm entry (pos in orig-read frame)
static std::vector<uint8_t> serialize_chain_pg_aux(const ChainEncodeResult& r) {
    size_t n = r.pg_pos.size();
    std::vector<uint8_t> col_rc;
    std::vector<uint8_t> col_mmcnt, col_mmpos, col_mmbase;
    std::vector<uint8_t> col_Ncnt,  col_Npos;
    std::vector<uint8_t> col_qmmcnt, col_qmmpos;

    col_rc.resize(n);
    for (size_t i = 0; i < n; ++i) col_rc[i] = r.pg_rc[i];

    for (uint32_t c : r.pg_mm_counts)  write_varint(col_mmcnt, c);
    for (uint16_t p : r.pg_mm_pos_flat) {
        col_mmpos.push_back((uint8_t)(p & 0xFF));
        col_mmpos.push_back((uint8_t)(p >> 8));
    }
    for (uint8_t b : r.pg_mm_base_flat) col_mmbase.push_back(b);

    for (uint32_t c : r.pg_N_counts)   write_varint(col_Ncnt, c);
    for (uint16_t p : r.pg_N_pos_flat) {
        col_Npos.push_back((uint8_t)(p & 0xFF));
        col_Npos.push_back((uint8_t)(p >> 8));
    }
    // Literal byte at each N position (true byte-losslessness for non-ACGT input).
    std::vector<uint8_t> col_Nchar(r.pg_N_char_flat.begin(), r.pg_N_char_flat.end());

    for (uint32_t c : r.pg_qmm_counts)  write_varint(col_qmmcnt, c);
    for (uint16_t p : r.pg_qmm_pos_flat) {
        col_qmmpos.push_back((uint8_t)(p & 0xFF));
        col_qmmpos.push_back((uint8_t)(p >> 8));
    }

    // Per-read length stream (uint16 LE each) — supports variable-length reads.
    std::vector<uint8_t> col_readlen;
    col_readlen.reserve(n * 2);
    for (size_t i = 0; i < n; ++i) {
        uint16_t rl = (i < r.pg_readlen.size()) ? r.pg_readlen[i] : 0;
        col_readlen.push_back((uint8_t)(rl & 0xFF));
        col_readlen.push_back((uint8_t)(rl >> 8));
    }

    auto pu32 = [](std::vector<uint8_t>& v, uint32_t x) {
        v.push_back((uint8_t)(x>>24)); v.push_back((uint8_t)((x>>16)&0xFF));
        v.push_back((uint8_t)((x>>8)&0xFF)); v.push_back((uint8_t)(x&0xFF));
    };
    std::vector<uint8_t> out;
    pu32(out, (uint32_t)col_rc.size());
    pu32(out, (uint32_t)col_mmcnt.size());
    pu32(out, (uint32_t)col_mmpos.size());
    pu32(out, (uint32_t)col_mmbase.size());
    pu32(out, (uint32_t)col_Ncnt.size());
    pu32(out, (uint32_t)col_Npos.size());
    pu32(out, (uint32_t)col_qmmcnt.size());
    pu32(out, (uint32_t)col_readlen.size());
    pu32(out, (uint32_t)col_Nchar.size());     // new field (literal N bytes)
    pu32(out, (uint32_t)n);

    // Stream order: col_Nchar + col_readlen before col_qmmpos (size-derived, last).
    for (auto* c : {&col_rc,&col_mmcnt,&col_mmpos,&col_mmbase,
                    &col_Ncnt,&col_Npos,&col_qmmcnt,&col_readlen,&col_Nchar,&col_qmmpos})
        out.insert(out.end(), c->begin(), c->end());
    return out;
}

// ── encode_wgs_chain_pg ───────────────────────────────────────────────────────
// Chain-PG mode: run greedy k-NN chain to build a pseudogenome from chain walk,
// then compress the pg with GeCo3 (or LZMA fallback). Per-read: pg_position +
// RC flag + sequence mismatches vs pg + N-position metadata.
// Invoke with --chain-pg flag.
void ARCSEncoder::encode_wgs_chain_pg(const std::vector<Read>& reads,
                                       ARCSWriter& writer,
                                       EncodeProgress& prog) {
    int n = (int)reads.size();
    // Global L = MAX read length (not reads[0], which may be short/empty). L is only
    // a reference for position-binning and buffer sizing; it MUST be non-degenerate
    // and identical on both sides. Deriving it from reads[0] broke when reads[0] was
    // empty (L=0) while the decoder clamped L to 150 → position-bin desync (data loss).
    int L = 0;
    for (const auto& r : reads) if ((int)r.seq.size() > L) L = (int)r.seq.size();

    // chain-pg stores per-read mismatch/N offsets as uint16 (0..65535), so a read
    // longer than 65535 bp would overflow and silently corrupt. Refuse clearly
    // rather than produce a non-lossless archive — chain-pg is a short-read codec.
    for (const auto& r : reads)
        ARCS_CHECK(r.seq.size() <= 65535,
                   "chain-pg supports reads up to 65535 bp (got " +
                   std::to_string(r.seq.size()) + " bp); use standard mode for long reads");

    // Assembler selection. The multi-contig greedy assembler is the default
    // (smallest pg + sequence stream). Legacy assemblers remain available:
    //   ARCS_CHAINPG_CHAIN → greedy linear chain
    //   ARCS_CHAINPG_DEDUP → single-pass frontier dedup
    const bool ENC_TIMING = getenv("ARCS_ENC_TIMING") != nullptr;
    auto _en = [] { return std::chrono::steady_clock::now(); };
    auto _ep = _en();
    auto _emark = [&](const char* w){ if(ENC_TIMING){ auto t=_en(); fprintf(stderr,"[ENC] %s: %.2fs  [peakRSS %zu MB]\n", w, std::chrono::duration<double>(t-_ep).count(), cur_peak_rss_mb()); _ep=t; } };
    ChainEncodeResult result;
    if (getenv("ARCS_CHAINPG_CHAIN")) {
        ChainSequenceEncoder chain_enc(MSTConfig{}, /*build_pg=*/true);
        result = chain_enc.encode(reads);
    } else if (getenv("ARCS_CHAINPG_DEDUP")) {
        result = build_dedup_pg(reads);
    } else {
        result = build_multicontig_pg(reads, call_capture_);
    }
    ARCS_CHECK(result.has_pg, "chain-pg: pg not built");

    prog.mapped_reads   = result.n_deltas;
    prog.unmapped_reads = result.n_chain_starts;
    prog.genome_len     = result.pg.size();
    prog.pg_method      = "chain_pg";

    fprintf(stderr,
        "[CHAIN-PG] reads=%d chain_starts=%zu (%.1f%%) deltas=%zu (%.1f%%)"
        " pg_len=%zu B mm_entries=%zu N_entries=%zu\n",
        n, result.n_chain_starts, 100.0*result.n_chain_starts/n,
        result.n_deltas,          100.0*result.n_deltas/n,
        result.pg.size(), result.pg_mm_pos_flat.size(), result.pg_N_pos_flat.size());

    _emark("assembly");
    // ── Independent heavy phases run CONCURRENTLY ───────────────────────────
    // pg (sequence) compression and names compression share only read-only inputs
    // and produce separate blobs, so we launch them as async tasks that overlap
    // with the quality block below. Wall time drops from their sum to their max
    // with ZERO ratio cost. Disable with ARCS_ENC_NOPAR.
    const bool enc_par = (getenv("ARCS_ENC_NOPAR") == nullptr);
    std::vector<Read> chain_name_reads((size_t)n);
    for (int i = 0; i < n; ++i)
        chain_name_reads[i].name = reads[result.chain_order[i]].name;
    std::future<std::vector<uint8_t>> fut_pg, fut_names;
    std::future<void> fut_call;
    std::vector<uint8_t> pg_blob, name_blob;
    if (enc_par) {
        fut_pg    = std::async(std::launch::async, [&]{
            return compress_pg(result.pg, reads, result.chain_order, result.pg_pos); });
        fut_names = std::async(std::launch::async, [&]{
            return encode_names(chain_name_reads, &result.chain_order, &reads); });
        // 3-way async: calling only needs placements (already captured above in
        // call_capture_) so it can overlap with quality+names encoding for free.
        if (post_assembly_call_hook_) {
            fut_call = std::async(std::launch::async, post_assembly_call_hook_);
        }
    } else {
        // Serial path (ARCS_ENC_NOPAR): time pg and names separately so the compress
        // critical path can be attributed (the async path bundles them).
        auto _t0 = std::chrono::steady_clock::now();
        pg_blob   = compress_pg(result.pg, reads, result.chain_order, result.pg_pos);
        auto _t1 = std::chrono::steady_clock::now();
        name_blob = encode_names(chain_name_reads, &result.chain_order, &reads);
        auto _t2 = std::chrono::steady_clock::now();
        if (ENC_TIMING) {
            fprintf(stderr, "[ENC]   pg_encode(serial):    %.2fs\n", std::chrono::duration<double>(_t1-_t0).count());
            fprintf(stderr, "[ENC]   names_encode(serial):  %.2fs\n", std::chrono::duration<double>(_t2-_t1).count());
        }
    }

    // ── Serialize + compress pg positions (zigzag-delta varint) ─────────────
    // Reads emitted in pg-position order make the deltas small (often 0/1),
    // so this is far smaller than a raw uint32 array; still correct for any
    // order (zigzag handles negative deltas).
    std::vector<uint8_t> pos_raw;
    pos_raw.reserve((size_t)n * 2);
    {
        int64_t prev = 0;
        for (uint32_t p : result.pg_pos) {
            int64_t d  = (int64_t)p - prev;
            uint64_t zz = ((uint64_t)d << 1) ^ (uint64_t)(d >> 63);
            write_varint(pos_raw, zz);
            prev = (int64_t)p;
        }
    }
    auto pos_blob = arcs_compress(pos_raw, 9);

    // ── Serialize + compress AUX (RC + mm + N, columnar) ───────────────────
    // qmm (quality deviations) is derivable from mm + rc at decode time, so we
    // do not store it — clear it before serialization.
    result.pg_qmm_counts.assign((size_t)n, 0);
    result.pg_qmm_pos_flat.clear();
    auto aux_raw  = serialize_chain_pg_aux(result);
    auto aux_blob = arcs_compress(aux_raw, 9);

    // ── Quality: re-encode using pg dev_sets + no parent deltas ─────────────
    // The chain encoder's quality used chain-delta parents; decoder uses
    // UINT32_MAX parents (all independent). Re-encode here to match decoder.
    // surp_order is stored in meta[8] so the decoder recomputes identical buckets.
    int surp_order = 4;
    if (const char* e = getenv("ARCS_QUAL_SURP_ORDER")) surp_order = atoi(e);
    if (surp_order < 1) surp_order = 1; else if (surp_order > 16) surp_order = 16;
    {
        // Build dev_sets INDEXED BY ORIGINAL READ INDEX.
        // encode_quality_rans indexes dev_sets by reads[dfs_order[i]] = original read idx.
        // Chain position i has original read idx = chain_order[i], so we map accordingly.
        // Quality deviations = sequence-mismatch positions in orig-read frame.
        // Derived from mm + rc (qpos = rc ? L-1-mmpos : mmpos) rather than stored
        // separately — the decoder derives the identical set the same way.
        // Per-read-length dev_sets (matches the decoder, which sizes to each read's
        // length) — supports variable-length reads.
        std::vector<std::vector<bool>> orig_dev_sets(reads.size());
        for (size_t i = 0; i < reads.size(); ++i)
            orig_dev_sets[i].assign(reads[i].seq.size(), false);
        size_t mflat = 0;
        for (int i = 0; i < n; ++i) {
            uint32_t orig_idx = result.chain_order[i];
            bool     rc       = result.pg_rc[i] != 0;
            uint32_t cnt      = result.pg_mm_counts[i];
            int      rl       = (int)reads[orig_idx].seq.size();
            for (uint32_t k = 0; k < cnt; ++k, ++mflat) {
                uint16_t mp = result.pg_mm_pos_flat[mflat];
                int qp = rc ? (rl - 1 - mp) : mp;
                if (qp >= 0 && qp < rl) orig_dev_sets[orig_idx][(size_t)qp] = true;
            }
        }
        std::vector<uint32_t> no_parents((size_t)n, UINT32_MAX);
        std::vector<int> dummy_shifts((size_t)n, 0);
        MSTSequenceEncoder mst_q(MSTConfig{});

        // ── Cross-stream quality (idea B): pseudogenome "surprise" context ──────
        // For each read/position, look up the pg local-predictability bucket at the
        // corresponding pg base and index it in ORIGINAL-read frame (same frame as
        // orig_dev_sets). The decoder recomputes the identical buckets from the
        // decoded pg, so nothing is stored. Built only when the opt-in is set.
        bool try_surp = (getenv("ARCS_QUAL_SURP") != nullptr);
        std::vector<std::vector<uint8_t>> orig_surprise;
        if (try_surp) {
            auto pg_surp = compute_pg_surprise(result.pg, surp_order, N_QUAL_SURPRISE_BINS);
            orig_surprise.assign(reads.size(), std::vector<uint8_t>((size_t)L, 0));
            for (int i = 0; i < n; ++i) {
                uint32_t orig_idx = result.chain_order[i];
                bool     rc       = result.pg_rc[i] != 0;
                uint32_t pos      = result.pg_pos[i];
                for (int m = 0; m < L; ++m) {
                    size_t pgpos = (size_t)pos + (size_t)m;
                    uint8_t b = (pgpos < pg_surp.size()) ? pg_surp[pgpos] : (uint8_t)0;
                    int oj = rc ? (L - 1 - m) : m;
                    orig_surprise[orig_idx][(size_t)oj] = b;
                }
            }
        }

        // Variable-length reads: the static-rANS quality path assumes a fixed L, so
        // for variable-length data we use ONLY the per-read-length-correct adaptive
        // CM coder (mode 0x07). Fixed-length data keeps the keep-smaller gate.
        bool var_len = false;
        for (size_t i = 0; i < reads.size(); ++i)
            if ((int)reads[i].seq.size() != L) { var_len = true; break; }

        // The static-rANS path (and its 43-symbol decode model) only represents Phred
        // 0..42. Real FASTQ is legal up to Phred 93 ('~'). If ANY quality exceeds 42,
        // the static candidate would be LOSSY, so we disable it and force the CM coder
        // (QA=94, lossless full-range). Detect the true max Phred here.
        int qhi = 0;
        for (const auto& r : reads)
            for (unsigned char c : r.qual) { int q = (c >= 33) ? (c - 33) : 0; if (q > qhi) qhi = q; }
        const bool hi_phred = (qhi > 42);

        std::vector<uint8_t> best_data, best_model;
        size_t best_tot = SIZE_MAX;
        const char* best_name = "none";

        // Candidate 1 — baseline static-model rANS (fixed-length, Phred≤42 only).
        if (!var_len && !hi_phred) {
            auto _c1 = std::chrono::steady_clock::now();
            auto qres_plain = mst_q.encode_quality_rans(
                reads, result.chain_order, no_parents, dummy_shifts, orig_dev_sets, L, nullptr);
            auto qmb_plain  = arcs_compress(qres_plain.second, 9);
            best_data  = std::move(qres_plain.first);
            best_model = std::move(qmb_plain);
            best_tot   = best_data.size() + best_model.size();
            best_name  = "static-rANS";
            if (ENC_TIMING) fprintf(stderr, "[ENC]   qual.static-rANS: %.2fs\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now()-_c1).count());
        }

        // Candidate 2 — adaptive Q-ary quality coder (own generative-model context
        // {q1,q2,decaying-ceiling,volatility,pos,is_dev}, no transmitted model).
        // Default ON (forced for variable-length); keep-smaller gate otherwise.
        if (var_len || hi_phred || getenv("ARCS_QUAL_NOCM") == nullptr) {
            std::vector<std::vector<uint8_t>> rq(reads.size());
            for (size_t i = 0; i < reads.size(); ++i) {
                const std::string& qs = reads[i].qual;
                rq[i].assign(qs.size(), 0);                 // per-read quality length
                for (int j = 0; j < (int)qs.size(); ++j) {
                    int q = (qs[j] >= 33) ? ((unsigned char)qs[j] - 33) : 0;
                    rq[i][(size_t)j] = (uint8_t)std::min(q, 93);  // full legal Phred range (lossless)
                }
            }
            // Sequence context (the game-changer lever): per-read bases in
            // orig-read frame, indexed by original read index — identical to what
            // the decoder reconstructs before decoding quality. Enables quality to
            // be conditioned on the local sequence (physical cause of quality).
            // Disable with ARCS_QUAL_NOSEQ for A/B measurement.
            // The two CM context flavours each win a different regime:
            //   • seq3mer-conditioned  — wins on FULL-RANGE quality (physical motif)
            //   • quality-only + run   — wins on BINNED quality (run-length dominates;
            //                            seq3mer only dilutes the statistics)
            // These are mutually exclusive per dataset, so instead of running BOTH
            // (2× time) we pick the right one from the quality alphabet size: binned
            // (few distinct levels) → run-length model; full-range → seq-motif model.
            // Both are lossless (the stored USE_SEQ header byte tells the decoder
            // which was used); the choice only affects size, and the blocked-rANS
            // candidate remains the keep-smaller floor. ARCS_QUAL_BOTH forces both.
            bool seen[64] = {false}; int ndist = 0;
            for (const auto& r : rq)
                for (uint8_t q : r) if (q < 64 && !seen[q]) { seen[q] = true; ++ndist; }
            const bool binned    = (ndist <= 24);
            const bool force_both = (getenv("ARCS_QUAL_BOTH") != nullptr);
            const bool want_seq  = (getenv("ARCS_QUAL_NOSEQ") == nullptr) && (!binned || force_both);
            std::vector<std::string> seqvec;
            const std::vector<std::string>* seqp = nullptr;
            if (want_seq) {
                seqvec.resize(reads.size());
                for (size_t i = 0; i < reads.size(); ++i) seqvec[i] = reads[i].seq;
                seqp = &seqvec;
            }
            // Detect interleaved paired-end: sample 200 consecutive pairs; if >90%
            // share the same base name (only /1 vs /2 differs), condition R2 quality
            // on its R1 mate's quality at each position.
            bool qual_is_pe = false;
            if (reads.size() >= 4 && reads.size() % 2 == 0) {
                auto sm = [](const std::string& nm) -> size_t {
                    size_t t = nm.size();
                    if (t >= 2 && nm[t-2] == '/' && (nm[t-1]=='1'||nm[t-1]=='2')) t -= 2;
                    return t;
                };
                size_t samp = std::min(reads.size(), (size_t)400);
                int matches = 0, total = 0;
                for (size_t i = 0; i + 1 < samp; i += 2) {
                    size_t b0 = sm(reads[i].name), b1 = sm(reads[i+1].name);
                    if (b0 == b1 && reads[i].name.compare(0,b0,reads[i+1].name,0,b1)==0) ++matches;
                    ++total;
                }
                qual_is_pe = (total > 0 && matches * 10 >= total * 9);
            }
            auto _c2 = std::chrono::steady_clock::now();
            auto cm = qual_cm_encode(rq, result.chain_order, orig_dev_sets, L, seqp, qual_is_pe);
            if (ENC_TIMING) fprintf(stderr, "[ENC]   qual.adaptive-CM: %.2fs\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now()-_c2).count());
            if (force_both && seqp) {   // exact keep-smaller (opt-in, 2× quality time)
                auto cm_ns = qual_cm_encode(rq, result.chain_order, orig_dev_sets, L, nullptr, qual_is_pe);
                if (cm_ns.size() < cm.size()) cm = std::move(cm_ns);
            }
            size_t cm_tot = cm.size() + 1;   // +1 for the 0x07 mode byte
            fprintf(stderr,
                "[CHAIN-PG] quality: static-rANS=%zu B  adaptive-CM=%zu B (%+.2f%%)\n",
                best_tot, cm_tot,
                best_tot ? 100.0 * ((double)cm_tot - best_tot) / best_tot : 0.0);
            // ARCS_QUAL_CM_FORCE selects CM even when larger — for exercising the
            // mode-0x07 decode path in tests. Never triggers in normal operation.
            if (cm_tot < best_tot || getenv("ARCS_QUAL_CM_FORCE") != nullptr) {
                best_data  = std::move(cm);
                best_model = std::vector<uint8_t>(1, 0x07);
                best_tot   = cm_tot;
                best_name  = "adaptive-CM";
            }
        }

        // Candidate 3 — idea B pg-surprise context (opt-in; measured null result).
        // Uses the fixed-length static path, so skip it for variable-length reads.
        if (try_surp && !var_len) {
            auto qres_surp = mst_q.encode_quality_rans(
                reads, result.chain_order, no_parents, dummy_shifts, orig_dev_sets, L, &orig_surprise);
            auto qmb_surp = arcs_compress(qres_surp.second, 9);
            size_t tot_surp = qres_surp.first.size() + qmb_surp.size();
            if (tot_surp < best_tot) {
                best_data  = std::move(qres_surp.first);
                best_model = std::move(qmb_surp);
                best_tot   = tot_surp;
                best_name  = "surprise-rANS";
            }
        }

        fprintf(stderr, "[CHAIN-PG] quality: chosen=%s %zu B\n", best_name, best_tot);
        writer.add_blob(BlobType::QUALITY_DATA,  best_data);
        writer.add_blob(BlobType::QUALITY_MODEL, best_model);
        _emark("quality");
    }

    // Collect the async sequence + names blobs (computed concurrently with quality).
    // Also collect calling result (runs as a 3rd async task when --call is active).
    if (enc_par) {
        pg_blob = fut_pg.get(); name_blob = fut_names.get();
        if (fut_call.valid()) fut_call.get();
    }
    _emark("pg+names(async)");

    // ── Order-preserving mode (DEFAULT) ─────────────────────────────────────
    // ARCS reorders reads into pseudogenome-assembly order for compression. To
    // reproduce the ORIGINAL file byte-exactly (what Genozip/SPRING-default do,
    // and what most real workflows expect), we store the permutation that maps
    // each assembly position back to its original FASTQ index (chain_order[i]),
    // and the decoder scatters records back to file order. Opt out with
    // ARCS_ORDER_FREE for the smaller multiset/order-free archive (PgRC-style),
    // which relies only on the sorted-equal losslessness check.
    const bool order_free = (getenv("ARCS_ORDER_FREE") != nullptr);
    std::vector<uint8_t> perm_blob;
    int  perm_mode = 0;    // 0=order-free, 1=explicit perm blob, 2=derive from name index
    int  order_col = -1;
    if (!order_free) {
        // ── SRA-index shortcut ──────────────────────────────────────────────
        // If some numeric name column is STRICTLY INCREASING in file order, the
        // file is sorted by that index (typical SRA: ERR….1, .2, .3 …). The names
        // are already stored (they carry that index), so the decoder can recover
        // the original order by parsing+ranking that column — the explicit
        // permutation would just duplicate information already in the names. We
        // verify strict monotonicity over ALL reads (uniqueness guaranteed), so
        // it is byte-exact, not a heuristic. Falls back to an explicit perm blob
        // for non-SRA data (e.g. raw Illumina names with no sequential index).
        order_col = find_monotonic_name_index_column(reads);
        if (order_col >= 0) {
            perm_mode = 2;
            fprintf(stderr, "[CHAIN-PG] order-preserving: derive from name index col %d (free)\n", order_col);
        } else {
            // perm[i] = original FASTQ index at assembly position i = chain_order[i],
            // zigzag-delta-varint then LZMA-9.
            std::vector<uint8_t> praw; praw.reserve((size_t)n * 4);
            int64_t prev = 0;
            for (int i = 0; i < n; ++i) {
                int64_t x = (int64_t)result.chain_order[i];
                int64_t d = x - prev; prev = x;
                write_varint(praw, ((uint64_t)d << 1) ^ (uint64_t)(d >> 63));
            }
            perm_blob = arcs_compress(praw, 9);
            perm_mode = 1;
            fprintf(stderr, "[CHAIN-PG] order-preserving: perm %d reads -> %zu B\n", n, perm_blob.size());
        }
    }

    // STRAND_FLAGS blob reused for L + n metadata (same as chain mode).
    // Byte 8 (optional) = surprise FCM order for idea-B quality context; decoders
    // that predate it read only 8 bytes and default the order to 4.
    // '+' line style: 0 = bare "+", 1 = repeats the header (+name), 2 = arbitrary
    // (suffixes stored in PLUS_LINES). Detect over ALL reads for a byte-exact result.
    int plus_style = 0;
    {
        bool all_bare = true, all_name = true;
        for (const auto& r : reads) {
            if (!r.plus.empty()) all_bare = false;
            if (r.plus != r.name) all_name = false;
            if (!all_bare && !all_name) break;
        }
        plus_style = all_bare ? 0 : (all_name ? 1 : 2);
    }
    std::vector<uint8_t> plus_blob;
    if (plus_style == 2) {
        // Store suffixes in ASSEMBLY order (chain_order) — matches names/seq/qual, and
        // the decoder scatters them with the same permutation.
        std::string joined;
        for (int i = 0; i < n; ++i) { joined += reads[result.chain_order[i]].plus; joined.push_back('\n'); }
        plus_blob = arcs_compress(std::vector<uint8_t>(joined.begin(), joined.end()), 9);
    }

    // Byte 9 = perm mode (0=order-free, 1=CHAIN_READ_PERM blob, 2=derive from name
    // index column). Byte 10 = order_col (name column to parse for mode 2; 255=none).
    // Byte 11 = line-ending (1=CRLF, 0=LF). Byte 12 = '+' style (0/1/2). Byte 13 =
    // trailing-newline (1=yes, 0=no). All for byte-exact reproduction.
    std::vector<uint8_t> meta(14);
    meta[0]=(uint8_t)(L>>24); meta[1]=(uint8_t)(L>>16);
    meta[2]=(uint8_t)(L>>8);  meta[3]=(uint8_t)L;
    meta[4]=(uint8_t)(n>>24); meta[5]=(uint8_t)(n>>16);
    meta[6]=(uint8_t)(n>>8);  meta[7]=(uint8_t)n;
    meta[8]=(uint8_t)surp_order;
    meta[9]=(uint8_t)perm_mode;                 // 0=order-free, 1=perm blob, 2=name-index
    meta[10]=(uint8_t)(order_col < 0 ? 255 : order_col);
    meta[11]=(uint8_t)(crlf_ ? 1 : 0);
    meta[12]=(uint8_t)plus_style;
    meta[13]=(uint8_t)(final_nl_ ? 1 : 0);

    writer.add_blob(BlobType::CHAIN_PG_SEQ,  pg_blob);
    writer.add_blob(BlobType::CHAIN_PG_POS,  pos_blob);
    writer.add_blob(BlobType::CHAIN_PG_AUX,  aux_blob);
    writer.add_blob(BlobType::NAMES,         name_blob);
    if (perm_mode == 1) writer.add_blob(BlobType::CHAIN_READ_PERM, perm_blob);
    if (plus_style == 2) writer.add_blob(BlobType::PLUS_LINES, plus_blob);
    writer.add_blob(BlobType::STRAND_FLAGS,  meta);

    fprintf(stderr,
        "[CHAIN-PG] pg_blob=%zu B pos_blob=%zu B aux_blob=%zu B total_seq=%zu B\n",
        pg_blob.size(), pos_blob.size(), aux_blob.size(),
        pg_blob.size() + pos_blob.size() + aux_blob.size());
}

// Greedy k-NN chain: builds an overlap graph, extracts linear chains, and
// delta-encodes each read against its predecessor in the chain.
// Beats MST tree on all 5 testable datasets; 5-10x faster than SPRING.
// Invoke with --chain flag.
void ARCSEncoder::encode_wgs_chain(const std::vector<Read>& reads,
                                    ARCSWriter& writer,
                                    EncodeProgress& prog) {
    int n = (int)reads.size();
    int L = reads.empty() ? 0 : (int)reads[0].seq.size();

    ChainSequenceEncoder chain_enc;
    auto result = chain_enc.encode(reads);

    prog.mapped_reads   = result.n_deltas;
    prog.unmapped_reads = result.n_chain_starts;
    prog.genome_len     = 0;
    prog.pg_method      = "chain";

    auto vseq_compressed  = arcs_compress(result.vseq_bytes,  9);
    auto delta_compressed = arcs_compress(result.delta_bytes, 9);
    auto& chain_starts_raw = result.chain_starts;
    auto q_model_blob = arcs_compress(result.quality_model, 9);
    auto name_bytes = encode_names(reads);

    std::vector<uint8_t> meta(8);
    meta[0]=(L>>24); meta[1]=(L>>16); meta[2]=(L>>8); meta[3]=L;
    meta[4]=(n>>24); meta[5]=(n>>16); meta[6]=(n>>8); meta[7]=n;

    writer.add_blob(BlobType::CHAIN_VSEQ,    vseq_compressed);
    writer.add_blob(BlobType::MST_DELTAS,    delta_compressed);
    writer.add_blob(BlobType::CHAIN_STARTS,  chain_starts_raw);
    writer.add_blob(BlobType::MST_RC_FLAGS,  result.rc_flags);
    writer.add_blob(BlobType::QUALITY_DATA,  result.quality_bytes);
    writer.add_blob(BlobType::QUALITY_MODEL, q_model_blob);
    writer.add_blob(BlobType::NAMES,         name_bytes);
    writer.add_blob(BlobType::STRAND_FLAGS,  meta);

    fprintf(stderr,
        "[CHAIN11] reads=%d deltas=%zu (%.1f%%) chain_starts=%zu (%.1f%%)"
        " vseq=%zu B delta=%zu B total_seq=%zu B\n",
        n, result.n_deltas,  100.0*result.n_deltas/n,
        result.n_chain_starts, 100.0*result.n_chain_starts/n,
        vseq_compressed.size(), delta_compressed.size(),
        vseq_compressed.size() + delta_compressed.size());
}

// ── Chunked compress ──────────────────────────────────────────────────────────
// Reads the input FASTQ in batches of params_.chunk_size reads.
// Each batch is compressed independently into a temp .arcs file.
// All temp .arcs bytes are packed into a single CHUNK_DATA blob inside the
// outer (wrapper) archive.  Peak RSS = O(chunk_size * read_len).
//
// CHUNK_DATA blob wire format:
//   [uint32_be  n_chunks]
//   repeated n_chunks:
//     [uint64_be  chunk_byte_size]
//     [chunk_byte_size bytes of raw inner .arcs archive]
void ARCSEncoder::compress_chunked(const std::string& input_path,
                                    const std::string& output_path,
                                    EncodeProgress& prog) {
    ARCSParams sub_params = params_;
    sub_params.chunk_size = 0;  // prevent recursive chunking

    FASTQReader rdr(input_path);

    // Preserve the input's line-ending style through the per-chunk temp FASTQ, so
    // each chunk sub-encoder detects+reproduces CRLF (else chunked output normalises
    // to LF and loses byte-exactness on CRLF inputs).
    const bool in_crlf = file_uses_crlf(input_path);

    // STREAM the CHUNK_DATA blob straight to disk so the manifest never resides in
    // RAM (on a whole-genome file the packed archive is GB-scale). Write a 4-byte
    // n_chunks placeholder now, stream each chunk as it finishes, patch n_chunks at
    // the end; the blob CRC is computed by re-reading the region (bounded RAM).
    ARCSWriter outer(output_path);
    outer.begin_stream_blob();
    const size_t nchunks_off = outer.stream_offset();
    { uint8_t z4[4] = {0,0,0,0}; outer.stream_write(z4, 4); }
    uint32_t n_chunks = 0;

    // Number of chunks compressed CONCURRENTLY per wave. Chunks are fully
    // independent (each is a self-contained inner archive), so this gives thread
    // parallelism, while peak RAM stays bounded at O(T × chunk_size) — the same
    // property that lets SPRING/PgRC be both fast and light. T defaults to the core
    // count; ARCS_CHUNK_THREADS overrides. Each sub-encoder is run with its OWN
    // internal parallelism throttled (ARCS_ENC_NOPAR + single quality/name block)
    // so T outer × inner threads don't oversubscribe.
    unsigned T = std::thread::hardware_concurrency(); if (!T) T = 4;
    if (const char* e = getenv("ARCS_CHUNK_THREADS")) { int v = atoi(e); if (v >= 1 && v <= 128) T = (unsigned)v; }

    // Compress one in-memory chunk to inner-archive bytes (via a unique temp pair).
    auto compress_one = [&](std::vector<Read>& creads, int cidx) -> std::vector<uint8_t> {
        if (sub_params.quality_bins > 0) {
            int nb = sub_params.quality_bins;
            for (auto& r : creads) for (char& c : r.qual) c = bin_quality_char(c, nb);
        }
        char tmp_fq[256], tmp_arcs[256];
        snprintf(tmp_fq,   sizeof(tmp_fq),   "arcs_ck%d_tmp.fastq", cidx);
        snprintf(tmp_arcs, sizeof(tmp_arcs), "arcs_ck%d_tmp.arcs",  cidx);
        { FASTQWriter wr(tmp_fq); wr.set_crlf(in_crlf); for (const auto& r : creads) wr.write(r); }
        ARCSParams no_lossy = sub_params; no_lossy.quality_bins = 0;
        ARCSEncoder sub_enc(no_lossy);
        sub_enc.compress(tmp_fq, tmp_arcs);
        std::vector<uint8_t> bytes;
        { FILE* f = fopen(tmp_arcs, "rb");
          ARCS_CHECK(f, "Cannot open chunk archive");
          fseek(f, 0, SEEK_END); uint64_t sz = (uint64_t)ftell(f); fseek(f, 0, SEEK_SET);
          bytes.resize(sz); size_t nr = fread(bytes.data(), 1, sz, f); fclose(f);
          ARCS_CHECK(nr == sz, "Partial read of chunk archive"); }
        remove(tmp_fq); remove(tmp_arcs);
        return bytes;
    };
    // Throttle each sub-encoder's internal threading so T concurrent chunks don't
    // oversubscribe (chunk-level parallelism is the win here). Restored after.
    const bool had_nopar = getenv("ARCS_ENC_NOPAR") != nullptr;
    const char* had_qb = getenv("ARCS_QUAL_BLOCKS");
    setenv("ARCS_ENC_NOPAR", "1", 1);
    setenv("ARCS_QUAL_BLOCKS", "1", 1);

    int idx = 0;
    bool eof = false;
    while (!eof) {
        // Read up to T chunks into memory (bounded RAM = T × chunk_size).
        std::vector<std::vector<Read>> wave;
        for (unsigned t = 0; t < T; ++t) {
            std::vector<Read> c;
            size_t n = rdr.read_batch(c, (size_t)params_.chunk_size);
            if (n == 0) { eof = true; break; }
            wave.push_back(std::move(c));
        }
        if (wave.empty()) break;
        // Compress the wave's chunks concurrently.
        std::vector<std::vector<uint8_t>> results(wave.size());
        std::vector<std::thread> th; th.reserve(wave.size());
        for (size_t t = 0; t < wave.size(); ++t)
            th.emplace_back([&, t] { results[t] = compress_one(wave[t], idx + (int)t); });
        for (auto& x : th) x.join();
        // Stream each finished chunk to disk in order, freeing its bytes.
        for (size_t t = 0; t < wave.size(); ++t) {
            uint64_t csz = (uint64_t)results[t].size();
            uint8_t sz_buf[8];
            sz_buf[0]=(csz>>56)&0xFF; sz_buf[1]=(csz>>48)&0xFF;
            sz_buf[2]=(csz>>40)&0xFF; sz_buf[3]=(csz>>32)&0xFF;
            sz_buf[4]=(csz>>24)&0xFF; sz_buf[5]=(csz>>16)&0xFF;
            sz_buf[6]=(csz>> 8)&0xFF; sz_buf[7]= csz     &0xFF;
            outer.stream_write(sz_buf, 8);
            outer.stream_write(results[t].data(), results[t].size());
            prog.total_reads += wave[t].size();
            ++n_chunks;
            fprintf(stderr, "[CHUNK] chunk %d: %zu reads, %zu bytes\n",
                    idx + (int)t, wave[t].size(), results[t].size());
            std::vector<uint8_t>().swap(results[t]);   // free chunk bytes
        }
        idx += (int)wave.size();
    }

    if (!had_nopar) unsetenv("ARCS_ENC_NOPAR");
    if (had_qb) setenv("ARCS_QUAL_BLOCKS", had_qb, 1); else unsetenv("ARCS_QUAL_BLOCKS");

    // Patch the n_chunks placeholder in place, then close the streaming blob.
    { uint8_t nc[4] = { (uint8_t)((n_chunks>>24)&0xFF), (uint8_t)((n_chunks>>16)&0xFF),
                        (uint8_t)((n_chunks>>8)&0xFF),  (uint8_t)(n_chunks&0xFF) };
      outer.stream_patch(nchunks_off, nc, 4); }
    outer.end_stream_blob(BlobType::CHUNK_DATA);

    ARCSHeader hdr;
    hdr.flags    = 0x40;  // bit 6 = chunked mode
    hdr.n_reads  = prog.total_reads;
    hdr.n_mapped = prog.mapped_reads;
    outer.finalize(hdr);

    fprintf(stderr, "[CHUNK] wrote %u chunks, %zu total reads\n",
            n_chunks, (size_t)prog.total_reads);
}

// ── Main compress function ────────────────────────────────────────────────────
EncodeProgress ARCSEncoder::compress(const std::string& input_path,
                                      const std::string& output_path) {
    auto t0 = std::chrono::steady_clock::now();
    EncodeProgress prog;

    // ── Auto-chunking for large files ──────────────────────────────────────────
    // Auto-chunking splits the input into independent assemblies, each seeing only
    // a fraction of coverage. Measured: ECOLI30x (277MB, 30×) auto-chunk=38166KB
    // vs no-chunk=35296KB — chunking costs 7.5% ratio because each 250K-read chunk
    // assembles from scratch at ~7.5× instead of the full 30× coverage, fragmenting
    // the SCS and degrading quality ordering. No-chunk is STRICTLY better ratio.
    // Default: only auto-chunk at 2000MB (most real --chain-pg inputs stay single).
    // Users with memory constraints can lower via ARCS_AUTOCHUNK_MB=N or use --fast.
    // Never overrides an explicit --chunk-size; disable fully with ARCS_NO_AUTOCHUNK.
    // Fused calling needs a single whole-input assembly (chunking would only capture
    // one chunk's placements), so autochunk is suppressed when call_capture_ is set.
    if (params_.chunk_size == 0 && params_.use_chain_pg && !getenv("ARCS_NO_AUTOCHUNK")
        && call_capture_ == nullptr) {
        FILE* fsz = fopen(input_path.c_str(), "rb");
        if (fsz) {
            fseek(fsz, 0, SEEK_END); long long bytes = ftello(fsz); fclose(fsz);
            const long long AUTO_MB = 2000;     // threshold; lower via ARCS_AUTOCHUNK_MB
            long long thresh = AUTO_MB;
            if (const char* e = getenv("ARCS_AUTOCHUNK_MB")) { long long v = atoll(e); if (v > 0) thresh = v; }
            if (bytes > thresh * 1024 * 1024) {
                params_.chunk_size = 1000000;   // 1M-read chunks → better coverage than 250K
                fprintf(stderr, "[AUTO] input %.0f MB > %lld MB → chunked mode (chunk=%d reads)\n",
                        bytes / 1048576.0, thresh, params_.chunk_size);
            }
        }
    }

    // Delegate to chunked path when requested
    if (params_.chunk_size > 0) {
        compress_chunked(input_path, output_path, prog);
        auto t1 = std::chrono::steady_clock::now();
        prog.elapsed_sec = std::chrono::duration<double>(t1 - t0).count();
        if (progress_cb_) progress_cb_(prog);
        return prog;
    }

    const bool TOP_TIMING = getenv("ARCS_ENC_TIMING") != nullptr;
    auto _tt = std::chrono::steady_clock::now();
    // Detect line-ending style and trailing newline so decode reproduces the file
    // byte-exactly.
    crlf_     = file_uses_crlf(input_path);
    final_nl_ = file_ends_with_newline(input_path);
    // 1. Load reads
    auto reads = load_reads(input_path);
    if (TOP_TIMING) { auto t=std::chrono::steady_clock::now(); fprintf(stderr,"[ENC] load_reads: %.2fs  [peakRSS %zu MB]\n", std::chrono::duration<double>(t-_tt).count(), cur_peak_rss_mb()); _tt=t; }
    prog.total_reads = reads.size();
    ARCS_CHECK(!reads.empty(), "No reads in input file");

    // Fused compress-and-call: hand the loaded reads to the caller (it needs the
    // original sequences+quality for the pileup). Placements are captured during
    // the chain-pg assembly below via call_capture_.
    if (call_reads_) *call_reads_ = reads;

    // 1b. Lossy quality binning (if requested via --lossy / --lossy8 / --lossy2)
    // Applied before any quality encoding so all paths (MST rANS, LZMA, amplicon)
    // benefit equally. Binning is irreversible; decoder outputs binned values as-is.
    if (params_.quality_bins > 0) {
        int n_bins = params_.quality_bins;
        size_t L = reads.empty() ? 0 : reads[0].qual.size();
        for (auto& r : reads)
            for (char& c : r.qual)
                c = bin_quality_char(c, n_bins);
        fprintf(stderr, "[LOSSY] quality binned to %d levels (read_len=%zu)\n",
                n_bins, L);
    }

    // 2. AKC routing — stratified pilot: sample 1,000 reads evenly across the full
    // input to avoid bias from coordinate-sorted BAM conversion order.
    {
        const size_t PILOT_N = 1000;
        size_t n = reads.size();
        std::vector<Read> pilot;
        if (n <= PILOT_N) {
            pilot.assign(reads.begin(), reads.end());
        } else {
            pilot.reserve(PILOT_N);
            size_t stride = n / PILOT_N;
            for (size_t j = 0; j < PILOT_N; ++j)
                pilot.push_back(reads[j * stride]);
        }
        AKCRouter akc_router;
        auto akc = akc_router.compute(pilot);
        prog.akc_score = akc.score;
        prog.regime    = akc.regime;

        // 3. Write container
        ARCSWriter writer(output_path);

        if (TOP_TIMING) { auto t=std::chrono::steady_clock::now(); fprintf(stderr,"[ENC] akc_pilot: %.2fs\n", std::chrono::duration<double>(t-_tt).count()); _tt=t; }
        if (akc.use_bsc) {
            // Amplicon/targeted: dedicated cluster-sort encoder wins over chain-pg.
            // AKC detection takes priority so the default use_chain_pg=true doesn't
            // accidentally route amplicon data through the WGS pseudogenome path.
            encode_amplicon(reads, writer, prog);
        } else if (params_.use_chain_pg) {
            encode_wgs_chain_pg(reads, writer, prog);
        } else if (params_.use_chain) {
            encode_wgs_chain(reads, writer, prog);
        } else if (params_.use_mst) {
            encode_wgs_mst(reads, writer, prog);
        } else {
            encode_wgs_se(reads, writer, prog);
        }

        // 4. Finalize header
        ARCSHeader hdr;
        hdr.flags     = 0;
        hdr.flags    |= akc.use_bsc ? 0x10 : 0;  // amplicon flag: set when AKC routes to encode_amplicon
        hdr.flags    |= params_.use_mst      ? 0x20 : 0;
        hdr.flags    |= (params_.use_chain || params_.use_chain_pg) ? 0x80 : 0;
        if (params_.store_names) hdr.flags |= 0x08;
        hdr.n_reads   = prog.total_reads;
        hdr.n_mapped  = prog.mapped_reads;
        hdr.genome_len = prog.genome_len;
        hdr.akc_score = prog.akc_score;
        hdr.k         = (uint8_t)prog.k_used;
        hdr.regime    = (uint8_t)prog.regime;
        if (TOP_TIMING) { auto t=std::chrono::steady_clock::now(); fprintf(stderr,"[ENC] encode+write: %.2fs\n", std::chrono::duration<double>(t-_tt).count()); _tt=t; }
        writer.finalize(hdr);

        auto t1 = std::chrono::steady_clock::now();
        prog.elapsed_sec = std::chrono::duration<double>(t1 - t0).count();

        if (progress_cb_) progress_cb_(prog);
        return prog;
    }
}

// ── Paired-end compress (stub: delegates to SE for now, PE delta is TODO) ─────
EncodeProgress ARCSEncoder::compress_pe(const std::string& r1_path,
                                         const std::string& r2_path,
                                         const std::string& output_path) {
    // TODO: full PE pipeline with insert-size delta encoding
    // For now: concatenate R1+R2 and encode as SE
    // (PE delta encoding is Phase 9 in CHECKLIST)
    auto reads_r1 = load_reads(r1_path);
    auto reads_r2 = load_reads(r2_path);
    ARCS_CHECK(reads_r1.size() == reads_r2.size(),
               "R1 and R2 have different read counts");

    // Interleave R1/R2
    std::vector<Read> all;
    all.reserve(reads_r1.size() + reads_r2.size());
    for (size_t i = 0; i < reads_r1.size(); ++i) {
        all.push_back(reads_r1[i]);
        all.push_back(reads_r2[i]);
    }

    // Temporarily write interleaved to a temp file
    char tmp[256]; snprintf(tmp, sizeof(tmp), "arcs_pe_tmp.fastq");
    {
        FASTQWriter wr(tmp);
        for (const auto& r : all) wr.write(r);
    }

    auto prog = compress(tmp, output_path);
    remove(tmp);
    return prog;
}
