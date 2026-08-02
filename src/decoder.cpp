#include <thread>
#include <future>
#include <chrono>
#include "decoder.h"
#include "dna_coder.h"
#include "qual_cm.h"
#include "fastq_io.h"
#include "mst_encoder.h"
#include "chain_encoder.h"
#include "name_num_codec.h"
#include "bsc_codec.h"
#include <cstring>
#include <sstream>
#include <algorithm>
#include <numeric>

static std::vector<std::string> decode_names_blob(
        const std::vector<uint8_t>& blob, size_t n_reads) {
    std::vector<std::string> names;
    if (blob.empty()) return names;
    if (blob[0] == 0x01 && blob.size() >= 7) {
        uint8_t plen = blob[1];
        if (2 + plen + 4 <= blob.size()) {
            std::string prefix(blob.begin() + 2, blob.begin() + 2 + plen);
            uint32_t start = (uint32_t)blob[2+plen]
                           | ((uint32_t)blob[3+plen] << 8)
                           | ((uint32_t)blob[4+plen] << 16)
                           | ((uint32_t)blob[5+plen] << 24);
            names.resize(n_reads);
            for (size_t i = 0; i < n_reads; ++i)
                names[i] = prefix + '.' + std::to_string(start + (uint32_t)i);
            return names;
        }
    }
    // Format 0x03: tokenized Illumina names = LZMA template stream + binary-packed
    // X,Y coordinates. Reconstruct each name as prefix + ":X:Y" + mate.
    if (blob[0] == 0x03 && blob.size() >= 5) {
        uint32_t t_len = (uint32_t)blob[1] | ((uint32_t)blob[2]<<8)
                       | ((uint32_t)blob[3]<<16) | ((uint32_t)blob[4]<<24);
        if (5 + (size_t)t_len <= blob.size()) {
            auto traw = arcs_decompress(blob.data()+5, t_len);
            auto xyraw = arcs_decompress(blob.data()+5+t_len, blob.size()-5-t_len);
            std::string tmpl(traw.begin(), traw.end());
            const uint8_t* xp = xyraw.data();
            size_t xn = xyraw.size();
            size_t pos = 0, xi = 0;
            while (pos < tmpl.size()) {
                size_t nl = tmpl.find('\n', pos);
                if (nl == std::string::npos) break;
                std::string line = tmpl.substr(pos, nl - pos);
                pos = nl + 1;
                size_t sep = line.find('\x01');
                std::string prefix = (sep==std::string::npos) ? line : line.substr(0, sep);
                std::string mate   = (sep==std::string::npos) ? std::string() : line.substr(sep+1);
                if (xi + 8 > xn) break;
                uint32_t X = (uint32_t)xp[xi] | ((uint32_t)xp[xi+1]<<8) | ((uint32_t)xp[xi+2]<<16) | ((uint32_t)xp[xi+3]<<24);
                uint32_t Y = (uint32_t)xp[xi+4] | ((uint32_t)xp[xi+5]<<8) | ((uint32_t)xp[xi+6]<<16) | ((uint32_t)xp[xi+7]<<24);
                xi += 8;
                names.push_back(prefix + ':' + std::to_string(X) + ':' + std::to_string(Y) + mate);
            }
            return names;
        }
    }
    // Format 0x04: columnar tokenized names (per-field delta). Reconstruct each
    // column (constant literal / numeric delta-zigzag-varint / NUL-separated literal)
    // then concatenate columns in order to rebuild each name exactly.
    if (blob[0] == 0x04 && blob.size() >= 2) {
        auto raw = arcs_decompress(blob.data() + 1, blob.size() - 1);
        const uint8_t* p = raw.data();
        const uint8_t* e = p + raw.size();
        auto gu32 = [&]() -> uint32_t {
            if (p + 4 > e) return 0;
            uint32_t v = (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24);
            p += 4; return v;
        };
        auto gvar = [&]() -> uint64_t {
            uint64_t v = 0; int s = 0;
            while (p < e) { uint8_t b = *p++; v |= (uint64_t)(b & 0x7f) << s; if (!(b & 0x80)) break; s += 7; }
            return v;
        };
        if (p + 6 > e) return names;
        uint32_t nn = gu32();
        uint8_t ncols = *p++;
        uint8_t order_key = *p++;                            // 255 = none
        // For key-order columns we need the order key decoded first. The key column
        // itself is stored as an ordinary numeric column (type 2/3) earlier or later
        // in the stream; we must decode it before scattering key-order fields. To keep
        // one linear pass, we record raw per-column payloads, decode the key, then
        // materialise. Simpler: decode all columns into their stored order, then fix
        // key-order ones using the decoded key.
        std::vector<std::vector<std::string>> col((size_t)ncols);
        std::vector<int> ctype((size_t)ncols, -1);
        std::vector<std::vector<int64_t>> keyvals((size_t)ncols);   // for type-5 (file order)
        for (int c = 0; c < ncols; ++c) {
            if (p >= e) return names;
            uint8_t type = *p++; ctype[(size_t)c] = type;
            if (type == 0) {                                 // constant literal
                uint32_t len = gu32();
                if (p + len > e) return names;
                std::string s((const char*)p, len); p += len;
                col[(size_t)c].assign(nn, s);
            } else if (type == 2) {                          // numeric delta (chain order)
                col[(size_t)c].resize(nn);
                int64_t prev = 0;
                for (uint32_t i = 0; i < nn; ++i) {
                    uint64_t zz = gvar();
                    int64_t d = (int64_t)(zz >> 1) ^ -(int64_t)(zz & 1);
                    prev += d; col[(size_t)c][i] = std::to_string(prev);
                }
            } else if (type == 3) {                          // numeric raw (chain order)
                col[(size_t)c].resize(nn);
                for (uint32_t i = 0; i < nn; ++i)
                    col[(size_t)c][i] = std::to_string((int64_t)gvar());
            } else if (type == 5) {                          // numeric delta in FILE order
                keyvals[(size_t)c].resize(nn);
                int64_t prev = 0;
                for (uint32_t i = 0; i < nn; ++i) {
                    uint64_t zz = gvar();
                    int64_t d = (int64_t)(zz >> 1) ^ -(int64_t)(zz & 1);
                    prev += d; keyvals[(size_t)c][i] = prev;   // value for file position i
                }
            } else {                                         // NUL-separated literal
                col[(size_t)c].resize(nn);
                for (uint32_t i = 0; i < nn; ++i) {
                    const uint8_t* q = p;
                    while (p < e && *p) ++p;
                    col[(size_t)c][i].assign((const char*)q, p - q);
                    if (p < e) ++p;
                }
            }
        }
        // Scatter key-order (type-5) columns back to record order using the key.
        if (order_key < ncols && !col[(size_t)order_key].empty()) {
            for (int c = 0; c < ncols; ++c) {
                if (ctype[(size_t)c] != 5) continue;
                col[(size_t)c].resize(nn);
                for (uint32_t i = 0; i < nn; ++i) {
                    int64_t k = (int64_t)strtoll(col[(size_t)order_key][i].c_str(), nullptr, 10);  // file pos of record i
                    if (k >= 1 && (uint32_t)k <= nn)
                        col[(size_t)c][i] = std::to_string(keyvals[(size_t)c][(size_t)k - 1]);
                }
            }
        }
        names.resize(nn);
        for (uint32_t i = 0; i < nn; ++i) {
            std::string& nm = names[i];
            for (int c = 0; c < ncols; ++c) nm += col[(size_t)c][i];
        }
        return names;
    }
    // Format 0x05: like 0x04, but numeric columns may be range-coded (sub-type 6,
    // our own order-0 entropy coder — no LZMA). Layout: [0x05][u32 struct_lzma_len]
    // [LZMA of structural buffer][rc_section]; sub-type-6 columns pull their
    // [u32 len][payload] from rc_section in column order.
    if (blob[0] == 0x05 && blob.size() >= 5) {
        uint32_t sl = (uint32_t)blob[1] | ((uint32_t)blob[2]<<8) | ((uint32_t)blob[3]<<16) | ((uint32_t)blob[4]<<24);
        if (5 + (size_t)sl > blob.size()) return names;
        auto raw = arcs_decompress(blob.data() + 5, sl);
        const uint8_t* rcp = blob.data() + 5 + sl;      // range-coded section
        const uint8_t* rce = blob.data() + blob.size();
        const uint8_t* p = raw.data();
        const uint8_t* e = p + raw.size();
        auto gu32 = [&]() -> uint32_t {
            if (p + 4 > e) return 0;
            uint32_t v = (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24);
            p += 4; return v;
        };
        auto gvar = [&]() -> uint64_t {
            uint64_t v = 0; int s = 0;
            while (p < e) { uint8_t b = *p++; v |= (uint64_t)(b & 0x7f) << s; if (!(b & 0x80)) break; s += 7; }
            return v;
        };
        if (p + 6 > e) return names;
        uint32_t nn = gu32();
        uint8_t ncols = *p++;
        uint8_t order_key = *p++;
        std::vector<std::vector<std::string>> col((size_t)ncols);
        std::vector<int> ctype((size_t)ncols, -1);
        std::vector<std::vector<int64_t>> keyvals((size_t)ncols);
        for (int c = 0; c < ncols; ++c) {
            if (p >= e) return names;
            uint8_t type = *p++; ctype[(size_t)c] = type;
            if (type == 0) {
                uint32_t len = gu32();
                if (p + len > e) return names;
                std::string s((const char*)p, len); p += len;
                col[(size_t)c].assign(nn, s);
            } else if (type == 2) {
                col[(size_t)c].resize(nn); int64_t prev = 0;
                for (uint32_t i = 0; i < nn; ++i) { uint64_t zz = gvar(); int64_t d = (int64_t)(zz>>1)^-(int64_t)(zz&1); prev += d; col[(size_t)c][i] = std::to_string(prev); }
            } else if (type == 3) {
                col[(size_t)c].resize(nn);
                for (uint32_t i = 0; i < nn; ++i) col[(size_t)c][i] = std::to_string((int64_t)gvar());
            } else if (type == 5) {
                keyvals[(size_t)c].resize(nn); int64_t prev = 0;
                for (uint32_t i = 0; i < nn; ++i) { uint64_t zz = gvar(); int64_t d = (int64_t)(zz>>1)^-(int64_t)(zz&1); prev += d; keyvals[(size_t)c][i] = prev; }
            } else if (type == 6) {                          // range-coded numeric column
                if (rcp + 4 > rce) return names;
                uint32_t L = (uint32_t)rcp[0] | ((uint32_t)rcp[1]<<8) | ((uint32_t)rcp[2]<<16) | ((uint32_t)rcp[3]<<24);
                rcp += 4;
                if (rcp + L > rce) return names;
                auto vals = namenum::decode_column(rcp, rcp + L);
                rcp += L;
                col[(size_t)c].resize(nn);
                for (uint32_t i = 0; i < nn && i < vals.size(); ++i) col[(size_t)c][i] = std::to_string(vals[i]);
            } else {
                col[(size_t)c].resize(nn);
                for (uint32_t i = 0; i < nn; ++i) {
                    const uint8_t* q = p; while (p < e && *p) ++p;
                    col[(size_t)c][i].assign((const char*)q, p - q);
                    if (p < e) ++p;
                }
            }
        }
        if (order_key < ncols && !col[(size_t)order_key].empty()) {
            for (int c = 0; c < ncols; ++c) {
                if (ctype[(size_t)c] != 5) continue;
                col[(size_t)c].resize(nn);
                for (uint32_t i = 0; i < nn; ++i) {
                    int64_t k = (int64_t)strtoll(col[(size_t)order_key][i].c_str(), nullptr, 10);
                    if (k >= 1 && (uint32_t)k <= nn) col[(size_t)c][i] = std::to_string(keyvals[(size_t)c][(size_t)k - 1]);
                }
            }
        }
        names.resize(nn);
        for (uint32_t i = 0; i < nn; ++i) { std::string& nm = names[i]; for (int c = 0; c < ncols; ++c) nm += col[(size_t)c][i]; }
        return names;
    }
    // Format 0x02: block-parallel LZMA (nb chunks split at newline boundaries).
    if (blob[0] == 0x02 && blob.size() >= 2) {
        int nb = blob[1];
        size_t hdr = 2 + (size_t)nb * 4;
        if (nb >= 1 && nb <= 64 && blob.size() >= hdr) {
            std::vector<uint32_t> clen((size_t)nb);
            size_t total = hdr;
            for (int b = 0; b < nb; ++b) {
                clen[(size_t)b] = (uint32_t)blob[2+b*4] | ((uint32_t)blob[3+b*4]<<8)
                                | ((uint32_t)blob[4+b*4]<<16) | ((uint32_t)blob[5+b*4]<<24);
                total += clen[(size_t)b];
            }
            if (total == blob.size()) {
                std::vector<size_t> off((size_t)nb + 1); off[0] = hdr;
                for (int b = 0; b < nb; ++b) off[(size_t)b+1] = off[(size_t)b] + clen[(size_t)b];
                std::vector<std::string> parts((size_t)nb);
                std::vector<std::thread> th; th.reserve((size_t)nb);
                for (int b = 0; b < nb; ++b)
                    th.emplace_back([&, b] {
                        auto r = arcs_decompress(blob.data() + off[(size_t)b], clen[(size_t)b]);
                        parts[(size_t)b].assign(r.begin(), r.end());
                    });
                for (auto& t : th) t.join();
                std::string s;
                for (auto& p : parts) s += p;
                std::istringstream ss(s);
                std::string line;
                while (std::getline(ss, line)) if (!line.empty()) names.push_back(line);
                return names;
            }
        }
    }
    // Format 0x06: paired-end name dedup (R1-rank implicit scheme).
    // Layout: [0x06][sfmt u8][n_pairs u32][part1_lz_len u32][part1_lz]
    //         [part2_lz_len u32][part2_lz][r2rank_lz_len u32][r2rank_lz]
    //         [mbits_lz_len u32][mbits_lz]
    // sfmt 0x01: part1=base template LZMA (R1-SCS-rank order), part2=XY binary LZMA
    // sfmt 0x00: part1=base names LZMA (R1-SCS-rank order), part2 unused (len=0)
    // r2rank: N/2 uint32_t LE — R1-SCS-rank of each R2's mate (near-monotonic → small LZMA)
    // mate_bits: LZMA of ceil(N/8) bytes; bit k = 1 if SCS pos k is R2.
    if (blob[0] == 0x06 && blob.size() >= 18) {
        const uint8_t* p = blob.data() + 1;
        const uint8_t* e = blob.data() + blob.size();
        auto gu32 = [&]() -> uint32_t {
            if (p + 4 > e) return 0;
            uint32_t v = (uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);
            p += 4; return v;
        };
        uint8_t sfmt = *p++;
        uint32_t n_pairs = gu32();
        uint32_t p1l = gu32(); if (p + p1l > e) return names;
        auto part1 = arcs_decompress(p, p1l); p += p1l;
        uint32_t p2l = gu32(); if (p + p2l > e) return names;
        auto part2 = (p2l > 0) ? arcs_decompress(p, p2l) : std::vector<uint8_t>{};
        p += p2l;
        uint32_t r2rl = gu32(); if (p + r2rl > e) return names;
        auto r2rank_raw = arcs_decompress(p, r2rl); p += r2rl;
        uint32_t mbl = gu32(); if (p + mbl > e) return names;
        auto mbits_raw = arcs_decompress(p, mbl); p += mbl;
        size_t n_total = (size_t)n_pairs * 2;
        size_t nbytes  = (n_total + 7) / 8;
        if (r2rank_raw.size() < n_pairs * 4 || mbits_raw.size() < nbytes) return names;
        const uint8_t* mbits = mbits_raw.data();
        // Decode r2rank into a flat array for fast lookup.
        std::vector<uint32_t> r2rank(n_pairs);
        for (size_t j = 0; j < n_pairs; ++j)
            r2rank[j] = (uint32_t)r2rank_raw[j*4]|((uint32_t)r2rank_raw[j*4+1]<<8)
                      | ((uint32_t)r2rank_raw[j*4+2]<<16)|((uint32_t)r2rank_raw[j*4+3]<<24);
        // Reconstruct N names in SCS order.
        names.resize(n_total);
        if (sfmt == 0x01) {
            // Illumina XY-packed (R1-SCS-rank order).
            std::string tmpl(part1.begin(), part1.end());
            std::vector<std::string> prefixes; prefixes.reserve(n_pairs);
            size_t pos = 0;
            while (prefixes.size() < n_pairs && pos < tmpl.size()) {
                size_t nl = tmpl.find('\n', pos);
                if (nl == std::string::npos) break;
                std::string line = tmpl.substr(pos, nl - pos);
                size_t sep = line.find('\x01');
                prefixes.push_back(sep == std::string::npos ? line : line.substr(0, sep));
                pos = nl + 1;
            }
            const uint8_t* xp = part2.data();
            size_t xn = part2.size();
            uint32_t r1_cnt = 0, r2_cnt = 0;
            for (size_t k = 0; k < n_total; ++k) {
                uint8_t mate = (mbits[k/8] >> (k%8)) & 1u;
                uint32_t pr = (mate == 0) ? r1_cnt++ : r2rank[r2_cnt++];
                if (pr >= (uint32_t)prefixes.size() || pr*8+8 > xn) continue;
                const uint8_t* xb = xp + pr * 8;
                uint32_t X=(uint32_t)xb[0]|((uint32_t)xb[1]<<8)|((uint32_t)xb[2]<<16)|((uint32_t)xb[3]<<24);
                uint32_t Y=(uint32_t)xb[4]|((uint32_t)xb[5]<<8)|((uint32_t)xb[6]<<16)|((uint32_t)xb[7]<<24);
                names[k] = prefixes[pr] + ':' + std::to_string(X) + ':' + std::to_string(Y)
                         + (mate ? "/2" : "/1");
            }
        } else {
            // Sub-format 0x00: plain LZMA base names (R1-SCS-rank order).
            std::vector<std::string> base_names; base_names.reserve(n_pairs);
            size_t pos = 0;
            while (base_names.size() < n_pairs && pos < part1.size()) {
                size_t nl = pos;
                while (nl < part1.size() && part1[nl] != '\n') ++nl;
                base_names.emplace_back((const char*)part1.data() + pos, nl - pos);
                pos = (nl < part1.size()) ? nl + 1 : nl + 1;
            }
            uint32_t r1_cnt = 0, r2_cnt = 0;
            for (size_t k = 0; k < n_total; ++k) {
                uint8_t mate = (mbits[k/8] >> (k%8)) & 1u;
                uint32_t pr = (mate == 0) ? r1_cnt++ : r2rank[r2_cnt++];
                if (pr < (uint32_t)base_names.size())
                    names[k] = base_names[pr] + (mate ? "/2" : "/1");
            }
        }
        return names;
    }
    auto raw = arcs_decompress(blob.data(), blob.size());
    std::string s(raw.begin(), raw.end());
    std::istringstream ss(s);
    std::string line;
    while (std::getline(ss, line)) if (!line.empty()) names.push_back(line);
    return names;
}

// ── reconstruct_seq ───────────────────────────────────────────────────────────
std::string ARCSDecoder::reconstruct_seq(
    const std::string& ref_substr,
    const std::vector<uint32_t>& mm_offsets,
    const std::vector<uint8_t>&  mm_bases,
    bool rc) const {

    std::string seq = ref_substr;

    // Apply mismatches: at each mm_offset, the original read had mm_bases[j]
    // (i.e., the pseudogenome has a different base there)
    for (size_t j = 0; j < mm_offsets.size(); ++j) {
        uint32_t off = mm_offsets[j];
        if (off < seq.size() && mm_bases[j] < 5)
            seq[off] = BASE_TO_CHAR[mm_bases[j]];
    }

    if (rc) seq = reverse_complement(seq);
    return seq;
}

// ── decode_quality ─────────────────────────────────────────────────────────────
std::string ARCSDecoder::decode_quality(
    const uint8_t*& ptr,
    const uint8_t* end,
    RansDecoder& dec,
    const ContextModel& q_model,
    const std::vector<bool>& is_dev,
    int read_len) const {

    std::string qual(read_len, '!');
    uint8_t q_prev = 30;

    for (int j = 0; j < read_len; ++j) {
        ContextKey ctx = quality_ctx(q_prev, is_dev.empty() ? false : is_dev[j], j, read_len);
        uint8_t q = const_cast<ContextModel&>(q_model).decode_sym(dec, ptr, end, ctx);
        qual[j]   = (char)(q + 33); // raw Phred → ASCII Phred+33
        q_prev    = q;
    }
    return qual;
}

// ── decode_wgs ────────────────────────────────────────────────────────────────
void ARCSDecoder::decode_wgs(const ARCSReader& rdr, FASTQWriter& out_writer) {
    const auto& hdr = rdr.header();
    size_t n_reads  = hdr.n_reads;

    // 1. Load genome
    auto genome_blob = rdr.read_blob(BlobType::GENOME);
    ARCS_CHECK(!genome_blob.empty(), "GENOME blob missing");
    auto genome_raw = arcs_decompress(genome_blob.data(), genome_blob.size());
    std::string genome(genome_raw.begin(), genome_raw.end());

    // 2. Load names
    std::vector<std::string> names;
    if (rdr.has_blob(BlobType::NAMES))
        names = decode_names_blob(rdr.read_blob(BlobType::NAMES), n_reads);
    while (names.size() < n_reads)
        names.push_back("read_" + std::to_string(names.size()));

    // 3. Load strand flags
    std::vector<bool> strand(n_reads, false);
    if (rdr.has_blob(BlobType::STRAND_FLAGS)) {
        auto sb = rdr.read_blob(BlobType::STRAND_FLAGS);
        for (size_t i = 0; i < n_reads; ++i)
            if (i/8 < sb.size())
                strand[i] = (sb[i/8] >> (i%8)) & 1;
    }

    // 4. Load positions
    auto pos_blob = rdr.read_blob(BlobType::SE_POSITIONS);
    ARCS_CHECK(!pos_blob.empty(), "SE_POSITIONS blob missing");
    auto pos_entries = decode_se_positions(pos_blob.data(), pos_blob.size(), n_reads);

    // Build lookup: read_idx → position
    std::vector<pos_t>  read_pos(n_reads, 0);
    std::vector<bool>   read_mapped(n_reads, false);
    for (const auto& e : pos_entries) {
        if (e.read_idx < n_reads) {
            read_pos[e.read_idx]    = e.ref_pos;
            read_mapped[e.read_idx] = true;
        }
    }

    // 5. Load and deserialize context models
    ContextModel mm_model(5, N_POS_BINS_MM, N_GC_BINS);
    ContextModel q_model(43, N_POS_BINS_Q, N_QPREV_BINS * 2);

    auto mm_model_blob = rdr.read_blob(BlobType::QUALITY_MODEL);
    ARCS_CHECK(!mm_model_blob.empty(), "QUALITY_MODEL blob missing");
    mm_model.deserialize(mm_model_blob.data(), mm_model_blob.size());

    // Quality model is appended after mismatch model in same blob
    // (encoder appended them; decoder needs separate blobs — TODO: split properly)
    // For now: use mm_model for both (will fix in v2 format)
    // q_model.deserialize(qmodel_bytes.data(), qmodel_bytes.size());

    // 6. Load mismatches stream
    auto mm_blob = rdr.read_blob(BlobType::MISMATCHES);
    ARCS_CHECK(!mm_blob.empty(), "MISMATCHES blob missing");

    const uint8_t* mm_ptr = mm_blob.data();
    const uint8_t* mm_end = mm_ptr + mm_blob.size();

    // 7. Load quality stream
    auto q_blob = rdr.read_blob(BlobType::QUALITY_DATA);
    ARCS_CHECK(!q_blob.empty(), "QUALITY_DATA blob missing");
    const uint8_t* q_ptr = q_blob.data();
    const uint8_t* q_end = q_ptr + q_blob.size();

    // 8. Decode each mapped read
    // read_len is not stored in v1 WGS-SE archives; hardcoded per-record below.
    (void)hdr.genome_len; // genome_len repurposed in MST archives; unused here

    // Decode mismatches and reconstruct sequences
    for (size_t i = 0; i < n_reads; ++i) {
        // Read mismatch count
        size_t n_mm = (size_t)read_varint(mm_ptr, mm_end);

        std::vector<uint32_t> mm_offsets;
        std::vector<uint8_t>  mm_bases;
        uint32_t prev_off = 0;
        for (size_t j = 0; j < n_mm; ++j) {
            uint32_t delta = (uint32_t)read_varint(mm_ptr, mm_end);
            mm_offsets.push_back(prev_off + delta);
            prev_off = mm_offsets.back();
        }

        // Decode base bytes
        if (n_mm > 0) {
            size_t base_len = (size_t)read_varint(mm_ptr, mm_end);
            ARCS_CHECK(mm_ptr + base_len <= mm_end, "MM base bytes overflow");

            RansDecoder dec;
            const uint8_t* bp  = mm_ptr;
            const uint8_t* ben = mm_ptr + base_len;
            dec.init(bp); bp += 4;

            // We need read_len to compute gc_frac... store in container TODO
            float gc_frac = 0.5f; // placeholder
            for (size_t j = 0; j < n_mm; ++j) {
                ContextKey ctx = mismatch_ctx((int)mm_offsets[j], 150, gc_frac);
                uint8_t base   = mm_model.decode_sym(dec, bp, ben, ctx);
                mm_bases.push_back(base);
            }
            mm_ptr += base_len;
        }

        // Reconstruct sequence
        std::string seq, qual;
        if (read_mapped[i]) {
            pos_t pos = read_pos[i];
            // read_len is unknown here — TODO: store in header
            // Workaround: read until mismatch positions give us a bound
            // This is a v1 limitation; will fix in v2
            int rlen = 150; // default; must be stored in container
            if (pos + rlen > genome.size()) rlen = (int)(genome.size() - pos);
            std::string ref_substr = genome.substr(pos, rlen);
            seq = reconstruct_seq(ref_substr, mm_offsets, mm_bases, strand[i]);
        }

        // Decode quality
        // For now: read quality size from stream and decode
        if (q_ptr < q_end) {
            size_t qlen = (size_t)read_varint(q_ptr, q_end);
            if (q_ptr + qlen <= q_end) {
                // Simple LZMA decode for quality (full block)
                // Context rANS quality decode is more complex — placeholder
                // TODO: implement full context rANS quality decode
                qual = std::string(seq.size(), 'I'); // placeholder quality
                q_ptr += qlen;
            }
        }

        if (seq.empty()) {
            seq  = std::string(150, 'N');
            qual = std::string(150, '!');
        }

        Read r;
        r.name = names[i];
        r.seq  = seq;
        r.qual = qual;
        out_writer.write(r);
    }

    // 9. Append unmapped reads
    if (rdr.has_blob(BlobType::UNMAPPED)) {
        auto unmapped_blob = rdr.read_blob(BlobType::UNMAPPED);
        if (!unmapped_blob.empty()) {
            auto raw = arcs_decompress(unmapped_blob.data(), unmapped_blob.size());
            // Write raw FASTQ directly
            std::string fastq_str(raw.begin(), raw.end());
            std::istringstream ss(fastq_str);
            std::string line;
            Read r;
            int line_num = 0;
            while (std::getline(ss, line)) {
                switch (line_num % 4) {
                    case 0: r.name = line.substr(1); break;
                    case 1: r.seq  = line; break;
                    case 3: r.qual = line;
                            out_writer.write(r); break;
                }
                ++line_num;
            }
        }
    }
}

// ── decode_amplicon ───────────────────────────────────────────────────────────
// Dispatches to the byte-exact V7* path (QUALITY_PERM present) or falls back
// to the legacy V7 path (COUNT_DATA only, quality-lossy) or legacy pre-dedup.
//
// V7* BYTE-EXACT PATH:
//   QUALITY_PERM present → perm[p] = original read index at cluster pos p.
//   QUALITY_DATA = flat concatenation of all quality strings in cluster order.
//   Decoder recovers output[perm[p]].qual = flat_qual[p*L .. (p+1)*L - 1].
//   Sequences: output[perm[p]].seq  = cluster_seqs[cluster_of(p)].
//   Names:     output[perm[p]].name = names[perm[p]] (stored in original order).
//
// V7 LEGACY PATH (COUNT_DATA only, no QUALITY_PERM):
//   All reads in a cluster share the first quality (quality-lossy).
//   Used for archives produced by the old encoder. Not used for new writes.
void ARCSDecoder::decode_amplicon(const ARCSReader& rdr, FASTQWriter& out_writer) {
    const auto& hdr  = rdr.header();
    size_t n_reads   = hdr.n_reads;
    bool   is_v7star = rdr.has_blob(BlobType::QUALITY_PERM);
    bool   is_v7     = rdr.has_blob(BlobType::COUNT_DATA);

    // ── Decompress unique sequences (BSC or LZMA-9, flagged) ─────────────────
    auto seq_blob = rdr.read_blob(BlobType::GENOME);
    std::vector<uint8_t> seq_raw;

    if ((is_v7 || is_v7star) && !seq_blob.empty()) {
        uint8_t flag = seq_blob[0];
        const uint8_t* payload = seq_blob.data() + 1;
        size_t payload_len = seq_blob.size() - 1;
        if (flag == 0x02) {
            // BSC decompress via subprocess
            char tmp_in[256], tmp_out[256];
            snprintf(tmp_in,  sizeof(tmp_in),  "arcs_dec_tmp.bsc");
            snprintf(tmp_out, sizeof(tmp_out), "arcs_dec_tmp.bin");
            FILE* f = fopen(tmp_in, "wb");
            if (f) { fwrite(payload, 1, payload_len, f); fclose(f); }
#ifdef _WIN32
            std::string bsc_null = " 2>nul";
#else
            std::string bsc_null = " 2>/dev/null";
#endif
            std::string cmd = std::string("bsc d ") + tmp_in + " " + tmp_out + bsc_null;
            if (system(cmd.c_str()) == 0) {
                FILE* g = fopen(tmp_out, "rb");
                if (g) {
                    fseek(g, 0, SEEK_END);
                    size_t sz = (size_t)ftell(g);
                    fseek(g, 0, SEEK_SET);
                    seq_raw.resize(sz);
                    if (fread(seq_raw.data(), 1, sz, g) != sz) seq_raw.clear();
                    fclose(g);
                }
                remove(tmp_out);
            }
            remove(tmp_in);
        } else {
            seq_raw = arcs_decompress(payload, payload_len);
        }
    } else if (!seq_blob.empty()) {
        seq_raw = arcs_decompress(seq_blob.data(), seq_blob.size());
    }

    // ── Parse unique sequences ────────────────────────────────────────────────
    std::vector<std::string> unique_seqs;
    {
        std::string seq_str(seq_raw.begin(), seq_raw.end());
        std::istringstream ss(seq_str);
        std::string line;
        while (std::getline(ss, line))
            if (!line.empty() && line[0] != '>') unique_seqs.push_back(line);
    }

    // ── Names (original input order) ──────────────────────────────────────────
    std::vector<std::string> names;
    if (rdr.has_blob(BlobType::NAMES))
        names = decode_names_blob(rdr.read_blob(BlobType::NAMES), n_reads);

    // ── Dedup counts ──────────────────────────────────────────────────────────
    std::vector<uint32_t> counts;
    if (is_v7 || is_v7star) {
        auto cb = rdr.read_blob(BlobType::COUNT_DATA);
        auto cr = arcs_decompress(cb.data(), cb.size());
        size_t nu = cr.size() / 4;
        counts.resize(nu);
        for (size_t i = 0; i < nu; ++i) {
            counts[i] = (uint32_t)cr[i*4+0]         |
                        ((uint32_t)cr[i*4+1] <<  8)  |
                        ((uint32_t)cr[i*4+2] << 16)  |
                        ((uint32_t)cr[i*4+3] << 24);
        }
    } else {
        counts.assign(unique_seqs.size(), 1);
    }

    // ── V7* BYTE-EXACT path ───────────────────────────────────────────────────
    if (is_v7star) {
        // Determine read length L from first unique sequence
        int L = unique_seqs.empty() ? 0 : (int)unique_seqs[0].size();

        // Read permutation: perm[p] = original read index at cluster pos p
        auto perm_blob = rdr.read_blob(BlobType::QUALITY_PERM);
        auto perm_raw  = arcs_decompress(perm_blob.data(), perm_blob.size());
        size_t n_perm  = perm_raw.size() / 4;
        std::vector<uint32_t> perm(n_perm);
        for (size_t i = 0; i < n_perm; ++i) {
            perm[i] = (uint32_t)perm_raw[i*4+0]         |
                      ((uint32_t)perm_raw[i*4+1] <<  8)  |
                      ((uint32_t)perm_raw[i*4+2] << 16)  |
                      ((uint32_t)perm_raw[i*4+3] << 24);
        }

        // Read flat quality: n_reads × L bytes in cluster order
        auto q_blob = rdr.read_blob(BlobType::QUALITY_DATA);
        auto q_raw  = arcs_decompress(q_blob.data(), q_blob.size());

        // Use n_perm as the authoritative read count (equals n_reads for valid archives)
        size_t n = (n_reads > 0) ? n_reads : n_perm;
        if (n_perm < n) n = n_perm; // safety: never exceed what we have in perm

        // Allocate output array — indexed by original read index
        std::vector<Read> output(n);
        // Assign names by original index
        for (size_t i = 0; i < n; ++i) {
            output[i].name = (i < names.size()) ? names[i] : ("read_" + std::to_string(i));
        }

        // Assign sequences and qualities via permutation
        // Also need to map cluster position p → cluster index → sequence
        size_t cluster_pos = 0; // current position in flat cluster-sorted sequence
        for (size_t ci = 0; ci < unique_seqs.size() && cluster_pos < n; ++ci) {
            uint32_t cnt = (ci < counts.size()) ? counts[ci] : 1;
            for (uint32_t j = 0; j < cnt && cluster_pos < n; ++j, ++cluster_pos) {
                size_t p    = cluster_pos;
                uint32_t oi = (p < perm.size()) ? perm[p] : (uint32_t)p;
                if (oi >= n) continue; // guard against corrupt archive

                output[oi].seq = unique_seqs[ci];

                // Quality slice: flat_qual[p * L .. (p+1)*L - 1]
                size_t q_off = p * (size_t)L;
                if (L > 0 && q_off + L <= q_raw.size()) {
                    output[oi].qual.assign(
                        (const char*)q_raw.data() + q_off, L);
                } else {
                    output[oi].qual.assign(L, '!');
                }
            }
        }

        // Write in original order
        for (size_t i = 0; i < n; ++i) {
            if (output[i].seq.empty()) {
                output[i].seq.assign(L > 0 ? L : 1, 'N');
                output[i].qual.assign(L > 0 ? L : 1, '!');
            }
            out_writer.write(output[i]);
        }
        return;
    }

    // ── V7 LEGACY path (quality-lossy, one quality per cluster) ──────────────
    // Used only when decoding archives produced by the old encoder.
    {
        auto q_blob = rdr.read_blob(BlobType::QUALITY_DATA);
        auto q_raw  = q_blob.empty() ? std::vector<uint8_t>{} :
                      arcs_decompress(q_blob.data(), q_blob.size());

        // Detect V7 residual format: first 2 bytes = read_len
        std::string qual_str;
        bool is_residual_fmt = false;
        std::vector<uint8_t> q_mean_arr;
        int residual_L = 0;
        size_t residual_data_off = 0;

        if (is_v7 && q_raw.size() >= 3) {
            int rl = (int)q_raw[0] | ((int)q_raw[1] << 8);
            if (rl >= 20 && rl <= 1000 && (size_t)(2 + rl) <= q_raw.size()) {
                residual_L        = rl;
                residual_data_off = 2 + rl;
                q_mean_arr.assign(q_raw.begin() + 2, q_raw.begin() + 2 + rl);
                is_residual_fmt   = true;
            }
        }
        if (!is_residual_fmt) qual_str.assign(q_raw.begin(), q_raw.end());

        size_t n_unique = unique_seqs.size();
        std::vector<std::string> canonical_quals;
        if (is_residual_fmt && !q_mean_arr.empty() && n_unique > 0) {
            canonical_quals.resize(n_unique, std::string(residual_L, 'I'));
            size_t offset = residual_data_off;
            for (int j = 0; j < residual_L; ++j) {
                for (size_t i = 0; i < n_unique; ++i) {
                    if (offset < q_raw.size() && j < (int)canonical_quals[i].size()) {
                        int res = (int)q_raw[offset++] - 32;
                        int qv  = (int)q_mean_arr[j] + res;
                        qv = std::min(std::max(qv, 0), 42);
                        canonical_quals[i][j] = (char)(qv + 33);
                    }
                }
            }
        }

        size_t qual_pos = 0;
        size_t name_idx = 0;
        for (size_t i = 0; i < n_unique; ++i) {
            uint32_t cnt  = (i < counts.size()) ? counts[i] : 1;
            int      qlen = (int)unique_seqs[i].size();
            for (uint32_t c = 0; c < cnt; ++c) {
                Read r;
                r.name = (name_idx < names.size()) ? names[name_idx] :
                         ("read_" + std::to_string(name_idx));
                ++name_idx;
                r.seq = unique_seqs[i];
                if (is_residual_fmt) {
                    r.qual = (i < canonical_quals.size()) ? canonical_quals[i]
                                                           : std::string(qlen, 'I');
                } else {
                    if (qual_pos + (size_t)qlen <= qual_str.size()) {
                        r.qual = qual_str.substr(qual_pos, qlen);
                        qual_pos += qlen;
                    } else {
                        r.qual = std::string(qlen, 'I');
                    }
                }
                out_writer.write(r);
            }
        }
    }
}

// ── decode_wgs_mst ────────────────────────────────────────────────────────────
// Decodes WGS archives written by encode_wgs_mst (MST_TREE + MST_DELTAS blobs).
// Fully lossless: sequences and quality reconstructed exactly.
void ARCSDecoder::decode_wgs_mst(const ARCSReader& rdr, FASTQWriter& out_writer) {
    const auto& hdr = rdr.header();

    // ── 1. Read metadata (read_len + n_reads from STRAND_FLAGS) ──────────────
    int    read_len = 150;
    size_t n_reads  = hdr.n_reads;
    if (rdr.has_blob(BlobType::STRAND_FLAGS)) {
        auto meta = rdr.read_blob(BlobType::STRAND_FLAGS);
        if (meta.size() >= 8) {
            read_len = (int)(((uint32_t)meta[0]<<24)|((uint32_t)meta[1]<<16)|
                              ((uint32_t)meta[2]<<8)|(uint32_t)meta[3]);
            n_reads  = (size_t)(((uint32_t)meta[4]<<24)|((uint32_t)meta[5]<<16)|
                                 ((uint32_t)meta[6]<<8)|(uint32_t)meta[7]);
        }
    }
    if (read_len <= 0) read_len = 150;
    if (n_reads == 0)  n_reads  = hdr.n_reads;

    // ── 2. Decode sequences ───────────────────────────────────────────────────
    // MST_TREE: BSC-compressed parent pointer array.
    // MST_DELTAS: LZMA-compressed metadata (flags, shifts, sub positions, nonoverlap lengths).
    // MST_SEQ_TEXT: BSC-compressed ACGT text stream (root seqs + sub bases + nonoverlap segs).
    auto tree_compressed  = rdr.read_blob(BlobType::MST_TREE);
    auto delta_compressed = rdr.read_blob(BlobType::MST_DELTAS);
    auto rc_blob          = rdr.has_blob(BlobType::MST_RC_FLAGS) ?
                            rdr.read_blob(BlobType::MST_RC_FLAGS) : std::vector<uint8_t>{};
    ARCS_CHECK(!tree_compressed.empty(),  "MST_TREE blob missing");
    ARCS_CHECK(!delta_compressed.empty(), "MST_DELTAS blob missing");

    auto tree_blob  = bsc_decompress_buf(tree_compressed.data(),  tree_compressed.size());
    auto delta_blob = arcs_decompress(delta_compressed.data(), delta_compressed.size());

    // Seq text: present in text-separated format (MST_SEQ_TEXT blob).
    std::vector<uint8_t> seq_text_blob;
    if (rdr.has_blob(BlobType::MST_SEQ_TEXT)) {
        auto stc = rdr.read_blob(BlobType::MST_SEQ_TEXT);
        seq_text_blob = bsc_decompress_buf(stc.data(), stc.size());
    }

    MSTSequenceDecoder mst_dec;
    std::vector<uint32_t>            dfs_order;
    std::vector<uint32_t>            parents;
    std::vector<std::vector<bool>>   dev_sets;
    const std::vector<uint8_t>* seq_text_ptr = seq_text_blob.empty() ? nullptr : &seq_text_blob;
    auto seqs = mst_dec.decode_sequences(
        tree_blob, delta_blob, rc_blob,
        n_reads, read_len, &dfs_order, &parents, &dev_sets, seq_text_ptr);

    // ── 3. Decode quality ─────────────────────────────────────────────────────
    auto q_data       = rdr.has_blob(BlobType::QUALITY_DATA)  ? rdr.read_blob(BlobType::QUALITY_DATA)  : std::vector<uint8_t>{};
    auto q_model_blob = rdr.has_blob(BlobType::QUALITY_MODEL) ? rdr.read_blob(BlobType::QUALITY_MODEL) : std::vector<uint8_t>{};

    // Quality model blob is LZMA-compressed (encoder compresses 210KB model → ~15KB).
    // Decompress if non-empty and not a bare mode byte.
    std::vector<uint8_t> q_model;
    if (q_model_blob.size() > 1) {
        auto decompressed = arcs_decompress(q_model_blob.data(), q_model_blob.size());
        q_model = decompressed.empty() ? q_model_blob : decompressed;
    } else {
        q_model = q_model_blob;
    }

    // Fallback: PE_R2_DELTAS stores quality model in WGS SCS path
    if (q_model.empty() && rdr.has_blob(BlobType::PE_R2_DELTAS))
        q_model = rdr.read_blob(BlobType::PE_R2_DELTAS);

    auto quals = mst_dec.decode_quality(q_data, q_model, dfs_order, parents, dev_sets, n_reads, read_len);

    // ── 4. Read names ─────────────────────────────────────────────────────────
    std::vector<std::string> names;
    if (rdr.has_blob(BlobType::NAMES))
        names = decode_names_blob(rdr.read_blob(BlobType::NAMES), n_reads);
    while (names.size() < n_reads)
        names.push_back("read_" + std::to_string(names.size()));

    // ── 5. Write FASTQ ────────────────────────────────────────────────────────
    for (size_t i = 0; i < n_reads; ++i) {
        Read r;
        r.name = names[i];
        r.seq  = (i < seqs.size()  && !seqs[i].empty())  ? seqs[i]  : std::string(read_len, 'N');
        r.qual = (i < quals.size() && !quals[i].empty()) ? quals[i] : std::string(read_len, 'I');
        out_writer.write(r);
    }
}

// ── decode_wgs_chain ──────────────────────────────────────────────────────────
// Decodes archives written by encode_wgs_chain (CHAIN_VSEQ + MST_DELTAS blobs).
void ARCSDecoder::decode_wgs_chain(const ARCSReader& rdr, FASTQWriter& out_writer) {
    const auto& hdr = rdr.header();

    int    read_len = 150;
    size_t n_reads  = hdr.n_reads;
    if (rdr.has_blob(BlobType::STRAND_FLAGS)) {
        auto meta = rdr.read_blob(BlobType::STRAND_FLAGS);
        if (meta.size() >= 8) {
            read_len = (int)(((uint32_t)meta[0]<<24)|((uint32_t)meta[1]<<16)|
                              ((uint32_t)meta[2]<<8)|(uint32_t)meta[3]);
            n_reads  = (size_t)(((uint32_t)meta[4]<<24)|((uint32_t)meta[5]<<16)|
                                 ((uint32_t)meta[6]<<8)|(uint32_t)meta[7]);
        }
    }
    if (read_len <= 0) read_len = 150;
    if (n_reads == 0)  n_reads  = hdr.n_reads;

    ARCS_CHECK(rdr.has_blob(BlobType::CHAIN_STARTS), "CHAIN_STARTS blob missing");
    auto chain_starts = rdr.read_blob(BlobType::CHAIN_STARTS);

    // Identity sorted_order — reads stored and output in chain order, no permutation.
    std::vector<uint32_t> sorted_order(n_reads);
    std::iota(sorted_order.begin(), sorted_order.end(), 0u);

    ARCS_CHECK(rdr.has_blob(BlobType::CHAIN_VSEQ), "CHAIN_VSEQ blob missing");
    auto vseq_compressed  = rdr.read_blob(BlobType::CHAIN_VSEQ);
    auto vseq_blob = arcs_decompress(vseq_compressed.data(), vseq_compressed.size());

    auto delta_compressed = rdr.read_blob(BlobType::MST_DELTAS);
    ARCS_CHECK(!delta_compressed.empty(), "MST_DELTAS blob missing");
    auto delta_blob = bsc_decompress_buf(delta_compressed.data(), delta_compressed.size());

    auto rc_blob = rdr.has_blob(BlobType::MST_RC_FLAGS)
                   ? rdr.read_blob(BlobType::MST_RC_FLAGS)
                   : std::vector<uint8_t>{};

    ChainSequenceDecoder chain_dec;
    std::vector<uint32_t>          parents;
    std::vector<std::vector<bool>> dev_sets;
    auto seqs = chain_dec.decode_sequences(
        vseq_blob, delta_blob, chain_starts, sorted_order, rc_blob,
        n_reads, read_len, &parents, &dev_sets);

    auto q_data       = rdr.has_blob(BlobType::QUALITY_DATA)
                        ? rdr.read_blob(BlobType::QUALITY_DATA)
                        : std::vector<uint8_t>{};
    auto q_model_blob = rdr.has_blob(BlobType::QUALITY_MODEL)
                        ? rdr.read_blob(BlobType::QUALITY_MODEL)
                        : std::vector<uint8_t>{};

    std::vector<uint8_t> q_model;
    if (q_model_blob.size() > 1) {
        auto decompressed = arcs_decompress(q_model_blob.data(), q_model_blob.size());
        q_model = decompressed.empty() ? q_model_blob : decompressed;
    } else {
        q_model = q_model_blob;
    }

    MSTSequenceDecoder mst_dec;
    auto quals = mst_dec.decode_quality(
        q_data, q_model, sorted_order, parents, dev_sets, n_reads, read_len);

    std::vector<std::string> names;
    if (rdr.has_blob(BlobType::NAMES))
        names = decode_names_blob(rdr.read_blob(BlobType::NAMES), n_reads);
    while (names.size() < n_reads)
        names.push_back("read_" + std::to_string(names.size()));

    for (size_t i = 0; i < n_reads; ++i) {
        Read r;
        r.name = names[i];
        r.seq  = (i < seqs.size()  && !seqs[i].empty())  ? seqs[i]  : std::string(read_len, 'N');
        r.qual = (i < quals.size() && !quals[i].empty()) ? quals[i] : std::string(read_len, 'I');
        out_writer.write(r);
    }
}

// ── decompress_chunked ────────────────────────────────────────────────────────
// Reads the CHUNK_DATA blob from the outer archive, parses the size-prefixed
// manifest, decompresses each inner .arcs archive to a temp file, then
// appends each temp output FASTQ to the final output in chunk order.
static void decompress_chunked_impl(const std::string& input_path,
                                     const std::string& output_path) {
    ARCSReader outer(input_path);
    auto manifest = outer.read_blob(BlobType::CHUNK_DATA);
    ARCS_CHECK(manifest.size() >= 4, "CHUNK_DATA blob too small");

    const uint8_t* p = manifest.data();
    const uint8_t* end = p + manifest.size();

    auto read_be32_p = [](const uint8_t* q) -> uint32_t {
        return ((uint32_t)q[0]<<24)|((uint32_t)q[1]<<16)|((uint32_t)q[2]<<8)|q[3];
    };
    auto read_be64_p = [](const uint8_t* q) -> uint64_t {
        uint64_t hi = ((uint32_t)q[0]<<24)|((uint32_t)q[1]<<16)|((uint32_t)q[2]<<8)|q[3];
        uint64_t lo = ((uint32_t)q[4]<<24)|((uint32_t)q[5]<<16)|((uint32_t)q[6]<<8)|q[7];
        return (hi << 32) | lo;
    };

    uint32_t n_chunks = read_be32_p(p);
    p += 4;

    // Collect each chunk's (pointer, size) — chunks are independent inner archives,
    // so they decode in any order; only the final concatenation must be in order.
    std::vector<const uint8_t*> cptr(n_chunks);
    std::vector<uint64_t>       csz(n_chunks);
    for (uint32_t ci = 0; ci < n_chunks; ++ci) {
        ARCS_CHECK((size_t)(end - p) >= 8, "CHUNK_DATA manifest truncated at chunk size");
        uint64_t sz = read_be64_p(p); p += 8;
        ARCS_CHECK((uint64_t)(end - p) >= sz, "CHUNK_DATA manifest truncated at chunk data");
        cptr[ci] = p; csz[ci] = sz; p += sz;
    }

    { FILE* f = fopen(output_path.c_str(), "wb");
      ARCS_CHECK(f, "Cannot create output FASTQ: " + output_path); fclose(f); }

    // Wave-parallel decode: decompress up to T chunks CONCURRENTLY to per-chunk temp
    // FASTQs, then append them to the output IN ORDER. Peak RAM = O(T × chunk) — the
    // chunked-mode property, now with decode parallelism (chunks are independent).
    unsigned T = std::thread::hardware_concurrency(); if (!T) T = 4;
    if (const char* e = getenv("ARCS_CHUNK_THREADS")) { int v = atoi(e); if (v >= 1 && v <= 128) T = (unsigned)v; }

    auto decode_one = [&](uint32_t ci) {
        char tmp_arcs[256], tmp_fq[256];
        snprintf(tmp_arcs, sizeof(tmp_arcs), "arcs_dec_ck%u_tmp.arcs", ci);
        snprintf(tmp_fq,   sizeof(tmp_fq),   "arcs_dec_ck%u_tmp.fastq", ci);
        { FILE* f = fopen(tmp_arcs, "wb"); ARCS_CHECK(f, "Cannot create temp chunk archive");
          fwrite(cptr[ci], 1, (size_t)csz[ci], f); fclose(f); }
        ARCSDecoder sub_dec;
        sub_dec.decompress(tmp_arcs, tmp_fq);
        remove(tmp_arcs);
    };

    for (uint32_t base = 0; base < n_chunks; base += T) {
        uint32_t hi = std::min<uint32_t>(base + T, n_chunks);
        std::vector<std::thread> th; th.reserve(hi - base);
        for (uint32_t ci = base; ci < hi; ++ci) th.emplace_back(decode_one, ci);
        for (auto& x : th) x.join();
        // Append this wave's temp FASTQs to the output in chunk order.
        for (uint32_t ci = base; ci < hi; ++ci) {
            char tmp_fq[256]; snprintf(tmp_fq, sizeof(tmp_fq), "arcs_dec_ck%u_tmp.fastq", ci);
            FILE* src = fopen(tmp_fq, "rb"); FILE* dst = fopen(output_path.c_str(), "ab");
            ARCS_CHECK(src && dst, "Cannot open files for chunk append");
            static thread_local char buf[65536]; size_t nr;
            while ((nr = fread(buf, 1, sizeof(buf), src)) > 0) fwrite(buf, 1, nr, dst);
            fclose(src); fclose(dst); remove(tmp_fq);
        }
        fprintf(stderr, "[CHUNK] decompressed %u/%u\n", hi, n_chunks);
    }
}

// ── Pseudogenome decompression (chain-pg) ────────────────────────────────────
// Codec flag byte: 0x01=LZMA, 0x03=GeCo3 subprocess, 0x04=ARCS-DNA (embedded).
static std::string decompress_pg(const std::vector<uint8_t>& blob) {
    if (blob.empty()) return {};
    uint8_t flag          = blob[0];
    const uint8_t* payload = blob.data() + 1;
    size_t  plen           = blob.size() - 1;

    if (flag == 0x08) {
        // 2-bit + LZMA-9 fast decode (format 0x08, trial24).
        // Non-ACTG fallback payload is plain LZMA-ASCII (same as vle_decode_pg).
        std::vector<uint8_t> data(payload, payload + plen);
        return pg_decode_2bit(data);
    }
    if (flag == 0x04) {
        // ARCS-DNA embedded decoder (no subprocess, no temp files).
        std::vector<uint8_t> data(payload, payload + plen);
        return dna_decode(data);
    }
    if (flag == 0x05) {
        // Block-parallel ARCS-DNA: [0x05][nb][nb × u32 len][payloads]. Each block is
        // an independent dna_encode stream; decode concurrently, concatenate → exact pg.
        if (blob.size() < 2) return {};
        int nb = blob[1];
        size_t hdr = 2 + (size_t)nb * 4;
        if (nb < 1 || nb > 64 || blob.size() < hdr) return {};
        std::vector<uint32_t> len((size_t)nb); size_t total = hdr;
        for (int b = 0; b < nb; ++b) {
            len[(size_t)b] = (uint32_t)blob[2+b*4] | ((uint32_t)blob[3+b*4]<<8)
                           | ((uint32_t)blob[4+b*4]<<16) | ((uint32_t)blob[5+b*4]<<24);
            total += len[(size_t)b];
        }
        if (total != blob.size()) return {};
        std::vector<size_t> off((size_t)nb + 1); off[0] = hdr;
        for (int b = 0; b < nb; ++b) off[(size_t)b+1] = off[(size_t)b] + len[(size_t)b];
        std::vector<std::string> parts((size_t)nb);
        std::vector<std::thread> th; th.reserve((size_t)nb);
        for (int b = 0; b < nb; ++b)
            th.emplace_back([&, b] {
                std::vector<uint8_t> chunk(blob.begin() + off[(size_t)b], blob.begin() + off[(size_t)b+1]);
                parts[(size_t)b] = dna_decode(chunk);
            });
        for (auto& t : th) t.join();
        std::string pg; for (auto& p : parts) pg += p;
        return pg;
    }
    if (flag == 0x06) {
        // Seeded block-parallel ARCS-DNA: [0x06][nb][nb×u32 len][nb×26 seed_bytes][payloads].
        // Each block's context seed (26 bases) seeds hist/rcomp before decoding.
        static const int SEED_LEN = 26;
        if (blob.size() < 2) return {};
        int nb = blob[1];
        size_t hdr = 2 + (size_t)nb * 4 + (size_t)nb * SEED_LEN;
        if (nb < 1 || nb > 64 || blob.size() < hdr) return {};
        std::vector<uint32_t> len((size_t)nb); size_t total = hdr;
        for (int b = 0; b < nb; ++b) {
            len[(size_t)b] = (uint32_t)blob[2+b*4] | ((uint32_t)blob[3+b*4]<<8)
                           | ((uint32_t)blob[4+b*4]<<16) | ((uint32_t)blob[5+b*4]<<24);
            total += len[(size_t)b];
        }
        if (total != blob.size()) return {};
        size_t seed_base = 2 + (size_t)nb * 4;
        std::vector<std::string> seeds((size_t)nb);
        for (int b = 0; b < nb; ++b)
            seeds[(size_t)b] = std::string(
                blob.begin() + seed_base + (size_t)b * SEED_LEN,
                blob.begin() + seed_base + (size_t)b * SEED_LEN + SEED_LEN);
        std::vector<size_t> off((size_t)nb + 1); off[0] = hdr;
        for (int b = 0; b < nb; ++b) off[(size_t)b+1] = off[(size_t)b] + len[(size_t)b];
        std::vector<std::string> parts((size_t)nb);
        std::vector<std::thread> th; th.reserve((size_t)nb);
        for (int b = 0; b < nb; ++b)
            th.emplace_back([&, b] {
                std::vector<uint8_t> chunk(blob.begin() + off[(size_t)b], blob.begin() + off[(size_t)b+1]);
                parts[(size_t)b] = dna_decode(chunk, seeds[(size_t)b]);
            });
        for (auto& t : th) t.join();
        std::string pg; for (auto& p : parts) pg += p;
        return pg;
    }
    if (flag == 0x07) {
        // Multi-codec multi-block: [0x07][nb][nb×u32 len][nb×26 seed][nb×1 codec_id][payloads]
        // codec_id 0x00 = ARCS adaptive FCM (seeded); 0x01 = VLE 2-bit + LZMA-9.
        static const int SEED_LEN = 26;
        if (blob.size() < 2) return {};
        int nb = blob[1];
        if (nb < 1 || nb > 64) return {};
        size_t hdr = 2 + (size_t)nb * 4 + (size_t)nb * SEED_LEN + (size_t)nb;
        if (blob.size() < hdr) return {};
        std::vector<uint32_t> len((size_t)nb); size_t total = hdr;
        for (int b = 0; b < nb; ++b) {
            len[(size_t)b] = (uint32_t)blob[2+b*4] | ((uint32_t)blob[3+b*4]<<8)
                           | ((uint32_t)blob[4+b*4]<<16) | ((uint32_t)blob[5+b*4]<<24);
            total += len[(size_t)b];
        }
        if (total != blob.size()) return {};
        size_t seed_base   = 2 + (size_t)nb * 4;
        size_t codec_base  = seed_base + (size_t)nb * SEED_LEN;
        std::vector<std::string> seeds((size_t)nb);
        std::vector<uint8_t>    codec_ids((size_t)nb);
        for (int b = 0; b < nb; ++b) {
            seeds[(size_t)b] = std::string(
                blob.begin() + seed_base + (size_t)b * SEED_LEN,
                blob.begin() + seed_base + (size_t)b * SEED_LEN + SEED_LEN);
            codec_ids[(size_t)b] = blob[codec_base + (size_t)b];
        }
        std::vector<size_t> off((size_t)nb + 1); off[0] = hdr;
        for (int b = 0; b < nb; ++b) off[(size_t)b+1] = off[(size_t)b] + len[(size_t)b];
        std::vector<std::string> parts((size_t)nb);
        std::vector<std::thread> th; th.reserve((size_t)nb);
        for (int b = 0; b < nb; ++b)
            th.emplace_back([&, b] {
                std::vector<uint8_t> chunk(blob.begin() + off[(size_t)b], blob.begin() + off[(size_t)b+1]);
                if (codec_ids[(size_t)b] == 0x01)
                    parts[(size_t)b] = vle_decode_pg(chunk);
                else
                    parts[(size_t)b] = dna_decode(chunk, seeds[(size_t)b]);
            });
        for (auto& t : th) t.join();
        std::string pg; for (auto& p : parts) pg += p;
        return pg;
    }
    if (flag == 0x01) {
        auto raw = arcs_decompress(payload, plen);
        return std::string(raw.begin(), raw.end());
    }
    if (flag == 0x03) {
        // GeDe3 subprocess: writes arcs_pg.seq.de (appends .de to .co name).
        static const char* TMP_CO  = "arcs_pg.seq.co";
        static const char* TMP_OUT = "arcs_pg.seq.de";
#ifdef _WIN32
        const char* NR = " 2>nul";
#else
        const char* NR = " 2>/dev/null";
#endif
        { FILE* f = fopen(TMP_CO, "wb"); if (!f) return {};
          fwrite(payload, 1, plen, f); fclose(f); }
        auto try_gede3 = [&](const std::string& bin) {
            return system((bin + " -v " + TMP_CO + NR).c_str()) == 0;
        };
        // Decoder counterpart of GeCo3 (only for archives written with ARCS_USE_GECO3).
        // Look up GeDe3 on PATH; ARCS_GEDE3_BIN overrides.
        const char* gede_bin = getenv("ARCS_GEDE3_BIN");
        if (!try_gede3(gede_bin ? gede_bin : "GeDe3")) {
            remove(TMP_CO); return {};
        }
        std::string pg;
        FILE* f = fopen(TMP_OUT, "rb");
        if (f) {
            fseek(f, 0, SEEK_END); size_t sz = (size_t)ftell(f); fseek(f, 0, SEEK_SET);
            pg.resize(sz);
            if (fread(&pg[0], 1, sz, f) != sz) pg.clear();
            fclose(f);
        }
        remove(TMP_CO); remove(TMP_OUT);
        return pg;
    }
    return {};
}

// ── decode_wgs_chain_pg ───────────────────────────────────────────────────────
// Decodes archives written by encode_wgs_chain_pg (CHAIN_PG_SEQ/POS/AUX blobs).
void ARCSDecoder::decode_wgs_chain_pg(const ARCSReader& rdr, FASTQWriter& out_writer) {
    const auto& hdr = rdr.header();

    // Read L and n from STRAND_FLAGS (same as chain mode)
    int    L = 150;
    size_t n = hdr.n_reads;
    int    surp_order = 4;   // idea-B surprise FCM order (meta[8]; default for old archives)
    int    perm_mode = 0;    // meta[9]: 0=order-free, 1=perm blob, 2=derive from name index
    int    order_col = -1;   // meta[10]: name column to parse for mode 2
    int    plus_style = 0;   // meta[12]: 0=bare '+', 1=+name, 2=PLUS_LINES blob
    bool   final_nl = true;  // meta[13]: input had a trailing newline
    if (rdr.has_blob(BlobType::STRAND_FLAGS)) {
        auto meta = rdr.read_blob(BlobType::STRAND_FLAGS);
        if (meta.size() >= 8) {
            L = (int)(((uint32_t)meta[0]<<24)|((uint32_t)meta[1]<<16)|
                      ((uint32_t)meta[2]<<8)|(uint32_t)meta[3]);
            n = (size_t)(((uint32_t)meta[4]<<24)|((uint32_t)meta[5]<<16)|
                          ((uint32_t)meta[6]<<8)|(uint32_t)meta[7]);
        }
        if (meta.size() >= 9 && meta[8] >= 1 && meta[8] <= 16) surp_order = (int)meta[8];
        if (meta.size() >= 10) perm_mode = (int)meta[9];
        if (meta.size() >= 11 && meta[10] != 255) order_col = (int)meta[10];
        if (meta.size() >= 12 && meta[11] == 1) out_writer.set_crlf(true);   // reproduce CRLF
        if (meta.size() >= 13) plus_style = (int)meta[12];
        if (meta.size() >= 14) final_nl = (meta[13] == 1);
    }
    if (L <= 0) L = 150;
    if (n == 0) return;

    const bool DEC_TIMING = getenv("ARCS_DEC_TIMING") != nullptr;
    auto _dn = [] { return std::chrono::steady_clock::now(); };
    auto _dp = _dn();
    auto _dmark = [&](const char* w){ if(DEC_TIMING){ auto t=_dn(); fprintf(stderr,"[DEC] %s: %.2fs\n", w, std::chrono::duration<double>(t-_dp).count()); _dp=t; } };

    // ── Pre-read all blobs (serial I/O), then decode pg + names concurrently. ──
    // Quality decode cannot overlap pg: it needs seqs/dev_sets from reconstruction.
    // Parallel gain: names (~0.03s) + pos/aux/qm LZMA (~0.1s) hidden behind pg decode.
    ARCS_CHECK(rdr.has_blob(BlobType::CHAIN_PG_SEQ), "chain-pg: missing CHAIN_PG_SEQ");
    auto _pg_raw   = rdr.read_blob(BlobType::CHAIN_PG_SEQ);
    ARCS_CHECK(rdr.has_blob(BlobType::CHAIN_PG_POS), "chain-pg: missing CHAIN_PG_POS");
    auto _pos_comp = rdr.read_blob(BlobType::CHAIN_PG_POS);
    ARCS_CHECK(rdr.has_blob(BlobType::CHAIN_PG_AUX), "chain-pg: missing CHAIN_PG_AUX");
    auto _aux_comp = rdr.read_blob(BlobType::CHAIN_PG_AUX);
    auto _q_data   = rdr.has_blob(BlobType::QUALITY_DATA)      ? rdr.read_blob(BlobType::QUALITY_DATA)      : std::vector<uint8_t>{};
    auto _qm_comp  = rdr.has_blob(BlobType::QUALITY_MODEL)     ? rdr.read_blob(BlobType::QUALITY_MODEL)     : std::vector<uint8_t>{};
    auto _nm_raw   = rdr.has_blob(BlobType::NAMES)             ? rdr.read_blob(BlobType::NAMES)             : std::vector<uint8_t>{};
    auto _perm_raw = (perm_mode == 1 && rdr.has_blob(BlobType::CHAIN_READ_PERM))
                     ? rdr.read_blob(BlobType::CHAIN_READ_PERM) : std::vector<uint8_t>{};
    auto _plus_raw = (plus_style == 2 && rdr.has_blob(BlobType::PLUS_LINES))
                     ? rdr.read_blob(BlobType::PLUS_LINES)      : std::vector<uint8_t>{};

    // Launch pg decode (dominant ~1.5s) and names decode (~0.03s) concurrently.
    auto pg_fut    = std::async(std::launch::async, [&]{ return decompress_pg(_pg_raw); });
    auto names_fut = std::async(std::launch::async, [&]{ return decode_names_blob(_nm_raw, n); });

    // Main thread: LZMA decompress fast blobs while pg runs in background.
    auto pos_raw = arcs_decompress(_pos_comp.data(), _pos_comp.size());
    auto aux_raw = arcs_decompress(_aux_comp.data(), _aux_comp.size());
    std::vector<uint8_t> q_model_raw;
    if (_qm_comp.size() > 1) {
        auto _dec = arcs_decompress(_qm_comp.data(), _qm_comp.size());
        q_model_raw = _dec.empty() ? _qm_comp : _dec;
    } else {
        q_model_raw = _qm_comp;
    }

    // Wait for pg (pos/aux/qm will have finished; names may already be done).
    std::string pg = pg_fut.get();
    ARCS_CHECK(!pg.empty(), "chain-pg: pg decompression failed");
    _dmark("pg_decode");
    if (DEC_TIMING) {
        fprintf(stderr, "[DEC] pg_blob_bytes: %zu  (FCM-compressed in archive)\n", _pg_raw.size());
        fprintf(stderr, "[DEC] pg_raw_chars:  %zu  (uncompressed pseudogenome)\n", pg.size());
    }
    if (getenv("ARCS_DUMP_PG")) {
        const char* path = getenv("ARCS_DUMP_PG");
        FILE* f = fopen(path, "wb");
        if (f) { fwrite(pg.data(), 1, pg.size(), f); fclose(f); }
        fprintf(stderr, "[DEC] pg dumped to %s (%zu bytes)\n", path, pg.size());
    }

    // ── 2. Parse pg positions ─────────────────────────────────────────────────
    std::vector<uint32_t> pg_pos(n);
    {
        const uint8_t* p = pos_raw.data();
        const uint8_t* e = pos_raw.data() + pos_raw.size();
        int64_t prev = 0;
        for (size_t i = 0; i < n; ++i) {
            uint64_t zz = read_varint(p, e);
            int64_t  d  = (int64_t)(zz >> 1) ^ -(int64_t)(zz & 1);
            prev += d;
            pg_pos[i] = (uint32_t)prev;
        }
    }

    // ── 3. Parse AUX data (already decompressed in parallel pre-read above) ─────
    ARCS_CHECK(aux_raw.size() >= 36, "chain-pg: AUX blob too small");

    auto ru32 = [](const uint8_t* p) -> uint32_t {
        return ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|((uint32_t)p[2]<<8)|(uint32_t)p[3];
    };
    const uint8_t* ah  = aux_raw.data();
    uint32_t sz_rc      = ru32(ah+ 0);
    uint32_t sz_mmcnt   = ru32(ah+ 4);
    uint32_t sz_mmpos   = ru32(ah+ 8);
    uint32_t sz_mmbase  = ru32(ah+12);
    uint32_t sz_Ncnt    = ru32(ah+16);
    uint32_t sz_Npos    = ru32(ah+20);
    uint32_t sz_qmmcnt  = ru32(ah+24);
    uint32_t sz_readlen = ru32(ah+28);
    // Format detection: new archives store sz_Nchar at ah+32 then n at ah+36;
    // old archives (no literal N bytes) store n at ah+32. Distinguish by whether
    // n matches at the new offset.
    uint32_t sz_Nchar = 0;
    size_t   hdr_bytes;
    if (ru32(ah+36) == (uint32_t)n) {          // new 10-field header
        sz_Nchar = ru32(ah+32);
        hdr_bytes = 40;
    } else {                                    // legacy 9-field header
        ARCS_CHECK(ru32(ah+32) == (uint32_t)n, "chain-pg: n_reads mismatch in AUX");
        hdr_bytes = 36;
    }

    const uint8_t* ap = aux_raw.data() + hdr_bytes;
    const uint8_t* col_rc      = ap;  ap += sz_rc;
    const uint8_t* col_mmcnt   = ap;  const uint8_t* col_mmcnt_e  = ap+sz_mmcnt;   ap += sz_mmcnt;
    const uint8_t* col_mmpos   = ap;  ap += sz_mmpos;
    const uint8_t* col_mmbase  = ap;  ap += sz_mmbase;
    const uint8_t* col_Ncnt    = ap;  const uint8_t* col_Ncnt_e   = ap+sz_Ncnt;    ap += sz_Ncnt;
    const uint8_t* col_Npos    = ap;  ap += sz_Npos;
    const uint8_t* col_qmmcnt  = ap;  const uint8_t* col_qmmcnt_e = ap+sz_qmmcnt;  ap += sz_qmmcnt;
    const uint8_t* col_readlen = ap;  ap += sz_readlen;
    const uint8_t* col_Nchar   = ap;  ap += sz_Nchar;
    const uint8_t* col_qmmpos  = ap;
    size_t nchar_i = 0;   // running index into col_Nchar (N entries are in read order)
    // Per-read length accessor (uint16 LE); falls back to global L if absent.
    auto readlen_of = [&](size_t i) -> int {
        if (sz_readlen >= 2*(i+1)) return (int)((uint32_t)col_readlen[2*i] | ((uint32_t)col_readlen[2*i+1]<<8));
        return L;
    };

    const uint8_t* p_mmcnt  = col_mmcnt;
    const uint8_t* p_mmpos  = col_mmpos;
    const uint8_t* p_mmbase = col_mmbase;
    const uint8_t* p_Ncnt   = col_Ncnt;
    const uint8_t* p_Npos   = col_Npos;
    const uint8_t* p_qmmcnt = col_qmmcnt;
    const uint8_t* p_qmmpos = col_qmmpos;

    // ── 4. Reconstruct sequences + build dev_sets for quality context ─────────
    std::vector<std::string>        seqs(n);
    std::vector<std::vector<bool>>  dev_sets(n, std::vector<bool>((size_t)L, false));

    // Idea-B surprise: recompute the pg predictability buckets from the DECODED
    // pg (identical to the encoder's) and index them per read in orig-read frame,
    // exactly mirroring encoder.cpp. Zero stored bytes; used only if the quality
    // stream was encoded with the surprise context (mode 0x05/0x06).
    auto pg_surp = compute_pg_surprise(pg, surp_order, N_QUAL_SURPRISE_BINS);
    std::vector<std::vector<uint8_t>> surprise(n, std::vector<uint8_t>((size_t)L, 0));

    for (size_t i = 0; i < n; ++i) {
        bool     rc      = col_rc[i] != 0;
        uint32_t pos_i   = pg_pos[i];
        const int Li = readlen_of(i);                     // per-read length
        ARCS_CHECK(pos_i + (uint32_t)Li <= (uint32_t)pg.size(),
                   "chain-pg: pg position out of bounds");
        if ((int)dev_sets[i].size() != Li) dev_sets[i].assign((size_t)Li, false);
        if ((int)surprise[i].size() != Li) surprise[i].assign((size_t)Li, 0);

        for (int m = 0; m < Li; ++m) {
            size_t pgpos = (size_t)pos_i + (size_t)m;
            uint8_t b = (pgpos < pg_surp.size()) ? pg_surp[pgpos] : (uint8_t)0;
            int oj = rc ? (Li - 1 - m) : m;
            surprise[i][(size_t)oj] = b;
        }

        // Extract pg slice → enc_seq
        std::string enc_seq(pg.data() + pos_i, (size_t)Li);

        // Apply sequence mismatches (enc_seq frame). Each mismatch is also a
        // quality-deviation position (orig-read frame) — derive dev_sets here
        // instead of reading a stored qmm stream.
        uint32_t mm_cnt = (uint32_t)read_varint(p_mmcnt, col_mmcnt_e);
        for (uint32_t k = 0; k < mm_cnt; ++k) {
            uint16_t mpos = (uint16_t)((uint32_t)p_mmpos[0] | ((uint32_t)p_mmpos[1]<<8));
            p_mmpos += 2;
            uint8_t  mbase = *p_mmbase++;
            if (mpos < (uint16_t)Li) {
                enc_seq[mpos] = BASE_TO_CHAR[mbase & 3];
                uint16_t qp = (uint16_t)(rc ? (Li - 1 - mpos) : mpos);
                if (qp < (uint16_t)Li) dev_sets[i][(size_t)qp] = true;
            }
        }

        // Restore non-ACGT positions (enc_seq frame). New archives store the literal
        // original byte (N, lowercase, IUPAC, …) in col_Nchar → byte-lossless. Legacy
        // archives (sz_Nchar==0) restore 'N' as before.
        uint32_t N_cnt = (uint32_t)read_varint(p_Ncnt, col_Ncnt_e);
        for (uint32_t k = 0; k < N_cnt; ++k) {
            uint16_t npos = (uint16_t)((uint32_t)p_Npos[0] | ((uint32_t)p_Npos[1]<<8));
            p_Npos += 2;
            char ch = (sz_Nchar && nchar_i < sz_Nchar) ? (char)col_Nchar[nchar_i] : 'N';
            ++nchar_i;
            if (npos < (uint16_t)Li) enc_seq[npos] = ch;
        }

        // Apply RC (enc_seq → orig_seq)
        seqs[i] = rc ? reverse_complement(enc_seq) : std::move(enc_seq);

        // Build dev_sets[i] from quality-mm positions (orig-read frame)
        uint32_t qmm_cnt = (uint32_t)read_varint(p_qmmcnt, col_qmmcnt_e);
        for (uint32_t k = 0; k < qmm_cnt; ++k) {
            uint16_t qpos = (uint16_t)((uint32_t)p_qmmpos[0] | ((uint32_t)p_qmmpos[1]<<8));
            p_qmmpos += 2;
            if (qpos < (uint16_t)Li) dev_sets[i][(size_t)qpos] = true;
        }
    }
    // Per-read lengths (chain-position order) for the quality path.
    std::vector<int> rlens(n);
    for (size_t i = 0; i < n; ++i) rlens[i] = readlen_of(i);

    _dmark("reconstruct");
    // ── 5. Decode quality (blobs pre-read; q_model_raw decompressed above) ──────
    const auto& q_data = _q_data;   // alias to pre-read blob

    // sorted_order = identity, parents = UINT32_MAX (no tree for chain-pg)
    std::vector<uint32_t> sorted_order(n);
    std::iota(sorted_order.begin(), sorted_order.end(), 0u);
    std::vector<uint32_t> parents(n, UINT32_MAX);

    std::vector<std::string> quals;
    uint8_t qmode = q_model_raw.empty() ? 0x00 : q_model_raw[0];
    if (qmode == 0x07) {
        // Adaptive CM quality coder (fqzcomp-style). No transmitted model.
        // Per-read length (rlens) supports variable-length reads.
        std::vector<std::vector<uint8_t>> rq_out(n);
        for (size_t i = 0; i < n; ++i) rq_out[i].assign((size_t)rlens[i], 0);
        bool ok = qual_cm_decode(q_data, sorted_order, dev_sets, L, rq_out, &seqs, &rlens);
        ARCS_CHECK(ok, "chain-pg: adaptive-CM quality decode failed");
        quals.assign(n, std::string());
        for (size_t i = 0; i < n; ++i) {
            quals[i].resize((size_t)rlens[i]);
            for (int j = 0; j < rlens[i]; ++j)
                quals[i][(size_t)j] = (char)(rq_out[i][(size_t)j] + 33);
        }
    } else {
        MSTSequenceDecoder mst_dec;
        quals = mst_dec.decode_quality(
            q_data, q_model_raw, sorted_order, parents, dev_sets, n, L, &surprise);
    }

    _dmark("quality_decode");
    // ── 6. Names (decoded concurrently with pg; join here) ───────────────────
    std::vector<std::string> names = names_fut.get();
    while (names.size() < n)
        names.push_back("read_" + std::to_string(names.size()));
    _dmark("names_decode");

    // ── 7. Write FASTQ ─────────────────────────────────────────────────────────
    // Order-preserving (default): recover perm[i] = original FASTQ index of assembly
    // position i, then scatter records into file order for a byte-exact result.
    //   mode 1: perm stored explicitly in CHAIN_READ_PERM.
    //   mode 2: file was sorted by a monotonic name-index column, so the order is
    //           the rank of that column's value — parse it from the names (which are
    //           already stored) and rank; no perm blob needed (free on SRA data).
    //   mode 0 / legacy: emit assembly order (only sorted-equal to input).
    std::vector<uint32_t> perm;
    if (perm_mode == 1 && !_perm_raw.empty()) {
        const auto& perm_blob = _perm_raw;
        auto praw = arcs_decompress(perm_blob.data(), perm_blob.size());
        perm.resize(n);
        size_t p = 0; int64_t prev = 0; bool ok = true;
        for (size_t i = 0; i < n; ++i) {
            uint64_t zz = 0; int shift = 0;
            for (;;) {
                if (p >= praw.size()) { ok = false; break; }
                uint8_t b = praw[p++];
                zz |= (uint64_t)(b & 0x7f) << shift;
                if (!(b & 0x80)) break;
                shift += 7;
            }
            if (!ok) break;
            int64_t d = (int64_t)(zz >> 1) ^ -(int64_t)(zz & 1);
            prev += d;
            perm[i] = (uint32_t)prev;
        }
        ARCS_CHECK(ok, "chain-pg: CHAIN_READ_PERM truncated");
    } else if (perm_mode == 2 && order_col >= 0) {
        // Parse the monotonic index column from each name (digit-run tokenization,
        // matching the encoder), then rank the values → original file position.
        std::vector<std::pair<int64_t,uint32_t>> kv(n);
        for (size_t i = 0; i < n; ++i) {
            const std::string& s = names[i];
            int col = 0; size_t j = 0; int64_t val = -1; bool got = false;
            while (j < s.size()) {
                bool dig = (s[j] >= '0' && s[j] <= '9');
                size_t k = j + 1;
                while (k < s.size() && ((s[k] >= '0' && s[k] <= '9') == dig)) ++k;
                if (col == order_col) { if (dig) { val = strtoll(s.substr(j, k-j).c_str(), nullptr, 10); got = true; } break; }
                ++col; j = k;
            }
            ARCS_CHECK(got, "chain-pg: order index column not found in name");
            kv[i] = { val, (uint32_t)i };
        }
        std::sort(kv.begin(), kv.end());
        perm.resize(n);
        for (size_t r = 0; r < n; ++r) perm[kv[r].second] = (uint32_t)r;   // rank = file pos
    }

    // Reconstruct the '+' line suffix per read (assembly order) from plus_style.
    std::vector<std::string> plus_asm;
    if (plus_style == 1) {
        plus_asm.resize(n);
        for (size_t i = 0; i < n; ++i) plus_asm[i] = names[i];       // +name (before names moved)
    } else if (plus_style == 2 && !_plus_raw.empty()) {
        auto raw = arcs_decompress(_plus_raw.data(), _plus_raw.size());
        plus_asm.resize(n);
        size_t pos = 0;
        for (size_t i = 0; i < n; ++i) {
            size_t nl = pos; while (nl < raw.size() && raw[nl] != '\n') ++nl;
            plus_asm[i].assign((const char*)raw.data() + pos, nl - pos);
            pos = (nl < raw.size()) ? nl + 1 : nl;
        }
    }

    if (!perm.empty()) {
        std::vector<Read> out((size_t)n);
        for (size_t i = 0; i < n; ++i) {
            uint32_t d = perm[i];
            ARCS_CHECK(d < n, "chain-pg: permutation index out of range");
            out[d].name = std::move(names[i]);
            out[d].seq  = std::move(seqs[i]);
            out[d].qual = std::move(quals[i]);
            if (!plus_asm.empty()) out[d].plus = std::move(plus_asm[i]);
        }
        for (size_t i = 0; i < n; ++i) out_writer.write(out[i], i + 1 < n ? true : final_nl);
    } else {
        for (size_t i = 0; i < n; ++i) {
            Read r;
            r.name = std::move(names[i]);
            r.seq  = std::move(seqs[i]);
            r.qual = std::move(quals[i]);
            if (!plus_asm.empty()) r.plus = std::move(plus_asm[i]);
            out_writer.write(r, i + 1 < n ? true : final_nl);
        }
    }
    _dmark("write");
}

// ── decompress ───────────────────────────────────────────────────────────────
void ARCSDecoder::decompress(const std::string& input_path,
                              const std::string& output_path) {
    // Detect chunked outer archive before opening a regular reader
    {
        FILE* probe = fopen(input_path.c_str(), "rb");
        ARCS_CHECK(probe, "Cannot open: " + input_path);
        uint8_t hdr_buf[64] = {0};
        fread(hdr_buf, 1, 64, probe);
        fclose(probe);
        bool chunked = (hdr_buf[5] & 0x40) != 0;
        if (chunked) {
            decompress_chunked_impl(input_path, output_path);
            return;
        }
    }

    ARCSReader rdr(input_path);
    FASTQWriter out_writer(output_path);

    bool is_chain_pg  = rdr.has_blob(BlobType::CHAIN_PG_SEQ);
    bool is_amplicon  = !is_chain_pg && (rdr.header().flags & 0x10) != 0;
    bool is_mst       = (rdr.header().flags & 0x20) != 0;
    bool is_chain     = !is_chain_pg && !is_amplicon &&
                        ((rdr.header().flags & 0x80) != 0 || rdr.has_blob(BlobType::CHAIN_STARTS));
    if (is_chain_pg)
        decode_wgs_chain_pg(rdr, out_writer);
    else if (is_amplicon)
        decode_amplicon(rdr, out_writer);
    else if (is_chain)
        decode_wgs_chain(rdr, out_writer);
    else if (is_mst)
        decode_wgs_mst(rdr, out_writer);
    else
        decode_wgs(rdr, out_writer);

    out_writer.flush();
}

// ── print_info ────────────────────────────────────────────────────────────────
void ARCSDecoder::print_info(const std::string& input_path) {
    ARCSReader rdr(input_path);
    const auto& hdr = rdr.header();

    static const char* regime_names[] = {"AMPLICON","TARGETED","WGS","METAGENOMIC"};
    printf("=== ARCS Archive Info ===\n");
    printf("Version:     %d\n",   hdr.version);
    printf("Reads:       %llu\n", (unsigned long long)hdr.n_reads);
    printf("Mapped:      %llu\n", (unsigned long long)hdr.n_mapped);
    printf("Genome len:  %llu bp\n", (unsigned long long)hdr.genome_len);
    printf("AKC score:   %.4f\n", hdr.akc_score);
    printf("Regime:      %s\n",   regime_names[std::min(hdr.regime, (uint8_t)3)]);
    printf("k-mer size:  %d\n",   hdr.k);
    printf("Flags:       0x%02X\n", hdr.flags);

    // Print blob sizes — enumerate all known BlobType values
    printf("\nBlob sizes:\n");
    static const struct { BlobType type; const char* name; } known_blobs[] = {
        { BlobType::GENOME,          "GENOME"          },
        { BlobType::NAMES,           "NAMES"           },
        { BlobType::SE_POSITIONS,    "SE_POSITIONS"    },
        { BlobType::PE_R2_DELTAS,    "PE_R2_DELTAS"    },
        { BlobType::MISMATCHES,      "MISMATCHES"      },
        { BlobType::QUALITY_MODEL,   "QUALITY_MODEL"   },
        { BlobType::QUALITY_DATA,    "QUALITY_DATA"    },
        { BlobType::UNMAPPED,        "UNMAPPED"        },
        { BlobType::GUTENSOR,        "GUTENSOR"        },
        { BlobType::STRAND_FLAGS,    "STRAND_FLAGS"    },
        { BlobType::PE_STATS,        "PE_STATS"        },
        { BlobType::MST_TREE,        "MST_TREE"        },
        { BlobType::MST_DELTAS,      "MST_DELTAS"      },
        { BlobType::MST_RC_FLAGS,    "MST_RC_FLAGS"    },
        { BlobType::COUNT_DATA,      "COUNT_DATA"      },
        { BlobType::QUALITY_PERM,    "QUALITY_PERM"    },
        { BlobType::CHUNK_DATA,      "CHUNK_DATA"      },
        { BlobType::CHAIN_STARTS,    "CHAIN_STARTS"    },
        { BlobType::CHAIN_READ_PERM, "CHAIN_READ_PERM" },
        { BlobType::MST_TREE_DFS,    "MST_TREE_DFS"    },
        { BlobType::CHAIN_VSEQ,      "CHAIN_VSEQ"      },
        { BlobType::CHAIN_PG_SEQ,    "CHAIN_PG_SEQ"    },
        { BlobType::CHAIN_PG_POS,    "CHAIN_PG_POS"    },
        { BlobType::CHAIN_PG_AUX,    "CHAIN_PG_AUX"    },
        { BlobType::PLUS_LINES,      "PLUS_LINES"      },
        { BlobType::MST_SEQ_TEXT,    "MST_SEQ_TEXT"    },
    };
    for (const auto& kb : known_blobs) {
        if (rdr.has_blob(kb.type)) {
            auto blob = rdr.read_blob(kb.type);
            printf("  %-20s %7zu bytes\n", kb.name, blob.size());
        }
    }
}
