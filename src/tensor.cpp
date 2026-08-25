#include "tensor.h"
#include "container.h"
#include "position_enc.h"
#include "dna_coder.h"
#include "common.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <thread>
#include <algorithm>
#include <numeric>

// ── decode_pg_blob: mirrors decompress_pg() from decoder.cpp ─────────────────
// Handles all chain-pg format flags: 0x01=LZMA, 0x04=ARCS-FCM, 0x05=block-par
// FCM, 0x06=seeded-block-par FCM, 0x07=multi-codec, 0x08=2bit+LZMA.
// Returns empty string on failure (flag 0x03=GeCo3 subprocess not supported here).
static std::string decode_pg_blob(const std::vector<uint8_t>& blob) {
    if (blob.empty()) return {};
    uint8_t flag = blob[0];
    const uint8_t* payload = blob.data() + 1;
    size_t  plen           = blob.size() - 1;

    if (flag == 0x08) {
        std::vector<uint8_t> data(payload, payload + plen);
        return pg_decode_2bit(data);
    }
    if (flag == 0x04) {
        std::vector<uint8_t> data(payload, payload + plen);
        return dna_decode(data);
    }
    static const int SEED_LEN = 26;
    auto parallel_blocks = [&](bool seeded, bool multi_codec) -> std::string {
        if (blob.size() < 2) return {};
        int nb = (int)blob[1];
        if (nb < 1 || nb > 64) return {};
        size_t extra = seeded ? (size_t)nb * SEED_LEN : 0;
        if (multi_codec) extra += (size_t)nb;
        size_t hdr = 2 + (size_t)nb * 4 + extra;
        if (blob.size() < hdr) return {};
        std::vector<uint32_t> len((size_t)nb); size_t total = hdr;
        for (int b = 0; b < nb; ++b) {
            len[(size_t)b] = (uint32_t)blob[2+b*4] | ((uint32_t)blob[3+b*4]<<8)
                           | ((uint32_t)blob[4+b*4]<<16) | ((uint32_t)blob[5+b*4]<<24);
            total += len[(size_t)b];
        }
        if (total != blob.size()) return {};
        size_t seed_base  = 2 + (size_t)nb * 4;
        size_t codec_base = seed_base + (seeded ? (size_t)nb * SEED_LEN : 0);
        std::vector<std::string> seeds;
        std::vector<uint8_t>    cids;
        if (seeded) {
            seeds.resize((size_t)nb);
            for (int b = 0; b < nb; ++b)
                seeds[(size_t)b] = std::string(blob.begin()+seed_base+(size_t)b*SEED_LEN,
                                               blob.begin()+seed_base+(size_t)b*SEED_LEN+SEED_LEN);
        }
        if (multi_codec) {
            cids.resize((size_t)nb);
            for (int b = 0; b < nb; ++b) cids[(size_t)b] = blob[codec_base + (size_t)b];
        }
        std::vector<size_t> off((size_t)nb + 1); off[0] = hdr;
        for (int b = 0; b < nb; ++b) off[(size_t)b+1] = off[(size_t)b] + len[(size_t)b];
        std::vector<std::string> parts((size_t)nb);
        std::vector<std::thread> th; th.reserve((size_t)nb);
        for (int b = 0; b < nb; ++b)
            th.emplace_back([&, b] {
                std::vector<uint8_t> chunk(blob.begin()+off[(size_t)b], blob.begin()+off[(size_t)b+1]);
                if (multi_codec && cids[(size_t)b] == 0x01)
                    parts[(size_t)b] = vle_decode_pg(chunk);
                else if (seeded)
                    parts[(size_t)b] = dna_decode(chunk, seeds[(size_t)b]);
                else
                    parts[(size_t)b] = dna_decode(chunk);
            });
        for (auto& t : th) t.join();
        std::string pg; for (auto& p : parts) pg += p;
        return pg;
    };
    if (flag == 0x05) return parallel_blocks(false, false);
    if (flag == 0x06) return parallel_blocks(true,  false);
    if (flag == 0x07) return parallel_blocks(true,  true);
    if (flag == 0x01) {
        auto raw = arcs_decompress(payload, plen);
        return std::string(raw.begin(), raw.end());
    }
    return {}; // flag 0x03 (GeCo3 subprocess) not supported in tensor mode
}

// ── AUX column parser (mirrors decoder.cpp parse logic exactly) ───────────────
struct AuxCols {
    std::vector<uint8_t>  rc;         // rc[i] = 0 or 1
    std::vector<uint32_t> mm_counts;  // mm_counts[i] = mismatches for read i
    std::vector<uint16_t> mm_pos;     // flat mismatch positions (enc_seq frame)
    std::vector<uint8_t>  mm_base;    // flat mismatch bases (2-bit ACGT code)
    std::vector<uint16_t> readlen;    // readlen[i] for read i
    size_t n = 0;
};

static AuxCols parse_aux(const std::vector<uint8_t>& raw, size_t n_reads) {
    AuxCols a;
    a.n = n_reads;
    if (raw.size() < 36) return a;

    auto ru32 = [](const uint8_t* p) -> uint32_t {
        return ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|((uint32_t)p[2]<<8)|(uint32_t)p[3];
    };
    const uint8_t* ah = raw.data();
    uint32_t sz_rc      = ru32(ah+ 0);
    uint32_t sz_mmcnt   = ru32(ah+ 4);
    uint32_t sz_mmpos   = ru32(ah+ 8);
    uint32_t sz_mmbase  = ru32(ah+12);
    uint32_t sz_Ncnt    = ru32(ah+16);
    uint32_t sz_Npos    = ru32(ah+20);
    uint32_t sz_qmmcnt  = ru32(ah+24);
    uint32_t sz_readlen = ru32(ah+28);
    uint32_t sz_Nchar   = 0;
    size_t   hdr_bytes;
    if (raw.size() >= 40 && ru32(ah+36) == (uint32_t)n_reads) {
        sz_Nchar  = ru32(ah+32);
        hdr_bytes = 40;
    } else {
        hdr_bytes = 36;
    }

    const uint8_t* ap = raw.data() + hdr_bytes;
    const uint8_t* col_rc_p    = ap;                           ap += sz_rc;
    const uint8_t* col_mmcnt_p = ap; const uint8_t* col_mmcnt_e = ap + sz_mmcnt; ap += sz_mmcnt;
    const uint8_t* col_mmpos_p = ap;                           ap += sz_mmpos;
    const uint8_t* col_mmbase_p= ap;                           ap += sz_mmbase;
    /* skip Ncnt/Npos */                                        ap += sz_Ncnt + sz_Npos;
    /* skip qmmcnt */                                           ap += sz_qmmcnt;
    const uint8_t* col_rl_p    = ap;                           ap += sz_readlen;
    /* skip Nchar / qmmpos — not needed for tensors */

    a.rc.resize(n_reads);
    for (size_t i = 0; i < n_reads && i < sz_rc; ++i) a.rc[i] = col_rc_p[i];

    a.mm_counts.resize(n_reads, 0);
    {
        const uint8_t* p = col_mmcnt_p;
        for (size_t i = 0; i < n_reads; ++i)
            a.mm_counts[i] = (uint32_t)read_varint(p, col_mmcnt_e);
    }

    size_t total_mm = sz_mmpos / 2;
    a.mm_pos.resize(total_mm);
    for (size_t i = 0; i < total_mm; ++i)
        a.mm_pos[i] = (uint16_t)((uint32_t)col_mmpos_p[2*i] | ((uint32_t)col_mmpos_p[2*i+1]<<8));

    a.mm_base.resize(sz_mmbase);
    for (size_t i = 0; i < sz_mmbase; ++i) a.mm_base[i] = col_mmbase_p[i];

    a.readlen.resize(n_reads, 0);
    for (size_t i = 0; i < n_reads && sz_readlen >= 2*(i+1); ++i)
        a.readlen[i] = (uint16_t)((uint32_t)col_rl_p[2*i] | ((uint32_t)col_rl_p[2*i+1]<<8));

    return a;
}

// ── K-mer frequency computation ───────────────────────────────────────────────
// Returns 4^k normalized frequency vector over the pseudogenome.
// Non-ACGT characters are skipped. Result length = 4^k.
static std::vector<double> compute_kmer_freq(const std::string& pg, int k) {
    int sz = 1;
    for (int i = 0; i < k; ++i) sz *= 4;
    std::vector<uint64_t> cnt(sz, 0);

    auto base_idx = [](char c) -> int {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    };

    int L = (int)pg.size();
    for (int i = 0; i <= L - k; ++i) {
        int idx = 0; bool ok = true;
        for (int j = 0; j < k; ++j) {
            int b = base_idx(pg[i+j]);
            if (b < 0) { ok = false; break; }
            idx = idx * 4 + b;
        }
        if (ok) cnt[(size_t)idx]++;
    }

    uint64_t total = 0; for (auto c : cnt) total += c;
    std::vector<double> freq(sz, 0.0);
    if (total > 0)
        for (int i = 0; i < sz; ++i) freq[i] = (double)cnt[i] / total;
    return freq;
}

// ── Contig structure from coverage array ─────────────────────────────────────
struct ContigStats {
    uint32_t n_contigs;   // number of contigs (zero-coverage gaps segment the pg)
    uint32_t N50;         // N50 of contig lengths (bp)
    uint32_t longest;     // longest contig (bp)
    double   cov_p10, cov_p25, cov_p50, cov_p75, cov_p90; // depth percentiles
};

static ContigStats compute_contig_stats(const std::vector<uint32_t>& cov,
                                        size_t pg_len) {
    ContigStats cs{};
    if (cov.empty() || pg_len == 0) return cs;

    // Segment by coverage: run of >0 = contig, 0 = gap
    std::vector<uint32_t> contig_lens;
    uint32_t cur = 0;
    for (size_t i = 0; i < pg_len; ++i) {
        if (cov[i] > 0) {
            ++cur;
        } else if (cur > 0) {
            contig_lens.push_back(cur);
            cur = 0;
        }
    }
    if (cur > 0) contig_lens.push_back(cur);

    cs.n_contigs = (uint32_t)contig_lens.size();
    if (!contig_lens.empty()) {
        std::sort(contig_lens.rbegin(), contig_lens.rend());
        cs.longest = contig_lens[0];
        uint64_t total = 0; for (auto l : contig_lens) total += l;
        uint64_t half = total / 2, acc = 0;
        for (auto l : contig_lens) {
            acc += l;
            if (acc >= half) { cs.N50 = l; break; }
        }
    }

    // Coverage percentiles
    std::vector<uint32_t> cv_sorted(cov.begin(), cov.begin() + pg_len);
    std::sort(cv_sorted.begin(), cv_sorted.end());
    auto pct = [&](double p) -> double {
        size_t idx = (size_t)(p * (pg_len - 1));
        return (double)cv_sorted[idx];
    };
    cs.cov_p10 = pct(0.10); cs.cov_p25 = pct(0.25); cs.cov_p50 = pct(0.50);
    cs.cov_p75 = pct(0.75); cs.cov_p90 = pct(0.90);
    return cs;
}

// ── Mismatch position histogram ───────────────────────────────────────────────
// 10-bin histogram of mismatch positions normalized by read length.
// Captures systematic bias patterns (library prep artifacts, cycle degradation).
static std::vector<double> compute_mm_pos_hist(const std::vector<uint16_t>& mm_pos,
                                               int readlen, int n_bins = 10) {
    std::vector<double> hist(n_bins, 0.0);
    if (mm_pos.empty() || readlen <= 0) return hist;
    for (uint16_t p : mm_pos) {
        int bin = (int)((double)p / readlen * n_bins);
        if (bin >= n_bins) bin = n_bins - 1;
        hist[bin]++;
    }
    double tot = (double)mm_pos.size();
    for (auto& v : hist) v /= tot;
    return hist;
}

// ── MinHash sketch from pseudogenome ─────────────────────────────────────────
// Computes a MinHash sketch of n_hashes canonical k-mers from pg.
// Canonical = min(forward, reverse-complement) for strand-invariant comparison.
// Returns sorted vector of n_hashes smallest hash values.
// Use with k=21 for Mash-compatible Jaccard estimation across samples.
static uint64_t fnv1a64(const char* s, int k) {
    uint64_t h = 14695981039346656037ULL;
    for (int i = 0; i < k; ++i) {
        h ^= (uint8_t)s[i];
        h *= 1099511628211ULL;
    }
    return h;
}

static std::vector<uint64_t> compute_minhash(const std::string& pg, int k = 21, int n_hashes = 1000) {

    int L = (int)pg.size();
    if (L < k) return {};

    // We use a fixed-size min-heap of n_hashes elements
    std::vector<uint64_t> heap;
    heap.reserve(n_hashes + 1);

    auto rc_base = [](char c) -> char {
        switch (c) {
            case 'A': return 'T'; case 'T': return 'A';
            case 'C': return 'G'; case 'G': return 'C';
            default:  return 'N';
        }
    };

    std::string rcbuf(k, ' ');
    for (int i = 0; i <= L - k; ++i) {
        // Check all bases are ACGT
        bool ok = true;
        for (int j = 0; j < k; ++j) {
            char c = pg[i+j];
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { ok = false; break; }
        }
        if (!ok) continue;

        uint64_t hf = fnv1a64(pg.data() + i, k);
        // Compute reverse-complement
        for (int j = 0; j < k; ++j) rcbuf[j] = rc_base(pg[i + k - 1 - j]);
        uint64_t hr = fnv1a64(rcbuf.data(), k);
        uint64_t h  = std::min(hf, hr);  // canonical

        if ((int)heap.size() < n_hashes) {
            heap.push_back(h);
            std::push_heap(heap.begin(), heap.end()); // max-heap
        } else if (h < heap[0]) {
            std::pop_heap(heap.begin(), heap.end());
            heap.back() = h;
            std::push_heap(heap.begin(), heap.end());
        }
    }
    std::sort(heap.begin(), heap.end());
    return heap;
}

// ── JSON writer (no external deps) ───────────────────────────────────────────
static void json_kv_f(FILE* f, const char* key, double val, bool last = false) {
    fprintf(f, "  \"%s\": %.6f%s\n", key, val, last ? "" : ",");
}
static void json_kv_u(FILE* f, const char* key, uint64_t val, bool last = false) {
    fprintf(f, "  \"%s\": %llu%s\n", key, (unsigned long long)val, last ? "" : ",");
}
static void json_kv_s(FILE* f, const char* key, const char* val, bool last = false) {
    fprintf(f, "  \"%s\": \"", key);
    for (const char* p = val; *p; ++p) {
        if (*p == '\\') fprintf(f, "\\\\");
        else if (*p == '"') fprintf(f, "\\\"");
        else fputc(*p, f);
    }
    fprintf(f, "\"%s\n", last ? "" : ",");
}
static void json_kv_arr(FILE* f, const char* key, const std::vector<double>& v, bool last = false) {
    fprintf(f, "  \"%s\": [", key);
    for (size_t i = 0; i < v.size(); ++i)
        fprintf(f, "%.8f%s", v[i], i+1 < v.size() ? "," : "");
    fprintf(f, "]%s\n", last ? "" : ",");
}

// ── Main tensor export ────────────────────────────────────────────────────────
int run_tensor_export(const std::string& archive_path,
                      const std::string& output_prefix,
                      const TensorOptions& opts)
{
    // ── 1. Open archive and read header ──────────────────────────────────────
    ARCSReader rdr(archive_path);
    const ARCSHeader& hdr = rdr.header();
    size_t n        = (size_t)hdr.n_reads;
    size_t n_mapped = (size_t)hdr.n_mapped;
    size_t pg_len   = (size_t)hdr.genome_len;

    fprintf(stderr, "[TENSOR] archive: %s\n", archive_path.c_str());
    fprintf(stderr, "[TENSOR] n_reads=%zu  n_mapped=%zu  pg_len=%zu\n", n, n_mapped, pg_len);

    if (n == 0) {
        fprintf(stderr, "[TENSOR] error: archive has 0 reads\n");
        return 1;
    }

    // ── 2. Read and decompress blobs ─────────────────────────────────────────
    if (!rdr.has_blob(BlobType::CHAIN_PG_SEQ) ||
        !rdr.has_blob(BlobType::CHAIN_PG_POS) ||
        !rdr.has_blob(BlobType::CHAIN_PG_AUX)) {
        fprintf(stderr, "[TENSOR] error: archive is not in chain-pg format (tensor requires chain-pg mode)\n");
        return 1;
    }

    auto pg_comp  = rdr.read_blob(BlobType::CHAIN_PG_SEQ);
    auto pos_comp = rdr.read_blob(BlobType::CHAIN_PG_POS);
    auto aux_comp = rdr.read_blob(BlobType::CHAIN_PG_AUX);

    fprintf(stderr, "[TENSOR] decoding pg (%zu bytes)...\n", pg_comp.size());

    std::string pg = decode_pg_blob(pg_comp);
    if (!pg.empty())
        fprintf(stderr, "[TENSOR] pg decoded: %zu bp\n", pg.size());
    else
        fprintf(stderr, "[TENSOR] pg decode failed (codec 0x%02x not supported); GC unknown\n",
                pg_comp.empty() ? 0 : (unsigned)pg_comp[0]);

    fprintf(stderr, "[TENSOR] decompressing positions and AUX...\n");
    auto pos_raw = arcs_decompress(pos_comp.data(), pos_comp.size());
    auto aux_raw = arcs_decompress(aux_comp.data(), aux_comp.size());

    // ── 3. Decode pg positions (zigzag varint delta) ─────────────────────────
    std::vector<uint32_t> pg_pos(n, 0);
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

    // ── 4. Parse AUX columns ─────────────────────────────────────────────────
    AuxCols aux = parse_aux(aux_raw, n);

    // ── 5. Compute sample-level features ─────────────────────────────────────

    // Read length: use aux.readlen if populated, else header read_len
    int hdr_readlen = (int)hdr.read_len;
    double mean_readlen = 0.0, std_readlen = 0.0;
    {
        size_t cnt = 0; double sum = 0, sum2 = 0;
        for (size_t i = 0; i < n; ++i) {
            int rl = (aux.readlen[i] > 0) ? (int)aux.readlen[i] : hdr_readlen;
            if (rl > 0) { sum += rl; sum2 += (double)rl * rl; ++cnt; }
        }
        if (cnt > 0) {
            mean_readlen = sum / cnt;
            std_readlen  = std::sqrt(std::max(0.0, sum2/cnt - mean_readlen*mean_readlen));
        }
    }

    // Total mismatches and mismatch rate
    uint64_t total_mm = 0;
    for (size_t i = 0; i < n; ++i) total_mm += aux.mm_counts[i];
    double total_bases_mapped = (double)n_mapped * mean_readlen;
    double mismatch_rate = (total_bases_mapped > 0) ? (double)total_mm / total_bases_mapped : 0.0;

    // Mapping rate
    double mapping_rate = (n > 0) ? (double)n_mapped / n : 0.0;

    // RC fraction
    size_t n_rc = 0;
    for (size_t i = 0; i < n; ++i) if (aux.rc[i]) ++n_rc;
    double rc_frac = (n > 0) ? (double)n_rc / n : 0.5;

    // GC content of pseudogenome (available if pg decoded successfully)
    double gc_content = -1.0;
    if (!pg.empty()) {
        size_t gc = 0;
        for (char c : pg) if (c == 'G' || c == 'C') ++gc;
        gc_content = (double)gc / pg.size();
    }

    // Estimated mean coverage: mapped_bases / pg_len
    double coverage_mean_est = (pg_len > 0) ? total_bases_mapped / pg_len : 0.0;

    // Coverage distribution (actual, if pg fits in RAM or binned otherwise)
    double coverage_std = 0.0, coverage_cv = 0.0;
    size_t n_zero_cov = 0;

    // Use actual pg length if decoded, else header genome_len
    if (!pg.empty()) pg_len = pg.size();

    static const size_t MAX_FULL_PG = 100UL * 1024 * 1024; // 100MB
    if (pg_len > 0 && pg_len <= MAX_FULL_PG) {
        // Full coverage array (small genome / bacterial)
        std::vector<uint32_t> cov(pg_len, 0);
        for (size_t i = 0; i < n; ++i) {
            if (!aux.readlen[i] && !hdr_readlen) continue;
            uint32_t pos = pg_pos[i];
            int rl = (aux.readlen[i] > 0) ? (int)aux.readlen[i] : hdr_readlen;
            uint32_t end = std::min((uint32_t)(pos + rl), (uint32_t)pg_len);
            for (uint32_t p = pos; p < end; ++p) cov[p]++;
        }

        // Coverage stats
        double sum = 0, sum2 = 0;
        for (uint32_t c : cov) {
            sum  += c; sum2 += (double)c * c;
            if (c == 0) ++n_zero_cov;
        }
        coverage_mean_est = sum / pg_len;
        coverage_std = std::sqrt(std::max(0.0, sum2/pg_len - coverage_mean_est*coverage_mean_est));
        coverage_cv  = (coverage_mean_est > 0) ? coverage_std / coverage_mean_est : 0.0;

        // Optionally write coverage binary
        if (opts.emit_coverage) {
            std::string cpath = output_prefix + "_coverage.bin";
            FILE* cf = fopen(cpath.c_str(), "wb");
            if (cf) {
                // Header: [uint32 pg_len][uint32 bin_size=1][uint32[] coverage values]
                uint32_t pg_len32 = (uint32_t)pg_len, bin1 = 1;
                fwrite(&pg_len32, 4, 1, cf);
                fwrite(&bin1, 4, 1, cf);
                fwrite(cov.data(), 4, cov.size(), cf);
                fclose(cf);
                fprintf(stderr, "[TENSOR] wrote %s (%zu positions)\n", cpath.c_str(), pg_len);
            }
        }

    } else if (pg_len > MAX_FULL_PG) {
        // Binned coverage for large genomes
        size_t bin_size = pg_len / opts.coverage_bins + 1;
        size_t n_bins   = (pg_len + bin_size - 1) / bin_size;
        std::vector<uint32_t> cov_binned(n_bins, 0);

        for (size_t i = 0; i < n; ++i) {
            uint32_t pos = pg_pos[i];
            int rl = (aux.readlen[i] > 0) ? (int)aux.readlen[i] : hdr_readlen;
            size_t bin = (size_t)(pos / bin_size);
            if (bin < n_bins) cov_binned[bin] += (uint32_t)rl; // approx: add rl to bin
        }
        // Convert read-length sums to approximate depth per bin
        double sum = 0, sum2 = 0;
        for (size_t b = 0; b < n_bins; ++b) {
            double depth = (double)cov_binned[b] / bin_size;
            sum  += depth;
            sum2 += depth * depth;
            if (cov_binned[b] == 0) ++n_zero_cov;
        }
        double mean_bin = sum / n_bins;
        coverage_std = std::sqrt(std::max(0.0, sum2/n_bins - mean_bin*mean_bin));
        coverage_cv  = (coverage_mean_est > 0) ? coverage_std / coverage_mean_est : 0.0;

        if (opts.emit_coverage) {
            std::string cpath = output_prefix + "_coverage.bin";
            FILE* cf = fopen(cpath.c_str(), "wb");
            if (cf) {
                uint32_t n_bins32 = (uint32_t)n_bins, bs32 = (uint32_t)bin_size;
                fwrite(&n_bins32, 4, 1, cf);
                fwrite(&bs32,     4, 1, cf);
                fwrite(cov_binned.data(), 4, cov_binned.size(), cf);
                fclose(cf);
                fprintf(stderr, "[TENSOR] wrote %s (%zu bins, bin_size=%zu)\n", cpath.c_str(), n_bins, bin_size);
            }
        }
    }

    double cov_zero_frac = (pg_len > 0) ? (double)n_zero_cov / pg_len : 0.0;

    // Mismatch position: mean + 10-bin histogram
    double mm_mean_pos = 0.0;
    std::vector<double> mm_pos_hist(10, 0.0);
    if (!aux.mm_pos.empty()) {
        double sum = 0;
        for (uint16_t p : aux.mm_pos) sum += p;
        mm_mean_pos = sum / aux.mm_pos.size();
        int rl_for_hist = (mean_readlen > 0) ? (int)mean_readlen : (hdr_readlen > 0 ? hdr_readlen : 150);
        mm_pos_hist = compute_mm_pos_hist(aux.mm_pos, rl_for_hist);
    }

    // Position spread (diversity of placement over pg)
    double position_entropy = 0.0;
    if (pg_len > 0 && n > 0) {
        uint32_t mn = *std::min_element(pg_pos.begin(), pg_pos.end());
        uint32_t mx = *std::max_element(pg_pos.begin(), pg_pos.end());
        position_entropy = (double)(mx - mn) / pg_len;
    }

    // Contig jump rate (assembly fragmentation proxy)
    size_t n_contig_jumps = 0;
    if (n > 1 && mean_readlen > 0) {
        for (size_t i = 1; i < n; ++i) {
            int64_t delta = (int64_t)pg_pos[i] - (int64_t)pg_pos[i-1];
            if (delta < 0 || delta > (int64_t)(mean_readlen * 3))
                ++n_contig_jumps;
        }
    }
    double contig_jump_rate = (n > 1) ? (double)n_contig_jumps / (n - 1) : 0.0;

    // ── 6. Compute pg-derived rich features (k-mers + contig structure) ───────
    std::vector<double> kmer3_freq, kmer4_freq;
    ContigStats cstats{};

    std::vector<uint64_t> minhash_sketch;
    if (!pg.empty()) {
        fprintf(stderr, "[TENSOR] computing k-mer frequencies and MinHash sketch...\n");
        kmer3_freq    = compute_kmer_freq(pg, 3);             // 64-dim
        kmer4_freq    = compute_kmer_freq(pg, 4);             // 256-dim
        minhash_sketch = compute_minhash(pg, 21, 1000);       // 1000-hash sketch k=21
    }

    // Contig stats: only available when we built full coverage array (pg ≤ 100MB)
    // Re-build coverage if pg is decoded and small enough
    if (!pg.empty() && pg_len <= 100UL * 1024 * 1024) {
        std::vector<uint32_t> cov2(pg_len, 0);
        for (size_t i = 0; i < n; ++i) {
            uint32_t pos = pg_pos[i];
            int rl = (aux.readlen[i] > 0) ? (int)aux.readlen[i] : hdr_readlen;
            if (rl <= 0) continue;
            uint32_t end = std::min((uint32_t)(pos + rl), (uint32_t)pg_len);
            for (uint32_t p2 = pos; p2 < end; ++p2) cov2[p2]++;
        }
        cstats = compute_contig_stats(cov2, pg_len);
    }

    // Write MinHash sketch (binary: n_hashes uint64s, little-endian)
    if (!minhash_sketch.empty()) {
        std::string spath = output_prefix + "_sketch.bin";
        FILE* sf = fopen(spath.c_str(), "wb");
        if (sf) {
            uint32_t nh = (uint32_t)minhash_sketch.size(), kk = 21;
            fwrite(&nh, 4, 1, sf);
            fwrite(&kk, 4, 1, sf);
            fwrite(minhash_sketch.data(), 8, minhash_sketch.size(), sf);
            fclose(sf);
            fprintf(stderr, "[TENSOR] wrote sketch: %s (%u hashes, k=21)\n",
                    spath.c_str(), nh);
        }
    }

    // ── 7. Write sample feature JSON ──────────────────────────────────────────
    std::string jpath = output_prefix + "_sample.json";
    FILE* jf = fopen(jpath.c_str(), "w");
    if (!jf) {
        fprintf(stderr, "[TENSOR] error: cannot write %s\n", jpath.c_str());
        return 1;
    }
    fprintf(jf, "{\n");
    json_kv_s(jf, "archive",          archive_path.c_str());
    json_kv_u(jf, "n_reads",          (uint64_t)n);
    json_kv_u(jf, "n_mapped",         (uint64_t)n_mapped);
    json_kv_u(jf, "genome_len",       (uint64_t)pg_len);
    json_kv_f(jf, "mapping_rate",     mapping_rate);
    json_kv_f(jf, "gc_content",       gc_content);
    json_kv_f(jf, "mismatch_rate",    mismatch_rate);
    json_kv_u(jf, "total_mm",         (uint64_t)total_mm);
    json_kv_f(jf, "mean_readlen",     mean_readlen);
    json_kv_f(jf, "std_readlen",      std_readlen);
    json_kv_f(jf, "rc_frac",          rc_frac);
    json_kv_f(jf, "coverage_mean",    coverage_mean_est);
    json_kv_f(jf, "coverage_std",     coverage_std);
    json_kv_f(jf, "coverage_cv",      coverage_cv);
    json_kv_f(jf, "cov_zero_frac",    cov_zero_frac);
    json_kv_f(jf, "cov_p10",          cstats.cov_p10);
    json_kv_f(jf, "cov_p25",          cstats.cov_p25);
    json_kv_f(jf, "cov_p50",          cstats.cov_p50);
    json_kv_f(jf, "cov_p75",          cstats.cov_p75);
    json_kv_f(jf, "cov_p90",          cstats.cov_p90);
    json_kv_u(jf, "n_contigs",        (uint64_t)cstats.n_contigs);
    json_kv_u(jf, "N50",              (uint64_t)cstats.N50);
    json_kv_u(jf, "longest_contig",   (uint64_t)cstats.longest);
    json_kv_f(jf, "mm_mean_pos",      mm_mean_pos);
    json_kv_f(jf, "position_entropy", position_entropy);
    json_kv_f(jf, "contig_jump_rate", contig_jump_rate);
    if (!mm_pos_hist.empty())
        json_kv_arr(jf, "mm_pos_hist", mm_pos_hist);
    if (!kmer3_freq.empty())
        json_kv_arr(jf, "kmer3_freq", kmer3_freq);
    if (!kmer4_freq.empty())
        json_kv_arr(jf, "kmer4_freq", kmer4_freq, true);
    else {
        // close JSON without trailing comma
        // rewrite last line: patch kmer3 or mm_pos_hist to be last=true
        // simplest: just close with empty last field placeholder
        fprintf(jf, "  \"_eof\": 1\n");
    }
    fprintf(jf, "}\n");
    fclose(jf);
    fprintf(stderr, "[TENSOR] wrote %s\n", jpath.c_str());

    // ── 7. Write per-read TSV ─────────────────────────────────────────────────
    if (opts.emit_reads) {
        std::string rpath = output_prefix + "_reads.tsv";
        FILE* rf = fopen(rpath.c_str(), "w");
        if (rf) {
            fprintf(rf, "read_idx\tpg_pos\trc\tmm_count\tmm_density\treadlen\n");
            for (size_t i = 0; i < n; ++i) {
                int rl = (aux.readlen[i] > 0) ? (int)aux.readlen[i] : hdr_readlen;
                double mm_dens = (rl > 0) ? (double)aux.mm_counts[i] / rl : 0.0;
                fprintf(rf, "%zu\t%u\t%u\t%u\t%.6f\t%d\n",
                        i, pg_pos[i], (unsigned)aux.rc[i],
                        (unsigned)aux.mm_counts[i], mm_dens, rl);
            }
            fclose(rf);
            fprintf(stderr, "[TENSOR] wrote %s (%zu reads)\n", rpath.c_str(), n);
        }
    }

    return 0;
}
