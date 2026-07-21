#include <chrono>
#include <thread>
#include "chain_encoder.h"
#include "hamming_mst.h"
#include "position_enc.h"
#include "container.h"
#include <algorithm>
#include <numeric>
#include <functional>
#include <cassert>
#include <unordered_map>
#include <cstdio>
#include <cstdlib>

// ── ChainSequenceEncoder::encode ─────────────────────────────────────────────
// Greedy linear chain: build k-NN graph, extract chains, delta-encode in chain order.
// NO permutation stored — reads output in chain order (like --mst-dfs in trial10).
ChainEncodeResult ChainSequenceEncoder::encode(const std::vector<Read>& reads) {
    ChainEncodeResult result;
    if (reads.empty()) return result;

    int n = (int)reads.size();
    int L = (int)reads[0].seq.size();
    result.n_reads  = n;
    result.read_len = L;

    int max_mm = cfg_.max_mm < 0 ? std::max(1, L / 10) : cfg_.max_mm;

    // Adaptive max_shift: for long reads the default 50 misses overlaps with
    // large shifts (e.g., 150 bp reads need shift up to ~112).  Use 3/4 of L.
    MSTConfig cfg_adj = cfg_;
    cfg_adj.max_shift = std::max(cfg_.max_shift, 3 * L / 4);

    // ── 1. Extract sequences ──────────────────────────────────────────────────
    std::vector<std::string> seqs(n);
    for (int i = 0; i < n; ++i) seqs[i] = reads[i].seq;

    // ── 2. Build k-NN graph ───────────────────────────────────────────────────
    MSTSequenceEncoder mst_enc(cfg_adj);
    auto knn_edges = mst_enc.get_knn_edges(seqs);

    // Sort edges ascending by distance: best overlaps first
    std::sort(knn_edges.begin(), knn_edges.end(),
              [](const KNNEdge& a, const KNNEdge& b){ return a.dist < b.dist; });

    // ── 3. Greedy chain extraction via Union-Find ─────────────────────────────
    // Each read gets at most 1 predecessor and 1 successor (linear chain path).
    // Union-Find prevents cycles: don't connect u->v if already in the same chain.
    std::vector<int> predecessor(n, -1);
    std::vector<int> successor(n, -1);
    std::vector<int> in_deg(n, 0);
    std::vector<int> out_deg(n, 0);
    UnionFind uf(n);

    float dist_threshold = (float)max_mm / (float)L;

    for (const auto& e : knn_edges) {
        if (e.dist > dist_threshold) break;
        int u = (int)e.u, v = (int)e.v;

        // Try u -> v
        if (out_deg[u] == 0 && in_deg[v] == 0 && uf.find(u) != uf.find(v)) {
            successor[u]   = v;
            predecessor[v] = u;
            out_deg[u]++;
            in_deg[v]++;
            uf.unite(u, v);
            continue;
        }
        // Try v -> u (edge is undirected in k-NN graph)
        if (out_deg[v] == 0 && in_deg[u] == 0 && uf.find(v) != uf.find(u)) {
            successor[v]   = u;
            predecessor[u] = v;
            out_deg[v]++;
            in_deg[u]++;
            uf.unite(v, u);
        }
    }

    // ── 4. Walk chains -> chain_order ─────────────────────────────────────────
    // chain_order[pos] = original read index at chain position pos.
    std::vector<uint32_t> chain_order;
    chain_order.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (predecessor[i] == -1) {
            int cur = i;
            while (cur != -1) {
                chain_order.push_back((uint32_t)cur);
                cur = successor[cur];
            }
        }
    }
    assert((int)chain_order.size() == n);

    // ── 5. Delta-encode in chain order (columnar layout) ─────────────────────
    // We write separate column buffers for each field type so that LZMA can
    // exploit the homogeneity within each column (e.g., shifts are mostly 0,
    // n_subs are mostly 0–2, etc.).  The columns are concatenated with a
    // 10×4-byte size header before LZMA compression.
    //
    // Columns:
    //   col_vseq     : 2-bit packed verbatim sequences (chain-start reads)
    //   col_vN       : N-position varint lists (chain-start reads)
    //   col_drc      : RC-flag bits for delta reads (1 bit per delta, packed)
    //   col_dshift   : shift zigzag-varints (one per delta)
    //   col_dnsubs   : n_subs varints (one per delta)
    //   col_dsubpos  : sub-position varints (flat, all deltas concatenated)
    //   col_dsubbases: 2-bit packed sub-bases (all deltas concatenated)
    //   col_dexlen   : extra-length varints (one per delta)
    //   col_dexbase  : 2-bit packed extra bases (shift != 0)
    //   col_dN       : N-position varint lists (delta reads)

    std::vector<uint32_t>          parents(n, UINT32_MAX);
    std::vector<std::vector<bool>> dev_sets(n, std::vector<bool>(L, false));
    std::vector<bool>              chain_starts_vec(n, false);
    std::vector<bool>              rc_flags_vec(n, false);

    std::vector<uint8_t> col_vseq, col_vN;
    std::vector<uint8_t> col_drc, col_dshift, col_dnsubs;
    std::vector<uint8_t> col_dsubpos, col_dsubbases;
    std::vector<uint8_t> col_dexlen, col_dexbase, col_dN;
    int drc_bits = 0;   // bit count in last col_drc byte

    // Write len bases from seq[offset..] packed 2 bits/base into buf
    auto wpb = [&](std::vector<uint8_t>& buf, const std::string& seq, int offset, int len) {
        int bytes_needed = (len + 3) / 4;
        for (int b = 0; b < bytes_needed; ++b) {
            uint8_t packed = 0;
            for (int j = 0; j < 4; ++j) {
                int rel = b * 4 + j;
                if (rel >= len) break;
                uint8_t base = encode_base(seq[offset + rel]);
                if (base > 3) base = 0;
                packed |= (base << (6 - j * 2));
            }
            buf.push_back(packed);
        }
    };

    auto wnp = [&](std::vector<uint8_t>& buf, const std::string& orig) {
        uint32_t cnt = 0;
        for (int p = 0; p < L; ++p) if (encode_base(orig[p]) >= 4) ++cnt;
        write_varint(buf, cnt);
        uint32_t prev = 0;
        for (int p = 0; p < L; ++p) {
            if (encode_base(orig[p]) >= 4) {
                write_varint(buf, (uint32_t)p - prev);
                prev = (uint32_t)p;
            }
        }
    };

    auto push_rc_bit = [&](bool val) {
        if (drc_bits % 8 == 0) col_drc.push_back(0);
        if (val) col_drc.back() |= (1u << (drc_bits % 8));
        ++drc_bits;
    };

    std::string enc_frame_prev;

    // PG construction state (used only when build_pg_=true)
    uint32_t pg_seg_base   = 0;
    uint32_t pg_seg_offset = 0;

    for (int i = 0; i < n; ++i) {
        uint32_t idx = chain_order[i];
        const std::string& orig_seq = reads[idx].seq;

        bool is_chain_start;
        ReadAlignResult aln;

        if (predecessor[idx] == -1 || enc_frame_prev.empty()) {
            is_chain_start = true;
        } else {
            aln = align_reads(enc_frame_prev, orig_seq, cfg_adj.max_shift, max_mm);
            is_chain_start = (aln.n_subs > max_mm);
        }

        chain_starts_vec[i] = is_chain_start;

        if (is_chain_start) {
            wpb(col_vseq, orig_seq, 0, L);
            wnp(col_vN, orig_seq);
            parents[idx]   = UINT32_MAX;
            enc_frame_prev = orig_seq;
            ++result.n_chain_starts;

        } else {
            bool use_rc   = aln.rc_child;
            int  shift    = aln.shift;
            int  ov_len   = aln.overlap_len;
            int  ch_start = std::max(0, -shift);
            int  extra_len = std::abs(shift);

            for (uint32_t sp : aln.sub_positions) {
                int pos_in_child = ch_start + (int)sp;
                int orig_pos = use_rc ? (L - 1 - pos_in_child) : pos_in_child;
                if (orig_pos >= 0 && orig_pos < L)
                    dev_sets[idx][orig_pos] = true;
            }

            push_rc_bit(use_rc);
            write_varint(col_dshift, zigzag_encode(shift));
            write_varint(col_dnsubs, aln.sub_positions.size());

            uint32_t prev_p = 0;
            for (uint32_t p : aln.sub_positions) {
                write_varint(col_dsubpos, p - prev_p);
                prev_p = p;
            }

            int n_subs = (int)aln.sub_positions.size();
            for (int si = 0; si < n_subs; si += 4) {
                uint8_t packed = 0;
                for (int j = 0; j < 4 && si + j < n_subs; ++j)
                    packed |= ((aln.sub_bases[si + j] & 0x03) << (6 - j * 2));
                col_dsubbases.push_back(packed);
            }

            write_varint(col_dexlen, extra_len);
            if (extra_len > 0) {
                const std::string& enc_seq = use_rc ? reverse_complement(orig_seq) : orig_seq;
                wpb(col_dexbase, enc_seq, shift > 0 ? ov_len : 0, extra_len);
            }

            wnp(col_dN, orig_seq);

            const std::string& enc_seq = use_rc ? reverse_complement(orig_seq) : orig_seq;
            enc_frame_prev = enc_seq;

            parents[idx]    = chain_order[i - 1];
            rc_flags_vec[i] = use_rc;
            ++result.n_deltas;
        }

        // ── PG construction (optional) ──────────────────────────────────────────
        // enc_frame_prev now holds current read's enc_seq (orig_seq for chain starts).
        if (build_pg_) {
            if (is_chain_start) {
                // New pg segment: write orig_seq (chain start = no RC, enc_frame_prev = orig_seq)
                pg_seg_base   = (uint32_t)result.pg.size();
                pg_seg_offset = 0;
                uint32_t N_cnt = 0;
                for (int j = 0; j < L; ++j) {
                    char b = orig_seq[j];
                    if (encode_base(b) >= 4) {
                        result.pg_N_pos_flat.push_back((uint16_t)j);
                        result.pg.push_back('A');
                        ++N_cnt;
                    } else {
                        result.pg.push_back(b);
                    }
                }
                result.pg_pos.push_back(pg_seg_base);
                result.pg_rc.push_back(0);
                result.pg_mm_counts.push_back(0);
                result.pg_N_counts.push_back(N_cnt);
                result.pg_qmm_counts.push_back(0);

            } else {
                bool use_rc = aln.rc_child;
                int  shift  = aln.shift;
                int  ch_start = std::max(0, -shift);
                // enc_frame_prev = enc_seq (already set above)
                const std::string& enc_seq = enc_frame_prev;

                if (shift >= 0) {
                    // Forward extension: overlap at pg[pg_start..pg_start+ov_len-1]
                    int ov_len = L - shift;
                    uint32_t pg_start_i = pg_seg_base + pg_seg_offset + (uint32_t)shift;

                    uint32_t mm_cnt = 0, N_cnt = 0, qmm_cnt = 0;

                    // Overlap region: compare enc_seq[0..ov_len-1] vs pg[pg_start_i..]
                    for (int j = 0; j < ov_len; ++j) {
                        char b = enc_seq[j];
                        if (encode_base(b) >= 4) {
                            result.pg_N_pos_flat.push_back((uint16_t)j);
                            ++N_cnt;
                        } else if (b != result.pg[pg_start_i + (uint32_t)j]) {
                            result.pg_mm_pos_flat.push_back((uint16_t)j);
                            result.pg_mm_base_flat.push_back(encode_base(b));
                            ++mm_cnt;
                            // Quality deviation: convert enc_seq pos → orig_read pos
                            uint16_t qpos = (uint16_t)(use_rc ? (L - 1 - j) : j);
                            result.pg_qmm_pos_flat.push_back(qpos);
                            ++qmm_cnt;
                        }
                    }
                    // Novel region: enc_seq[ov_len..L-1] → append to pg
                    for (int j = ov_len; j < L; ++j) {
                        char b = enc_seq[j];
                        if (encode_base(b) >= 4) {
                            result.pg_N_pos_flat.push_back((uint16_t)j);
                            result.pg.push_back('A');
                            ++N_cnt;
                        } else {
                            result.pg.push_back(b);
                        }
                    }

                    result.pg_pos.push_back(pg_start_i);
                    result.pg_rc.push_back(use_rc ? 1 : 0);
                    result.pg_mm_counts.push_back(mm_cnt);
                    result.pg_N_counts.push_back(N_cnt);
                    result.pg_qmm_counts.push_back(qmm_cnt);
                    pg_seg_offset += (uint32_t)shift;

                } else {
                    // Negative shift: treat as new pg segment (write full enc_seq).
                    // Sequence has no mismatches vs pg (we write it directly).
                    // For quality: record mismatches vs chain parent (aln.sub_positions).
                    pg_seg_base   = (uint32_t)result.pg.size();
                    pg_seg_offset = 0;
                    uint32_t N_cnt = 0, qmm_cnt = 0;
                    for (int j = 0; j < L; ++j) {
                        char b = enc_seq[j];
                        if (encode_base(b) >= 4) {
                            result.pg_N_pos_flat.push_back((uint16_t)j);
                            result.pg.push_back('A');
                            ++N_cnt;
                        } else {
                            result.pg.push_back(b);
                        }
                    }
                    // Quality deviation positions (from chain alignment, orig frame)
                    for (uint32_t sp : aln.sub_positions) {
                        int pos_in_child = ch_start + (int)sp;
                        uint16_t qpos = (uint16_t)(use_rc ? (L - 1 - pos_in_child) : pos_in_child);
                        if (qpos < (uint16_t)L) {
                            result.pg_qmm_pos_flat.push_back(qpos);
                            ++qmm_cnt;
                        }
                    }
                    result.pg_pos.push_back(pg_seg_base);
                    result.pg_rc.push_back(use_rc ? 1 : 0);
                    result.pg_mm_counts.push_back(0);
                    result.pg_N_counts.push_back(N_cnt);
                    result.pg_qmm_counts.push_back(qmm_cnt);
                }
            }
        }
    }

    // ── 5b. Assemble two separate blobs ─────────────────────────────────────
    // BLOB 1 (CHAIN_VSEQ): verbatim sequences + verbatim N-positions.
    //   Random 2-bit DNA does not compress well — isolating it gives LZMA a
    //   clean context for the homogeneous delta fields in blob 2.
    //   Format: [4B sz_vseq][4B sz_vN][col_vseq][col_vN]
    {
        auto push_u32 = [](std::vector<uint8_t>& v, uint32_t x) {
            v.push_back(x>>24); v.push_back((x>>16)&0xFF);
            v.push_back((x>>8)&0xFF); v.push_back(x&0xFF);
        };
        std::vector<uint8_t> vs;
        push_u32(vs, (uint32_t)col_vseq.size()); push_u32(vs, (uint32_t)col_vN.size());
        vs.insert(vs.end(), col_vseq.begin(), col_vseq.end());
        vs.insert(vs.end(), col_vN.begin(),   col_vN.end());
        result.vseq_bytes = std::move(vs);
    }

    // BLOB 2 (MST_DELTAS): all delta fields — RC flags, shifts, n_subs, sub_pos,
    //   sub_bases, extra_len, extra_bases, delta N-positions.
    //   Format: [4B×8 sizes][col_drc][col_dshift][col_dnsubs][col_dsubpos]
    //           [col_dsubbases][col_dexlen][col_dexbase][col_dN]
    {
        auto push_u32 = [](std::vector<uint8_t>& v, uint32_t x) {
            v.push_back(x>>24); v.push_back((x>>16)&0xFF);
            v.push_back((x>>8)&0xFF); v.push_back(x&0xFF);
        };
        std::vector<uint8_t> db;
        push_u32(db, (uint32_t)col_drc.size());      push_u32(db, (uint32_t)col_dshift.size());
        push_u32(db, (uint32_t)col_dnsubs.size());   push_u32(db, (uint32_t)col_dsubpos.size());
        push_u32(db, (uint32_t)col_dsubbases.size());push_u32(db, (uint32_t)col_dexlen.size());
        push_u32(db, (uint32_t)col_dexbase.size());  push_u32(db, (uint32_t)col_dN.size());
        for (auto* col : {&col_drc,&col_dshift,&col_dnsubs,&col_dsubpos,
                          &col_dsubbases,&col_dexlen,&col_dexbase,&col_dN})
            db.insert(db.end(), col->begin(), col->end());
        result.delta_bytes = std::move(db);
    }

    // ── 6. Pack bitsets ───────────────────────────────────────────────────────
    result.chain_starts.assign((n + 7) / 8, 0);
    for (int i = 0; i < n; ++i)
        if (chain_starts_vec[i])
            result.chain_starts[i / 8] |= (uint8_t)(1 << (i % 8));

    result.rc_flags.assign((n + 7) / 8, 0);
    for (int i = 0; i < n; ++i)
        if (rc_flags_vec[i])
            result.rc_flags[i / 8] |= (uint8_t)(1 << (i % 8));

    // ── Export chain_order (needed by chain-pg encoder for name ordering) ───────
    result.chain_order = chain_order;
    if (build_pg_) result.has_pg = true;

    // ── 7. Encode quality ────────────────────────────────────────────────────
    // chain_order serves as dfs_order: reads[chain_order[i]] processed at position i.
    MSTSequenceEncoder mst_enc2(cfg_);
    std::vector<int> dummy_shifts(n, 0);
    auto qres = mst_enc2.encode_quality_rans(
        reads, chain_order, parents, dummy_shifts, dev_sets, L);
    result.quality_bytes = std::move(qres.first);
    result.quality_model = std::move(qres.second);

    fprintf(stderr,
        "[CHAIN11] reads=%d chain_starts=%zu (%.1f%%) deltas=%zu (%.1f%%)"
        " raw_delta_stream=%zu B\n",
        n,
        result.n_chain_starts, 100.0 * result.n_chain_starts / n,
        result.n_deltas,       100.0 * result.n_deltas / n,
        result.delta_bytes.size());

    return result;
}

// ── ChainSequenceDecoder::decode_sequences ────────────────────────────────────
// sorted_order is identity {0,1,...,n-1} in trial11 (chain order = output order).
// Logic unchanged from trial10 — works with any sorted_order.
std::vector<std::string> ChainSequenceDecoder::decode_sequences(
    const std::vector<uint8_t>&  vseq_bytes,
    const std::vector<uint8_t>&  delta_bytes,
    const std::vector<uint8_t>&  chain_starts_packed,
    const std::vector<uint32_t>& sorted_order,
    const std::vector<uint8_t>&  rc_flags_packed,
    size_t n_reads,
    int    L,
    std::vector<uint32_t>*           parents_out,
    std::vector<std::vector<bool>>*  dev_sets_out) const {

    (void)rc_flags_packed;  // rc info is in delta blob (col_drc)

    std::vector<std::string> decoded(n_reads);
    if (parents_out)   parents_out->assign(n_reads, UINT32_MAX);
    if (dev_sets_out)  dev_sets_out->assign(n_reads, std::vector<bool>(L, false));

    if (delta_bytes.empty() || sorted_order.size() != n_reads) return decoded;

    // ── Parse CHAIN_VSEQ blob (verbatim seqs + verbatim N-positions) ─────────
    ARCS_CHECK(vseq_bytes.size() >= 8, "CHAIN_VSEQ blob too small");
    auto ru32v = [](const uint8_t* p) -> uint32_t {
        return ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|((uint32_t)p[2]<<8)|p[3];
    };
    uint32_t sz_vseq = ru32v(vseq_bytes.data() + 0);
    uint32_t sz_vN   = ru32v(vseq_bytes.data() + 4);
    const uint8_t* pvseq   = vseq_bytes.data() + 8;
    const uint8_t* pvN     = pvseq + sz_vseq;
    const uint8_t* pvN_end = pvN + sz_vN;

    // ── Parse MST_DELTAS blob (delta fields only) ────────────────────────────
    ARCS_CHECK(delta_bytes.size() >= 32, "MST_DELTAS blob too small for delta header");
    const uint8_t* dhdr = delta_bytes.data();
    uint32_t sz_drc     = ru32v(dhdr +  0);
    uint32_t sz_dshift  = ru32v(dhdr +  4);
    uint32_t sz_dnsubs  = ru32v(dhdr +  8);
    uint32_t sz_dsubpos = ru32v(dhdr + 12);
    uint32_t sz_dsub    = ru32v(dhdr + 16);
    uint32_t sz_dexlen  = ru32v(dhdr + 20);
    uint32_t sz_dexbase = ru32v(dhdr + 24);
    uint32_t sz_dN      = ru32v(dhdr + 28);

    const uint8_t* base = delta_bytes.data() + 32;
    const uint8_t* cdrc     = base;                  base += sz_drc;
    const uint8_t* cdshift  = base; const uint8_t* cdshift_end  = base + sz_dshift;  base += sz_dshift;
    const uint8_t* cdnsubs  = base; const uint8_t* cdnsubs_end  = base + sz_dnsubs;  base += sz_dnsubs;
    const uint8_t* cdsubpos = base; const uint8_t* cdsubpos_end = base + sz_dsubpos; base += sz_dsubpos;
    const uint8_t* cdsub    = base;                  base += sz_dsub;
    const uint8_t* cdexlen  = base; const uint8_t* cdexlen_end  = base + sz_dexlen;  base += sz_dexlen;
    const uint8_t* cdexbase = base;                  base += sz_dexbase;
    const uint8_t* cdN      = base; const uint8_t* cdN_end      = base + sz_dN;      base += sz_dN;
    const uint8_t* pshift   = cdshift;
    const uint8_t* pnsubs   = cdnsubs;
    const uint8_t* psubpos  = cdsubpos;
    const uint8_t* psub     = cdsub;
    const uint8_t* pexlen   = cdexlen;
    const uint8_t* pexbase  = cdexbase;
    const uint8_t* pdN      = cdN;
    int drc_bit_idx = 0;  // next bit in cdrc to read

    auto unpack_col = [&](const uint8_t*& p, int len) -> std::string {
        std::string s(len, 'A');
        int bytes_needed = (len + 3) / 4;
        for (int b = 0; b < bytes_needed; ++b) {
            uint8_t packed = *p++;
            for (int j = 0; j < 4; ++j) {
                int pos = b * 4 + j;
                if (pos >= len) break;
                s[pos] = BASE_TO_CHAR[(packed >> (6 - j * 2)) & 3];
            }
        }
        return s;
    };

    auto restore_N_col = [&](std::string& s, const uint8_t*& p, const uint8_t* pend) {
        int n_N = (int)read_varint(p, pend);
        uint32_t prev = 0;
        for (int i = 0; i < n_N; ++i) {
            uint32_t pos = prev + (uint32_t)read_varint(p, pend);
            prev = pos;
            if (pos < (uint32_t)s.size()) s[pos] = 'N';
        }
    };

    auto read_rc_bit = [&]() -> bool {
        bool val = (cdrc[drc_bit_idx / 8] >> (drc_bit_idx % 8)) & 1;
        ++drc_bit_idx;
        return val;
    };

    std::string enc_frame_prev;

    for (size_t i = 0; i < n_reads; ++i) {
        uint32_t idx = sorted_order[i];
        bool is_chain_start = (chain_starts_packed[i / 8] >> (i % 8)) & 1;

        if (is_chain_start) {
            std::string seq = unpack_col(pvseq, L);
            restore_N_col(seq, pvN, pvN_end);
            decoded[idx]   = seq;
            enc_frame_prev = seq;

        } else {
            bool    rc    = read_rc_bit();
            int32_t shift = zigzag_decode((uint32_t)read_varint(pshift, cdshift_end));
            int     n_subs = (int)read_varint(pnsubs, cdnsubs_end);

            std::vector<uint32_t> sub_pos(n_subs);
            uint32_t prev_p = 0;
            for (int j = 0; j < n_subs; ++j) {
                sub_pos[j] = prev_p + (uint32_t)read_varint(psubpos, cdsubpos_end);
                prev_p = sub_pos[j];
            }

            std::vector<uint8_t> sub_bases(n_subs);
            for (int j = 0; j < n_subs; j += 4) {
                uint8_t packed = *psub++;
                for (int k = 0; k < 4 && j + k < n_subs; ++k)
                    sub_bases[j + k] = (packed >> (6 - k * 2)) & 3;
            }

            int extra_len = (int)read_varint(pexlen, cdexlen_end);
            std::string extra;
            if (extra_len > 0) extra = unpack_col(pexbase, extra_len);

            int pa_start = std::max(0, shift);
            int ch_start = std::max(0, -shift);
            int ov_len   = L - std::abs(shift);

            std::string child(L, 'N');
            for (int j = 0; j < ov_len && (pa_start + j) < (int)enc_frame_prev.size(); ++j)
                child[ch_start + j] = enc_frame_prev[pa_start + j];

            for (int j = 0; j < n_subs; ++j) {
                int pos_in_child = ch_start + (int)sub_pos[j];
                if (pos_in_child >= 0 && pos_in_child < L) {
                    child[pos_in_child] = BASE_TO_CHAR[sub_bases[j]];
                    if (dev_sets_out) {
                        int orig_pos = rc ? (L - 1 - pos_in_child) : pos_in_child;
                        if (orig_pos >= 0 && orig_pos < L)
                            (*dev_sets_out)[idx][orig_pos] = true;
                    }
                }
            }

            if (shift > 0 && extra_len > 0) {
                for (int j = 0; j < extra_len && (ov_len + j) < L; ++j)
                    child[ov_len + j] = extra[j];
            } else if (shift < 0 && extra_len > 0) {
                for (int j = 0; j < extra_len && j < ch_start; ++j)
                    child[j] = extra[j];
            }

            enc_frame_prev = child;

            std::string output_seq = rc ? reverse_complement(child) : child;
            restore_N_col(output_seq, pdN, cdN_end);
            decoded[idx] = output_seq;

            if (parents_out)
                (*parents_out)[idx] = sorted_order[i - 1];
        }
    }

    return decoded;
}

// ── Dedup pseudogenome assembler ─────────────────────────────────────────────
// See chain_encoder.h for the design + losslessness argument.
namespace {

// Strict uppercase-ACGT test. The pg/contigs store ONLY strict A/C/G/T (everything
// else — 'N', lowercase acgt, IUPAC codes, any other byte — is masked to 'A' and its
// literal byte saved in pg_N_char_flat for exact restoration). encode_base() is
// case-insensitive (good for k-mer seeding/mapping) but MUST NOT decide what goes
// into the pg, or lowercase would silently upcase (a losslessness bug).
static inline bool is_acgt_strict(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

constexpr int      DEDUP_K        = 16;    // seed length
int MAX_CAND   = 96;   // candidate positions verified per read (env ARCS_DEDUP_CAND)
int MAX_BUCKET = 200;  // cap index bucket size (env ARCS_DEDUP_BUCKET)

// Record a read that maps onto existing pg at (pos, rc). target = RC^rc(orig).
// Uses the read's ACTUAL length (target.size()) so variable-length reads work.
void record_mapped(ChainEncodeResult& r, uint32_t oi, uint32_t pos, int rc,
                   const std::string& target, int /*L_ignored*/) {
    const int L = (int)target.size();
    uint32_t mm_cnt = 0, N_cnt = 0, qmm_cnt = 0;
    for (int j = 0; j < L; ++j) {
        char tj = target[(size_t)j];
        if (!is_acgt_strict(tj)) {                  // N / lowercase / IUPAC / other
            r.pg_N_pos_flat.push_back((uint16_t)j);
            r.pg_N_char_flat.push_back((uint8_t)tj); // store literal byte (lossless)
            ++N_cnt;
        } else if (r.pg[pos + (uint32_t)j] != tj) { // mismatch vs pg
            r.pg_mm_pos_flat.push_back((uint16_t)j);
            r.pg_mm_base_flat.push_back(encode_base(tj));
            ++mm_cnt;
            uint16_t qpos = (uint16_t)(rc ? (L - 1 - j) : j);  // orig-read frame
            r.pg_qmm_pos_flat.push_back(qpos);
            ++qmm_cnt;
        }
    }
    r.chain_order.push_back(oi);
    r.pg_pos.push_back(pos);
    r.pg_rc.push_back((uint8_t)rc);
    r.pg_readlen.push_back((uint16_t)L);
    r.pg_mm_counts.push_back(mm_cnt);
    r.pg_N_counts.push_back(N_cnt);
    r.pg_qmm_counts.push_back(qmm_cnt);
}

// Append a read as a new forward reference segment (N → 'A' in pg).
void record_append(ChainEncodeResult& r, uint32_t oi, const std::string& orig, int /*L_ignored*/) {
    const int L = (int)orig.size();
    uint32_t pos = (uint32_t)r.pg.size();
    uint32_t N_cnt = 0;
    for (int j = 0; j < L; ++j) {
        char b = orig[(size_t)j];
        if (!is_acgt_strict(b)) {
            r.pg_N_pos_flat.push_back((uint16_t)j);
            r.pg_N_char_flat.push_back((uint8_t)b);  // store literal byte (lossless)
            r.pg.push_back('A');
            ++N_cnt;
        } else {
            r.pg.push_back(b);
        }
    }
    r.chain_order.push_back(oi);
    r.pg_pos.push_back(pos);
    r.pg_rc.push_back(0);
    r.pg_readlen.push_back((uint16_t)L);
    r.pg_mm_counts.push_back(0);
    r.pg_N_counts.push_back(N_cnt);
    r.pg_qmm_counts.push_back(0);
}

} // namespace

ChainEncodeResult build_dedup_pg(const std::vector<Read>& reads) {
    ChainEncodeResult r;
    const int n = (int)reads.size();
    const int L = reads.empty() ? 0 : (int)reads[0].seq.size();
    r.read_len = L;
    r.n_reads  = (size_t)n;
    r.has_pg   = true;
    if (n == 0) return r;
    // NOTE: do NOT bail when reads[0] is short. Per-read lengths are supported; the
    // placement loop bounds-guards seeds and falls through to record_append for any
    // read shorter than K (incl. 0-length), so short/sub-seed reads are handled
    // losslessly as their own contigs instead of crashing the downstream pipeline.

    const int      K     = DEDUP_K;
    const uint64_t KMASK = (K < 32) ? ((1ULL << (2 * K)) - 1) : ~0ULL;
    // Map only near-exact re-observations: appending ~L bases costs ~L·bpb/8
    // compressed, each mismatch ~2 B, so mapping wins only for a few mismatches.
    int    MAXMM  = 4;
    double OV_ERR = 0.10;               // max mismatch rate in an overlap region
    if (const char* s = getenv("ARCS_DEDUP_MAXMM")) MAXMM  = atoi(s);
    if (const char* s = getenv("ARCS_DEDUP_OVERR")) OV_ERR = atof(s);

    // Processing order for the single-pass frontier assembly. File order works
    // best when the FASTQ has genome locality (frontier extension builds long
    // contigs). Greedy-chain order is available for experimentation but tends to
    // fragment into one contig per chain.
    std::vector<uint32_t> order(n);
    std::iota(order.begin(), order.end(), 0u);
    if (const char* s = getenv("ARCS_DEDUP_ORDER")) {
        if (std::string(s) == "chain") {
            ChainSequenceEncoder chain(MSTConfig{}, /*build_pg=*/false);
            auto ch = chain.encode(reads);
            if (ch.chain_order.size() == (size_t)n) order = ch.chain_order;
        }
    }

    std::string& pg = r.pg;
    pg.reserve((size_t)n * 4);

    std::unordered_map<uint64_t, std::vector<uint32_t>> index;
    index.reserve((size_t)n * 2);

    // Pack the K-mer at s[off..off+K); returns false if it contains an N.
    auto pack_kmer = [&](const std::string& s, int off, uint64_t& out) -> bool {
        uint64_t v = 0;
        for (int j = 0; j < K; ++j) {
            uint8_t b = encode_base(s[(size_t)(off + j)]);
            if (b >= 4) return false;
            v = (v << 2) | b;
        }
        out = v & KMASK;
        return true;
    };
    // Index every K-mer of a newly appended pg region [start, start+len).
    auto index_region = [&](uint32_t start, int len) {
        uint64_t v = 0; int valid = 0;
        for (int j = 0; j < len; ++j) {
            uint8_t b = encode_base(pg[start + (uint32_t)j]);   // pg has no N
            v = ((v << 2) | b) & KMASK;
            if (++valid >= K) {
                auto& bucket = index[v];
                if ((int)bucket.size() < MAX_BUCKET)
                    bucket.push_back(start + (uint32_t)j - (uint32_t)(K - 1));
            }
        }
    };

    const int seed_offs[5] = {0, (L - K) / 4, (L - K) / 2, 3 * (L - K) / 4, L - K};
    const int MIN_OV = K;                // min overlap to extend the frontier

    // Index the K-mers ending in [from, to) of the current pg (handles the
    // junction when extending: pass a start a few bases before the append point).
    auto index_range = [&](uint32_t from, uint32_t to) {
        if (to < (uint32_t)K) return;
        uint32_t start = (from >= (uint32_t)(K - 1)) ? from - (uint32_t)(K - 1) : 0;
        uint64_t v = 0; int valid = 0;
        for (uint32_t p = start; p < to; ++p) {
            v = ((v << 2) | encode_base(pg[p])) & KMASK;
            if (++valid >= K) {
                uint32_t kpos = p - (uint32_t)(K - 1);
                if (kpos >= from - std::min(from, (uint32_t)(K - 1))) {
                    auto& bucket = index[v];
                    if ((int)bucket.size() < MAX_BUCKET) bucket.push_back(kpos);
                }
            }
        }
    };
    (void)index_region;  // superseded by index_range

    std::vector<uint32_t> cands;
    for (int oidx = 0; oidx < n; ++oidx) {
        int oi = (int)order[(size_t)oidx];
        const std::string& orig = reads[(size_t)oi].seq;

        // Best internal full match (novel = 0).
        int bestI_pos = -1, bestI_rc = 0, bestI_mm = MAXMM + 1;
        std::string bestI_target;
        // Best frontier extension (append only L-ov novel bases).
        int bestE_pos = -1, bestE_rc = 0, bestE_ov = MIN_OV - 1, bestE_mm = 0;
        std::string bestE_target;

        for (int rc = 0; rc < 2; ++rc) {
            std::string target = rc ? reverse_complement(orig) : orig;
            uint32_t pgsize = (uint32_t)pg.size();

            cands.clear();
            for (int si = 0; si < 5; ++si) {
                int off = seed_offs[si];
                if (off < 0 || off + K > L) continue;   // bounds-guard (short reads)
                uint64_t km;
                if (!pack_kmer(target, off, km)) continue;
                auto it = index.find(km);
                if (it == index.end()) continue;
                for (uint32_t kpos : it->second) {
                    if (kpos < (uint32_t)off) continue;
                    uint32_t cs = kpos - (uint32_t)off;
                    if (cs < pgsize) cands.push_back(cs);   // internal or frontier
                }
            }
            std::sort(cands.begin(), cands.end());
            cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

            int checked = 0;
            for (uint32_t cs : cands) {
                if (checked++ >= MAX_CAND) break;
                if (cs + (uint32_t)L <= pgsize) {
                    // Internal full match.
                    int mm = 0;
                    for (int j = 0; j < L && mm < bestI_mm; ++j) {
                        char tj = target[(size_t)j];
                        if (encode_base(tj) >= 4) continue;
                        if (pg[cs + (uint32_t)j] != tj) ++mm;
                    }
                    if (mm < bestI_mm) {
                        bestI_mm = mm; bestI_pos = (int)cs; bestI_rc = rc; bestI_target = target;
                    }
                } else {
                    // Frontier extension: overlap = pgsize - cs.
                    int ov = (int)(pgsize - cs);
                    if (ov < MIN_OV || ov >= L) continue;
                    int mm = 0;
                    for (int j = 0; j < ov; ++j) {
                        char tj = target[(size_t)j];
                        if (encode_base(tj) >= 4) continue;
                        if (pg[cs + (uint32_t)j] != tj) ++mm;
                    }
                    if (mm <= (int)(ov * OV_ERR) && ov > bestE_ov) {
                        bestE_ov = ov; bestE_pos = (int)cs; bestE_rc = rc;
                        bestE_mm = mm; bestE_target = target;
                    }
                }
            }
        }

        if (bestI_pos >= 0 && bestI_mm <= MAXMM) {
            // Re-observation: map with no new pg.
            record_mapped(r, (uint32_t)oi, (uint32_t)bestI_pos, bestI_rc, bestI_target, L);
            ++r.n_deltas;
        } else if (bestE_pos >= 0) {
            // Frontier extension: append novel suffix target[ov..L], map at bestE_pos.
            uint32_t pos    = (uint32_t)bestE_pos;
            uint32_t oldsize = (uint32_t)pg.size();
            const std::string& t = bestE_target;
            for (int j = bestE_ov; j < L; ++j) {
                char b = t[(size_t)j];
                pg.push_back(encode_base(b) >= 4 ? 'A' : b);
            }
            record_mapped(r, (uint32_t)oi, pos, bestE_rc, t, L);
            index_range(oldsize, (uint32_t)pg.size());
            ++r.n_deltas;
        } else {
            uint32_t pos = (uint32_t)pg.size();
            record_append(r, (uint32_t)oi, orig, L);
            index_range(pos, (uint32_t)pg.size());
            ++r.n_chain_starts;
        }
    }

    fprintf(stderr,
        "[CHAIN-PG-DEDUP] reads=%d mapped=%zu (%.1f%%) appended=%zu (%.1f%%) "
        "pg_len=%zu B (%.2fx smaller than raw reads)\n",
        n, r.n_deltas, 100.0 * (double)r.n_deltas / n,
        r.n_chain_starts, 100.0 * (double)r.n_chain_starts / n,
        pg.size(), (double)((size_t)n * L) / (double)std::max<size_t>(1, pg.size()));

    return r;
}

// ── Flat k-mer → postings index (drop-in for unordered_map<u64,vector<u64>>) ──
// The merge passes hash every contig k-mer into buckets of (contig,pos,strand)
// postings. std::unordered_map<u64,vector<u64>> allocates a vector per distinct
// k-mer (~24M insertions on DS2) — the dominant merge cost. This replaces it with
// open-addressing slots (kmer→head/tail posting) over a single contiguous postings
// pool (append-at-tail preserves INSERTION order, so the greedy sees postings in
// the exact same order → bit-identical merges → archives unchanged). Grows like a
// hash map so it self-sizes without pre-allocating for the 4^k key space.
struct KmerIndex {
    static constexpr uint64_t EMPTY = ~0ULL;    // real k-mers ≤ 2*MK bits (MK≤31) < ~0ULL
    std::vector<uint64_t> hkey;
    std::vector<int32_t>  hhead, htail, hcnt;    // after finalize(): hhead = start into cval
    std::vector<uint64_t> pval;                  // postings pool (build phase)
    std::vector<int32_t>  pnext;
    std::vector<uint64_t> cval;                  // contiguous postings (after finalize)
    size_t mask, count = 0;
    int max_bucket;
    KmerIndex(size_t expect_post, int maxb) : max_bucket(maxb) {
        // Pre-size the key table from the expected posting count so it almost never
        // grow()s during build. grow() rehashes the entire key table on every doubling
        // (~13 doublings from the old fixed 4096-slot start on a multi-M-posting merge),
        // and that repeated rehash was a large fraction of the merge's index-build cost.
        // grow() only remaps the KEY table — it never reorders the postings pool — so
        // starting larger is BIT-IDENTICAL to growing into the same size: finalize()
        // emits the exact same cval order, merges stay byte-for-byte unchanged. Distinct
        // k-mers are usually well under the posting count (genomic repeats), so target
        // ~expect_post/2 slots: enough to hold the keys at <0.7 load without a grow,
        // while also avoiding grow()'s transient old+new double-allocation (small RAM win).
        // Cap the initial size at 2^24 slots (~320 MB of key table) so a huge pg on a
        // whole-genome file can't over-allocate — beyond that it grows a couple of times
        // from 2^24, which is cheap, while bounding peak RAM. Byte output is identical
        // regardless of how many grows occur.
        unsigned bits = 12;
        size_t target = expect_post / 2 + 1024;
        while (((size_t)1 << bits) < target && bits < 24) ++bits;
        init(bits);
        pval.reserve(expect_post + 16); pnext.reserve(expect_post + 16);
    }
    void init(unsigned bits) {
        size_t cap = (size_t)1 << bits;
        hkey.assign(cap, EMPTY); hhead.assign(cap, -1); htail.assign(cap, -1); hcnt.assign(cap, 0);
        mask = cap - 1;
    }
    static inline size_t mix(uint64_t k) { return (size_t)(((k + 1) * 0x9E3779B97F4A7C15ULL) >> 32); }
    inline void insert(uint64_t k, uint64_t v) {
        size_t i = mix(k) & mask;
        while (hkey[i] != EMPTY && hkey[i] != k) i = (i + 1) & mask;
        if (hkey[i] == EMPTY) {
            hkey[i] = k;
            if (++count * 10 >= (mask + 1) * 7) { grow(); i = mix(k) & mask; while (hkey[i] != k) i = (i + 1) & mask; }
        }
        if (hcnt[i] >= max_bucket) return;
        int32_t pi = (int32_t)pval.size(); pval.push_back(v); pnext.push_back(-1);
        if (hhead[i] < 0) hhead[i] = pi; else pnext[htail[i]] = pi;
        htail[i] = pi; ++hcnt[i];
    }
    // Compact each slot's linked postings into one contiguous run (insertion order
    // preserved) so the hot lookup loops read sequential memory, matching the cache
    // behaviour of unordered_map's contiguous vector buckets — critical on datasets
    // with large contigs/buckets (human WGS) where lookups dominate the merge. After
    // this, hhead[slot] = start index into cval, hcnt[slot] = length.
    void finalize() {
        cval.resize(pval.size());
        size_t w = 0;
        for (size_t s = 0; s <= mask; ++s) {
            if (hkey[s] == EMPTY) continue;
            int32_t start = (int32_t)w;
            for (int32_t pi = hhead[s]; pi >= 0; pi = pnext[pi]) cval[w++] = pval[pi];
            hhead[s] = start;
        }
        std::vector<uint64_t>().swap(pval);       // free the build-phase pool
        std::vector<int32_t>().swap(pnext);
    }
    // Contiguous postings range [start, start+out_len) into cval, or out_len 0.
    inline int32_t find_run(uint64_t k, int32_t& out_len) const {
        size_t i = mix(k) & mask;
        while (hkey[i] != EMPTY) {
            if (hkey[i] == k) { out_len = hcnt[i]; return hhead[i]; }
            i = (i + 1) & mask;
        }
        out_len = 0; return -1;
    }
    void grow() {
        std::vector<uint64_t> ok; std::vector<int32_t> oh, ot, oc;
        ok.swap(hkey); oh.swap(hhead); ot.swap(htail); oc.swap(hcnt);
        init((unsigned)__builtin_ctzll(ok.size()) + 1);
        for (size_t s = 0; s < ok.size(); ++s) {
            if (ok[s] == EMPTY) continue;
            size_t i = mix(ok[s]) & mask;
            while (hkey[i] != EMPTY) i = (i + 1) & mask;
            hkey[i] = ok[s]; hhead[i] = oh[s]; htail[i] = ot[s]; hcnt[i] = oc[s];
        }
    }
};

// ── Flat multi-posting placement index (drop-in for unordered_map<u64,vector<u64>>) ─
// Replaces the node-based map for read placement while keeping the SAME information:
// up to MAX_BUCKET positions per k-mer, keep-first (identical to the map's push-back-
// until-full). Storage: open-addressed key table (keys/head/cnt) + one CONTIGUOUS
// postings pool (pval/pnext linked). This kills the map's per-node (~48 B) and per-key
// vector (24 B object + heap) overhead — the dominant peak-RAM structure during
// assembly — WITHOUT dropping any posting, so placement decisions and the archive are
// BYTE-IDENTICAL (probe_seed's candidates are sorted+uniqued, so the linked postings'
// order is irrelevant; only the first-MAX_BUCKET SET matters, which is preserved).
// The compact contiguous layout also lifts memory-bandwidth pressure. Insight behind
// keeping this lean: ARCS placement only needs each k-mer's occurrence list, not a
// general overlap graph. KEY INSIGHT stays: assembly is prediction-like — this is the
// same "k-mer → occurrences" shape as the FCM, now stored without container bloat.
struct FlatPlaceIndex {
    static constexpr uint64_t EMPTY = ~0ULL;
    std::vector<uint64_t> keys;   // open-addressed; EMPTY = free
    std::vector<int32_t>  head;   // slot → first posting index (-1 = none)
    std::vector<int32_t>  cnt;    // postings stored for this slot (for MAX_BUCKET cap)
    std::vector<uint64_t> pval;   // contiguous postings pool (packed cid<<32|pos)
    std::vector<int32_t>  pnext;  // linked list within the pool
    size_t mask = 0, count = 0;
    int max_bucket = 200;
    FlatPlaceIndex() { init(12); }
    void set_cap(int c) { max_bucket = c; }
    void init(unsigned bits) {
        size_t c = (size_t)1 << bits;
        keys.assign(c, EMPTY); head.assign(c, -1); cnt.assign(c, 0); mask = c - 1; count = 0;
    }
    static inline size_t mix(uint64_t k) { return (size_t)(((k + 1) * 0x9E3779B97F4A7C15ULL) >> 32); }
    inline void insert(uint64_t k, uint64_t v) {
        if (count * 10 >= (mask + 1) * 7) grow();
        size_t i = mix(k) & mask;
        while (keys[i] != EMPTY && keys[i] != k) i = (i + 1) & mask;
        if (keys[i] == EMPTY) { keys[i] = k; ++count; }
        if (cnt[i] >= max_bucket) return;                 // keep-first MAX_BUCKET (matches the map)
        int32_t pi = (int32_t)pval.size();
        pval.push_back(v); pnext.push_back(head[i]); head[i] = pi; ++cnt[i];
    }
    // Iterate: for (int32_t pi = idx.find(k); pi >= 0; pi = idx.pnext[pi]) use idx.pval[pi];
    inline int32_t find(uint64_t k) const {
        size_t i = mix(k) & mask;
        while (keys[i] != EMPTY) { if (keys[i] == k) return head[i]; i = (i + 1) & mask; }
        return -1;
    }
    void grow() {
        std::vector<uint64_t> ok; std::vector<int32_t> oh, oc;
        ok.swap(keys); oh.swap(head); oc.swap(cnt);
        init((unsigned)__builtin_ctzll(ok.size()) + 1);
        for (size_t s = 0; s < ok.size(); ++s) {
            if (ok[s] == EMPTY) continue;
            size_t i = mix(ok[s]) & mask;
            while (keys[i] != EMPTY) i = (i + 1) & mask;
            keys[i] = ok[s]; head[i] = oh[s]; cnt[i] = oc[s]; ++count;
        }
    }
};

// ── Multi-contig greedy assembler ────────────────────────────────────────────
// See chain_encoder.h. Two phases: (1) place every read into an independently
// growing contig; (2) concatenate contigs into the pg and emit per-read streams.
ChainEncodeResult build_multicontig_pg(const std::vector<Read>& reads, CallData* call_out) {
    ChainEncodeResult r;
    const int n = (int)reads.size();
    const int L = reads.empty() ? 0 : (int)reads[0].seq.size();
    r.read_len = L;
    r.n_reads  = (size_t)n;
    r.has_pg   = true;
    if (n == 0) return r;
    // NOTE: do NOT bail when reads[0] is short. Per-read lengths are supported; the
    // placement loop bounds-guards seeds and falls through to record_append for any
    // read shorter than K (incl. 0-length), so short/sub-seed reads are handled
    // losslessly as their own contigs instead of crashing the downstream pipeline.

    const int      K     = DEDUP_K;
    const uint64_t KMASK = (K < 32) ? ((1ULL << (2 * K)) - 1) : ~0ULL;
    const int      MIN_OV = K;
    int    MAXMM  = 4;
    double OV_ERR = 0.06;
    if (const char* s = getenv("ARCS_DEDUP_MAXMM"))  MAXMM      = atoi(s);
    if (const char* s = getenv("ARCS_DEDUP_OVERR"))  OV_ERR     = atof(s);
    if (const char* s = getenv("ARCS_DEDUP_CAND"))   MAX_CAND   = atoi(s);
    if (const char* s = getenv("ARCS_DEDUP_BUCKET")) MAX_BUCKET = atoi(s);

    // ── Sparse placement index — OPT-IN (ARCS_PLACE_SPARSE=w) ─────────────────────
    // On errored/high-coverage data the placement index balloons with spurious distinct
    // k-mers (error k-mers are singletons that never match anything — pure RAM waste).
    // Indexing only the ~1/w of k-mers whose hash is sampled (and probing with the SAME
    // predicate) shrinks the index ~w× — the located bacterial-RAM residual. Placement is
    // lossless by construction (it only chooses where reads go; the codec stores reads+
    // mismatches vs the resulting contigs), so this is a pure RATIO trade, never a
    // losslessness risk — VALIDATED byte-exact on clean-ecoli / errored-ecoli / human,
    // ratio neutral-to-better. place_smask==0 (default) = index/query EVERY position =
    // BIT-IDENTICAL to the original placement. Best with w a power of two. Stacks with
    // sparse-merge for −52% assembly RAM on the errored bacterial proxy. See §8d.
    uint64_t place_smask = 0;
    if (const char* s = getenv("ARCS_PLACE_SPARSE")) { int v = atoi(s); if (v > 1) place_smask = (uint64_t)v - 1; }

    // High-quality read pre-filter (PgRC2 principle): reads with N symbols or
    // low-quality tails are excluded from contig extension and seeding. They are
    // still placed on existing contigs and encoded losslessly — only the assembly
    // growth step is gated. On high-coverage data (e.g. 30× E.coli) this removes
    // error k-mers from the dBG → fewer dead ends → more contiguous pg → better ratio.
    // Disable with ARCS_NO_HQ_FILTER=1.
    const bool hq_filter = (getenv("ARCS_NO_HQ_FILTER") == nullptr);
    std::vector<bool> is_hq((size_t)n, true);
    if (hq_filter) {
        for (int i = 0; i < n; ++i) {
            const std::string& s = reads[(size_t)i].seq;
            const std::string& q = reads[(size_t)i].qual;
            const int rl = (int)s.size();
            bool hq = true;
            for (char c : s) { if (c == 'N' || c == 'n') { hq = false; break; } }
            if (hq && !q.empty()) {
                int tail = std::max(1, (int)(0.12 * rl));
                for (int j = rl - tail; j < rl; ++j)
                    if ((q[(size_t)j] - 33) <= 2) { hq = false; break; }
            }
            is_hq[(size_t)i] = hq;
        }
    }

    std::vector<std::string> contigs;
    contigs.reserve((size_t)n / 4 + 1);

    // Global k-mer index (serial path): kmer → up to MAX_BUCKET contig positions,
    // flat/contiguous (byte-identical to the old map, less RAM — see FlatPlaceIndex).
    FlatPlaceIndex index; index.set_cap(MAX_BUCKET);

    // Per-read placement (in file order → chain_order is identity).
    std::vector<uint32_t> pl_cid(n), pl_pos(n);
    std::vector<uint8_t>  pl_rc(n, 0);
    std::vector<uint8_t>  pl_new(n, 0);   // 1 if this read started a new contig

    auto pack_kmer = [&](const std::string& s, int off, uint64_t& out) -> bool {
        uint64_t v = 0;
        for (int j = 0; j < K; ++j) {
            uint8_t b = encode_base(s[(size_t)(off + j)]);
            if (b >= 4) return false;
            v = (v << 2) | b;
        }
        out = v & KMASK;
        return true;
    };

    std::vector<uint64_t> cands;   // (cid<<32|cs) candidate read-starts

    // Minimizer anchors: content-based seed positions that coincide between two
    // overlapping strings, so a read shares an anchor with its true neighbour even
    // when the overlap is small or shifted (where the 5 fixed seeds miss). We probe
    // these IN ADDITION to the fixed seeds (a strict superset of candidates), so
    // placement can only rise — never regress. Disable with ARCS_MINZ_OFF.
    const bool use_minz = (getenv("ARCS_MINZ_ON") != nullptr);  // OFF by default: seed-sensitivity is not the bottleneck (measured)
    int MINZ_W = 8;
    if (const char* s = getenv("ARCS_MINZ_W")) { int v = atoi(s); if (v>=2 && v<=64) MINZ_W = v; }
    auto hash64 = [](uint64_t x) -> uint64_t {
        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
        x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x;
    };
    std::vector<uint64_t> kh_buf;
    std::vector<int> minz_off;
    auto compute_minimizers = [&](const std::string& s, int Ls, std::vector<int>& out) {
        out.clear();
        int nk = Ls - K + 1;
        if (nk <= 0) return;
        kh_buf.assign((size_t)nk, ~0ULL);
        for (int p = 0; p + K <= Ls; ++p) {
            uint64_t km;
            if (pack_kmer(s, p, km)) kh_buf[(size_t)p] = hash64(km);
        }
        int last = -1;
        int W = (MINZ_W < nk) ? MINZ_W : nk;
        for (int i = 0; i + W <= nk; ++i) {
            int mpos = i; uint64_t mh = kh_buf[(size_t)i];
            for (int j = i + 1; j < i + W; ++j)
                if (kh_buf[(size_t)j] < mh) { mh = kh_buf[(size_t)j]; mpos = j; }
            if (mh != ~0ULL && mpos != last) { out.push_back(mpos); last = mpos; }
        }
    };

    const bool ASM_TIMING = getenv("ARCS_ASM_TIMING") != nullptr;
    auto _tnow = [] { return std::chrono::steady_clock::now(); };
    auto _tp0 = _tnow();

    // ── Phase 1: placement ────────────────────────────────────────────────────
    // ── PROTOTYPE: shard-parallel assembly test (ARCS_SHARD_TEST=T) ─────────────
    // Clearing the k-mer index at shard boundaries makes reads match ONLY contigs
    // from their own shard → simulates independent per-shard assembly (SERIALLY, so
    // ratio == the parallel version). The existing end-of-assembly MERGE then
    // reconciles the union of all shards' contigs — this is the "merge as
    // parallelization primitive" hypothesis. If final pg ≈ global pg, the idea holds.
    // place_range: run Phase-1 placement over reads [begin,end) into an INDEPENDENT
    // (contigs, index) pair with its own scratch. Called once over all reads (serial),
    // or once per contiguous read-shard on separate threads (shard-parallel assembly).
    // Writes pl_cid (LOCAL to the given contig set), pl_pos, pl_rc, pl_new for its read
    // range only. Threads use disjoint read ranges + private contigs/index → race-free.
    // (Local index_contig / minimizer scratch shadow the outer ones so the outer copies
    // are unused for the serial path but kept for source stability.)
    auto place_range = [&](const std::vector<int>& ids,
                           std::vector<std::string>& contigs,
                           FlatPlaceIndex& index) {
        std::vector<uint64_t> cands;
        std::vector<uint64_t> kh_buf;
        std::vector<int> minz_off;
        auto index_contig = [&](uint32_t cid, uint32_t from, uint32_t to) {
            const std::string& c = contigs[cid];
            if (to < (uint32_t)K) return;
            uint32_t start = (from >= (uint32_t)(K - 1)) ? from - (uint32_t)(K - 1) : 0;
            uint64_t v = 0; int valid = 0;
            for (uint32_t p = start; p < to; ++p) {
                v = ((v << 2) | encode_base(c[p])) & KMASK;
                if (++valid >= K && (place_smask == 0 || (hash64(v) & place_smask) == 0)) {
                    uint32_t kpos = p - (uint32_t)(K - 1);
                    index.insert(v, ((uint64_t)cid << 32) | kpos);   // keep-first, single slot
                }
            }
        };
        auto compute_minimizers = [&](const std::string& s, int Ls, std::vector<int>& out) {
            out.clear();
            int nk = Ls - K + 1;
            if (nk <= 0) return;
            kh_buf.assign((size_t)nk, ~0ULL);
            for (int p = 0; p + K <= Ls; ++p) {
                uint64_t km;
                if (pack_kmer(s, p, km)) kh_buf[(size_t)p] = hash64(km);
            }
            int last = -1;
            int W = (MINZ_W < nk) ? MINZ_W : nk;
            for (int i = 0; i + W <= nk; ++i) {
                int mpos = i; uint64_t mh = kh_buf[(size_t)i];
                for (int j = i + 1; j < i + W; ++j)
                    if (kh_buf[(size_t)j] < mh) { mh = kh_buf[(size_t)j]; mpos = j; }
                if (mh != ~0ULL && mpos != last) { out.push_back(mpos); last = mpos; }
            }
        };
      for (int oi : ids) {
        const std::string& orig = reads[(size_t)oi].seq;
        const int L = (int)orig.size();          // per-read length (variable-length safe)
        // Reads shorter than K cannot seed -> fall through to new-contig append
        // (record_append handles any length). Seed offsets are bounds-guarded below.
        const int seed_offs[5] = {0, (L - K) / 4, (L - K) / 2, 3 * (L - K) / 4, L - K};

        int bestI_cid = -1, bestI_pos = 0, bestI_rc = 0, bestI_mm = MAXMM + 1;
        int bestE_cid = -1, bestE_pos = 0, bestE_rc = 0, bestE_ov = MIN_OV - 1;
        std::string bestE_target;

        for (int rc = 0; rc < 2; ++rc) {
            std::string target = rc ? reverse_complement(orig) : orig;

            cands.clear();
            auto probe_seed = [&](int off) {
                if (off < 0 || off + K > L) return;     // bounds-guard (short reads)
                uint64_t km;
                if (!pack_kmer(target, off, km)) return;
                for (int32_t pi = index.find(km); pi >= 0; pi = index.pnext[pi]) {
                    uint64_t packed = index.pval[pi];
                    uint32_t cid  = (uint32_t)(packed >> 32);
                    uint32_t kpos = (uint32_t)(packed & 0xFFFFFFFFu);
                    if (kpos < (uint32_t)off) continue;
                    uint32_t cs = kpos - (uint32_t)off;
                    cands.push_back(((uint64_t)cid << 32) | cs);
                }
            };
            if (place_smask == 0) {
                for (int si = 0; si < 5; ++si) probe_seed(seed_offs[si]);   // fixed anchors
                if (use_minz) {                                            // + minimizer anchors
                    compute_minimizers(target, L, minz_off);
                    for (int off : minz_off) probe_seed(off);
                }
            } else {
                // Sparse mode: only sampled k-mers are indexed, so probe every position
                // whose k-mer is sampled (same predicate as the insert). A real overlap
                // ≥ K contains sampled seeds w.h.p.; acceptance is still verified below.
                for (int off = 0; off + K <= L; ++off) {
                    uint64_t km;
                    if (pack_kmer(target, off, km) && (hash64(km) & place_smask) == 0) probe_seed(off);
                }
            }
            std::sort(cands.begin(), cands.end());
            cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

            int checked = 0;
            for (uint64_t packed : cands) {
                if (checked++ >= MAX_CAND) break;
                uint32_t cid = (uint32_t)(packed >> 32);
                uint32_t cs  = (uint32_t)(packed & 0xFFFFFFFFu);
                uint32_t clen = (uint32_t)contigs[cid].size();
                const std::string& c = contigs[cid];
                if (cs + (uint32_t)L <= clen) {
                    int mm = 0;
                    for (int j = 0; j < L && mm < bestI_mm; ++j) {
                        char tj = target[(size_t)j];
                        if (encode_base(tj) >= 4) continue;
                        if (c[cs + (uint32_t)j] != tj) ++mm;
                    }
                    if (mm < bestI_mm) {
                        bestI_mm = mm; bestI_cid = (int)cid; bestI_pos = (int)cs; bestI_rc = rc;
                    }
                } else if (cs < clen) {
                    int ov = (int)(clen - cs);
                    if (ov < MIN_OV || ov >= L) continue;
                    int mm = 0;
                    for (int j = 0; j < ov; ++j) {
                        char tj = target[(size_t)j];
                        if (encode_base(tj) >= 4) continue;
                        if (c[cs + (uint32_t)j] != tj) ++mm;
                    }
                    if (mm <= (int)(ov * OV_ERR) && ov > bestE_ov) {
                        bestE_ov = ov; bestE_cid = (int)cid; bestE_pos = (int)cs;
                        bestE_rc = rc; bestE_target = target;
                    }
                }
            }
        }

        const bool read_is_hq = is_hq[(size_t)oi];
        if (bestI_cid >= 0 && bestI_mm <= MAXMM) {
            // Internal placement — always allowed for both HQ and LQ reads.
            pl_cid[oi] = (uint32_t)bestI_cid; pl_pos[oi] = (uint32_t)bestI_pos;
            pl_rc[oi]  = (uint8_t)bestI_rc;
        } else if (read_is_hq && bestE_cid >= 0) {
            // Contig extension — HQ reads only (LQ reads carry error k-mers that
            // would create spurious branches and fragment the assembly).
            uint32_t cid = (uint32_t)bestE_cid;
            uint32_t oldlen = (uint32_t)contigs[cid].size();
            const std::string& t = bestE_target;
            for (int j = bestE_ov; j < L; ++j) {
                char b = t[(size_t)j];
                contigs[cid].push_back(is_acgt_strict(b) ? b : 'A');
            }
            pl_cid[oi] = cid; pl_pos[oi] = (uint32_t)bestE_pos; pl_rc[oi] = (uint8_t)bestE_rc;
            index_contig(cid, oldlen, (uint32_t)contigs[cid].size());
        } else if (read_is_hq) {
            // New contig — HQ reads only, added to index so future reads can place on it.
            uint32_t cid = (uint32_t)contigs.size();
            std::string c((size_t)L, 'A');
            for (int j = 0; j < L; ++j) {
                char b = orig[(size_t)j];
                if (is_acgt_strict(b)) c[(size_t)j] = b;
            }
            contigs.push_back(std::move(c));
            pl_cid[oi] = cid; pl_pos[oi] = 0; pl_rc[oi] = 0; pl_new[oi] = 1;
            index_contig(cid, 0, (uint32_t)L);
        } else {
            // LQ read: no valid internal placement found. Store as singleton contig
            // but do NOT add to the placement index — prevents error k-mers from
            // attracting future reads and fragmenting the assembly.
            uint32_t cid = (uint32_t)contigs.size();
            std::string c((size_t)L, 'A');
            for (int j = 0; j < L; ++j) {
                char b = orig[(size_t)j];
                if (is_acgt_strict(b)) c[(size_t)j] = b;
            }
            contigs.push_back(std::move(c));
            pl_cid[oi] = cid; pl_pos[oi] = 0; pl_rc[oi] = 0; pl_new[oi] = 1;
            // intentionally no index_contig — LQ singletons are opaque to placement
        }
      }
    };

    // Dispatch: serial (ARCS_SHARD<=1) or shard-parallel assembly (ARCS_SHARD=T).
    // Each shard assembles a contiguous read range into its own (contigs,index) on a
    // thread; shards are concatenated (deterministic order + cid offsets) and the
    // Phase-1.6 MERGE below reconciles cross-shard overlaps. Deterministic and
    // ratio-preserving (validated ±0.3%). n<4*T falls back to serial (too few reads).
    // Default: auto-enable shard-parallel assembly for large inputs (it strictly
    // dominates serial — same ratio ±0.3 %, faster placement, lower RAM). Small inputs
    // stay serial (placement is a negligible fraction). ARCS_SHARD=N overrides;
    // ARCS_SHARD=1 or ARCS_NO_SHARD forces serial.
    int SHARD = 1;
    bool minz_route = (getenv("ARCS_SHARD_MINZ") != nullptr);   // route by content (minimizer)
    if (const char* e = getenv("ARCS_SHARD_MINZ")) { int v = atoi(e); if (v >= 2 && v <= 256) SHARD = v; }
    else if (const char* e2 = getenv("ARCS_SHARD")) { int v = atoi(e2); if (v >= 1 && v <= 256) SHARD = v; }
    else if (getenv("ARCS_NO_SHARD") == nullptr && n >= 200000) {
        unsigned hc = std::thread::hardware_concurrency(); if (!hc) hc = 4;
        SHARD = (int)std::min<unsigned>(hc, 8);
    }
    if (SHARD <= 1 || (size_t)n < 4 * (size_t)SHARD) {
        std::vector<int> ids((size_t)n); std::iota(ids.begin(), ids.end(), 0);
        // HQ reads process first: they build clean assembly; LQ reads find placements
        // on already-built contigs rather than seeding fragmented singleton contigs.
        std::stable_partition(ids.begin(), ids.end(), [&](int i){ return is_hq[(size_t)i]; });
        place_range(ids, contigs, index);
    } else {
        const int T = SHARD;
        // Build per-shard read-id lists (ascending oi → deterministic placement order).
        //  • minz_route: shard = f(read's MIN k-mer hash) → reads sharing that k-mer (≈
        //    same genomic locus) co-locate, so shards are genomically DISJOINT and the
        //    downstream MERGE has little cross-shard redundancy to reconcile.
        //  • else: contiguous file-order ranges (each shard is a whole-genome sample →
        //    merge grows ~T×; kept for comparison).
        std::vector<std::vector<int>> shard_ids((size_t)T);
        if (minz_route) {
            for (int oi = 0; oi < n; ++oi) {
                const std::string& s = reads[(size_t)oi].seq;
                const int Ls = (int)s.size();
                uint64_t mn = ~0ULL;
                for (int p = 0; p + K <= Ls; ++p) { uint64_t km; if (pack_kmer(s, p, km)) { uint64_t h = hash64(km); if (h < mn) mn = h; } }
                int sh = (mn == ~0ULL) ? (oi % T) : (int)((mn >> 13) % (uint64_t)T);
                shard_ids[(size_t)sh].push_back(oi);
            }
        } else {
            const size_t span = ((size_t)n + (size_t)T - 1) / (size_t)T;
            for (int oi = 0; oi < n; ++oi) {
                size_t sh = (size_t)oi / span; if (sh >= (size_t)T) sh = (size_t)T - 1;
                shard_ids[sh].push_back(oi);
            }
        }
        // Within each shard: HQ reads first so they build clean assembly before LQ reads place.
        for (auto& sh : shard_ids)
            std::stable_partition(sh.begin(), sh.end(), [&](int i){ return is_hq[(size_t)i]; });
        std::vector<std::vector<std::string>> sc((size_t)T);
        std::vector<FlatPlaceIndex> si((size_t)T);
        for (auto& x : si) x.set_cap(MAX_BUCKET);
        std::vector<std::thread> th; th.reserve((size_t)T);
        for (int s = 0; s < T; ++s) {
            if (shard_ids[(size_t)s].empty()) continue;
            th.emplace_back([&, s] { place_range(shard_ids[(size_t)s], sc[(size_t)s], si[(size_t)s]); });
        }
        for (auto& t : th) t.join();
        // Concatenate shard contigs deterministically; offset each shard's pl_cid.
        uint32_t off = 0;
        for (int s = 0; s < T; ++s) {
            for (int oi : shard_ids[(size_t)s]) pl_cid[oi] += off;
            for (auto& c : sc[(size_t)s]) contigs.push_back(std::move(c));
            off += (uint32_t)sc[(size_t)s].size();
        }
        std::vector<FlatPlaceIndex>().swap(si);
    }

    if (ASM_TIMING) { auto t=_tnow(); fprintf(stderr,"[ASM] Phase1 placement: %.2fs\n", std::chrono::duration<double>(t-_tp0).count()); _tp0=t; }

    // ── Phase 1.6: contig MERGE (overlap-layout) — collapse redundant contigs ──
    // Phase 1 only extends contig ends, so one genome locus spawns many partially
    // overlapping contigs (≈5× redundant pg). Here we link contigs whose SUFFIX
    // overlaps another contig's PREFIX (exact-ish), walk the links into long
    // super-contigs, and remap every read's (cid,pos). Shrinks the pg toward genome
    // size → smaller sequence blob AND unique read placement for self-alignment.
    // Lossless: reads still record mismatches vs the (re-polished) super-contig.
    // Gated by ARCS_MERGE_CONTIGS. Passes controllable via ARCS_MERGE_PASSES.
    // DEFAULT ON: the merge collapses the redundant contigs Phase-1's online greedy
    // produced (fixing the ordering-fragmentation that costs sequence ratio). 1 pass
    // captures ~96% of the gain and already beats SPRING's sequence; further passes
    // add <1% for 3× the time. Disable with ARCS_MERGE_OFF; tune with ARCS_MERGE_PASSES.
    if (getenv("ARCS_MERGE_OFF") == nullptr) {
        int MERGE_PASSES = 1;
        if (const char* s = getenv("ARCS_MERGE_PASSES")) MERGE_PASSES = atoi(s);
        // Merge SEED size — smaller than K so a seed survives the read error rate
        // (P(clean k-mer)=(1-e)^k; at e=3% a 16-mer is clean only 61% of the time,
        // an 11-mer 72%). Overlaps are still VERIFIED full-length, so a small seed
        // just finds more candidate overlaps without loosening acceptance.
        int MK = 11;
        if (const char* s = getenv("ARCS_MERGE_K")) MK = atoi(s);
        if (MK < 6) MK = 6; if (MK > 31) MK = 31;
        const uint64_t MKMASK = (MK < 32) ? ((1ULL << (2 * MK)) - 1) : ~0ULL;
        auto mpack = [&](const std::string& s, int off, uint64_t& out) -> bool {
            uint64_t v = 0;
            for (int j = 0; j < MK; ++j) { uint8_t b = encode_base(s[(size_t)(off + j)]); if (b >= 4) return false; v = (v << 2) | b; }
            out = v & MKMASK; return true;
        };
        const int MIN_OV_C = K;                 // min contig-contig overlap (verified)
        // ── Sparse (sampled) merge index — OPT-IN (ARCS_MERGE_SPARSE=w) ────────────
        // The ext-grow index-build inserts every k-mer of every contig on both strands
        // (millions of serial inserts) → the merge's dominant cost. Sampling only the
        // k-mers whose hash lands in 1/w of the space cuts inserts ~w× (build + RAM),
        // and the QUERY samples with the SAME predicate so a real overlap ≥ K still
        // contains sampled seeds w.h.p. Overlaps are still VERIFIED full-length, so this
        // never loosens acceptance — it can only miss a few marginal overlaps (a pure
        // RATIO trade, never a losslessness risk: the codec stores reads+mismatches vs
        // whatever contigs result). smask==0 (default) = insert/query EVERY position =
        // BIT-IDENTICAL to the original merge. Best with w a power of two.
        uint64_t smask = 0;
        if (const char* s = getenv("ARCS_MERGE_SPARSE")) { int v = atoi(s); if (v > 1) smask = (uint64_t)v - 1; }
        auto sampled = [&](uint64_t kk) -> bool { return smask == 0 || (KmerIndex::mix(kk) & smask) == 0; };
        const bool MG_TIMING = getenv("ARCS_MERGE_TIMING") != nullptr;
        auto _mn = []{ return std::chrono::steady_clock::now(); };
        auto _mp = _mn();
        auto _mmark = [&](const char* w){ if(MG_TIMING){ auto t=_mn(); fprintf(stderr,"[MERGE-T] %s: %.2fs\n", w, std::chrono::duration<double>(t-_mp).count()); _mp=t; } };
        // ── RC-aware iterative greedy super-contig growth ─────────────────────
        // Grow each super-contig by re-seeding its ACTUAL suffix each step: the
        // real string carries orientation, so RC extensions compose correctly
        // (no orientation-propagation bookkeeping). Index both strands of every
        // contig; on an RC extension, remap that contig's reads with the flip
        // transform (pos → off + lenD − pos − L, rc → 1−rc). Lossless.
        for (int pass = 0; pass < MERGE_PASSES; ++pass) {
            const int MC = (int)contigs.size();
            if (MC < 2) break;
            std::vector<std::string> rcstr(MC);
            // est_post drives the KmerIndex table pre-size. With sparse sampling only
            // ~1/(smask+1) of k-mers are actually inserted, so sizing for the FULL count
            // over-allocates the hash table by that factor — the dominant merge RAM on
            // fragmented (bacterial) data. Size for the sampled count instead. Output is
            // byte-identical (table capacity never affects posting order).
            size_t est_post = 0; for (int c = 0; c < MC; ++c) est_post += contigs[c].size() * 2;
            if (smask) est_post = est_post / (smask + 1) + 1024;
            KmerIndex cidx(est_post, MAX_BUCKET);      // kmer → (cid<<33|pos<<1|strand)
            // RC strings are independent per contig → compute them in parallel (the
            // single biggest serial cost of the merge index build).
            {
                unsigned hc = std::thread::hardware_concurrency(); if (!hc) hc = 4;
                int nt = (int)std::min<size_t>(hc, std::max(1, MC / 2000));
                if (nt <= 1) {
                    for (int c = 0; c < MC; ++c) rcstr[c] = reverse_complement(contigs[c]);
                } else {
                    std::vector<std::thread> th; th.reserve((size_t)nt);
                    for (int t = 0; t < nt; ++t)
                        th.emplace_back([&, t]{
                            for (int c = t; c < MC; c += nt) rcstr[c] = reverse_complement(contigs[c]);
                        });
                    for (auto& x : th) x.join();
                }
            }
            for (int c = 0; c < MC; ++c) {
                for (int st = 0; st < 2; ++st) {
                    const std::string& s = st ? rcstr[c] : contigs[c];
                    if ((int)s.size() < MK) continue;
                    uint64_t v = 0; int val = 0;
                    for (uint32_t p = 0; p < s.size(); ++p) {
                        v = ((v << 2) | encode_base(s[p])) & MKMASK;
                        if (++val >= MK && sampled(v))
                            cidx.insert(v, ((uint64_t)c << 33) | ((uint64_t)(p - (uint32_t)MK + 1) << 1) | (uint64_t)st);
                    }
                }
            }
            cidx.finalize();
            _mmark("ext-grow index build");
            std::vector<char> used(MC, 0);
            std::vector<int> o2n(MC, -1); std::vector<uint32_t> o2off(MC, 0); std::vector<uint8_t> o2rc(MC, 0);
            std::vector<int> ord(MC); std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(), [&](int a, int b){
                return contigs[a].size() != contigs[b].size() ? contigs[a].size() > contigs[b].size() : a < b; });
            std::vector<std::string> supers;
            for (int seed : ord) {
                if (used[seed]) continue;
                int nid = (int)supers.size();
                std::string sup = contigs[seed];
                used[seed] = 1; o2n[seed] = nid; o2off[seed] = 0; o2rc[seed] = 0;
                bool grew = true; int guard = 0;
                while (grew && guard++ < MC + 4) {
                    grew = false;
                    int bestD = -1; uint32_t bestOv = 0; int bestSt = 0;
                    // Look up one seed anchored at sup position `ap` and update the best
                    // overlap. With smask==0 the sampled() guard is inert, so this is
                    // exactly the original inner body → default merge stays bit-identical.
                    auto probe = [&](int ap) {
                        if (ap < 0) return;
                        uint64_t km; if (!mpack(sup, ap, km)) return;
                        if (smask && (KmerIndex::mix(km) & smask) != 0) return;   // only sampled seeds are indexed
                        int32_t _len; int32_t _st = cidx.find_run(km, _len);
                        for (int32_t _q = 0; _q < _len; ++_q) {
                            uint64_t pk = cidx.cval[_st + _q];
                            int d = (int)(pk >> 33); if (used[d]) continue;
                            uint32_t dp = (uint32_t)((pk >> 1) & 0x7FFFFFFFu); int st = (int)(pk & 1);
                            const std::string& Ds = st ? rcstr[d] : contigs[d];
                            uint32_t ov = dp + ((uint32_t)sup.size() - (uint32_t)ap);
                            if (ov < (uint32_t)MIN_OV_C || ov > sup.size() || ov >= Ds.size()) continue;
                            uint32_t astart = (uint32_t)sup.size() - ov;
                            int mm = 0, lim = (int)(ov * OV_ERR) + 1;
                            for (uint32_t j = 0; j < ov && mm < lim; ++j) if (sup[astart + j] != Ds[j]) ++mm;
                            if (mm < lim && ov > bestOv) { bestOv = ov; bestD = d; bestSt = st; }
                        }
                    };
                    if (smask == 0) {
                        // ORIGINAL 7-probe schedule — preserved exactly for the default path.
                        for (int back = 0; back <= 6; ++back) {
                            int ap = (int)sup.size() - MK - back * (MK / 2);
                            if (ap < 0) break;
                            probe(ap);
                        }
                    } else {
                        // Sparse mode: scan a suffix window densely and let probe() skip
                        // the non-sampled positions. The window spans ~2 sampling periods
                        // worth of k-mers so it reliably contains sampled seeds of any
                        // overlapping neighbour.
                        int top = (int)sup.size() - MK;
                        int scanlen = (int)(smask + 1) * MK * 2;
                        for (int ap = top; ap >= 0 && ap >= top - scanlen; --ap) probe(ap);
                    }
                    if (bestD >= 0) {
                        const std::string& Ds = bestSt ? rcstr[bestD] : contigs[bestD];
                        uint32_t S0 = (uint32_t)sup.size();
                        sup += Ds.substr(bestOv);
                        used[bestD] = 1; o2n[bestD] = nid; o2off[bestD] = S0 - bestOv; o2rc[bestD] = (uint8_t)bestSt;
                        grew = true;
                    }
                }
                supers.push_back(std::move(sup));
            }
            for (int c = 0; c < MC; ++c) if (!used[c]) { o2n[c] = (int)supers.size(); o2off[c] = 0; o2rc[c] = 0; supers.push_back(contigs[c]); used[c] = 1; }
            for (int oi = 0; oi < n; ++oi) {
                int c = (int)pl_cid[oi]; uint32_t p = pl_pos[oi]; uint8_t rc = pl_rc[oi];
                uint32_t rlen = (uint32_t)reads[(size_t)oi].seq.size();   // per-read length
                uint32_t np; uint8_t nrc;
                if (o2rc[c]) { uint32_t base = o2off[c] + (uint32_t)contigs[c].size(); uint32_t sub = p + rlen;
                               np = (base >= sub) ? (base - sub) : 0; nrc = (uint8_t)(rc ^ 1); }
                else         { np = o2off[c] + p; nrc = rc; }
                pl_cid[oi] = (uint32_t)o2n[c]; pl_pos[oi] = np; pl_rc[oi] = nrc;
            }
            _mmark("ext-grow greedy+remap");
            size_t before = contigs.size(); contigs.swap(supers);
            fprintf(stderr, "[MERGE] ext-grow pass %d: %zu → %zu contigs\n", pass, before, contigs.size());
            if (contigs.size() == before) break;
        }

        // ── Containment absorption: drop contigs fully inside another ──────────
        for (int cpass = 0; cpass < MERGE_PASSES; ++cpass) {
            const int MC = (int)contigs.size();
            if (MC < 2) break;
            size_t est_post = 0; for (int c = 0; c < MC; ++c) est_post += contigs[c].size();
            KmerIndex cidx(est_post, MAX_BUCKET);
            for (int c = 0; c < MC; ++c) {
                const std::string& s = contigs[c];
                if ((int)s.size() < MK) continue;
                uint64_t v = 0; int val = 0;
                for (uint32_t p = 0; p < s.size(); ++p) {
                    v = ((v << 2) | encode_base(s[p])) & MKMASK;
                    if (++val >= MK) cidx.insert(v, ((uint64_t)c << 32) | (p - (uint32_t)MK + 1));
                }
            }
            cidx.finalize();
            std::vector<int> ord(MC); std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(), [&](int x, int y){
                return contigs[x].size() != contigs[y].size() ? contigs[x].size() < contigs[y].size() : x < y; });
            std::vector<int> host(MC, -1); std::vector<uint32_t> hoff(MC, 0);
            std::vector<uint32_t> hlen(MC, 0); std::vector<uint8_t> hrc(MC, 0);
            std::vector<char> absorbed(MC, 0);
            // Try to place `probe` (B forward or RC(B)) fully inside some larger A.
            auto try_fit = [&](int b, const std::string& probe, int& bestA, uint32_t& bestOff, int& bestMm){
                for (int si = 0; si <= 8; ++si) {
                    int off = si * ((int)probe.size() - MK) / 8;
                    if (off < 0 || off + MK > (int)probe.size()) continue;
                    uint64_t km; if (!mpack(probe, off, km)) continue;
                    int32_t _len; int32_t _st = cidx.find_run(km, _len);
                    for (int32_t _q = 0; _q < _len; ++_q) {
                        uint64_t pk = cidx.cval[_st + _q];
                        int a = (int)(pk >> 32); uint32_t ap = (uint32_t)(pk & 0xFFFFFFFFu);
                        if (a == b || absorbed[a] || ap < (uint32_t)off) continue;
                        // A must be strictly bigger, or equal-size with smaller id (deterministic ties).
                        if (contigs[a].size() < probe.size()) continue;
                        if (contigs[a].size() == probe.size() && a >= b) continue;
                        uint32_t boff = ap - (uint32_t)off;
                        if (boff + probe.size() > contigs[a].size()) continue;
                        int mm = 0, lim = (int)(probe.size() * OV_ERR) + 1;
                        for (uint32_t j = 0; j < probe.size() && mm < lim; ++j)
                            if (contigs[a][boff + j] != probe[j]) ++mm;
                        if (mm < lim && mm < bestMm) { bestMm = mm; bestA = a; bestOff = boff; }
                    }
                }
            };
            for (int b : ord) {
                const std::string& B = contigs[b];
                if ((int)B.size() < MK) continue;
                int bestA = -1; uint32_t bestOff = 0; int bestMm = 1 << 30, bestRc = 0;
                try_fit(b, B, bestA, bestOff, bestMm);                 // forward
                int fA = bestA; uint32_t fOff = bestOff; int fMm = bestMm;
                int rA = -1; uint32_t rOff = 0; int rMm = 1 << 30;
                std::string rcB = reverse_complement(B);              // RC-aware
                try_fit(b, rcB, rA, rOff, rMm);
                if (fA >= 0 && (rA < 0 || fMm <= rMm)) { bestA = fA; bestOff = fOff; bestRc = 0; }
                else if (rA >= 0)                      { bestA = rA; bestOff = rOff; bestRc = 1; }
                else                                    bestA = -1;
                if (bestA >= 0) { host[b] = bestA; hoff[b] = bestOff; hrc[b] = (uint8_t)bestRc;
                                  hlen[b] = (uint32_t)B.size(); absorbed[b] = 1; }
            }
            // Rebuild survivors; resolve host chains (applying RC transforms).
            std::vector<int> old2new(MC, -1);
            std::vector<std::string> keep;
            for (int c = 0; c < MC; ++c) if (!absorbed[c]) { old2new[c] = (int)keep.size(); keep.push_back(contigs[c]); }
            auto resolve = [&](int c, uint32_t p, uint8_t rc, uint32_t rlen, int& nid, uint32_t& np, uint8_t& nrc){
                int cur = c; uint32_t pos = p; uint8_t r = rc;
                while (absorbed[cur]) {
                    if (hrc[cur]) {                       // B absorbed RC → flip
                        uint32_t base = hoff[cur] + hlen[cur];
                        uint32_t sub  = pos + rlen;       // per-read length
                        pos = (base >= sub) ? (base - sub) : 0;
                        r ^= 1;
                    } else pos = hoff[cur] + pos;
                    cur = host[cur];
                }
                nid = old2new[cur]; np = pos; nrc = r;
            };
            int nabs = 0;
            for (int oi = 0; oi < n; ++oi) {
                int nid; uint32_t np; uint8_t nrc;
                resolve((int)pl_cid[oi], pl_pos[oi], pl_rc[oi],
                        (uint32_t)reads[(size_t)oi].seq.size(), nid, np, nrc);
                pl_cid[oi] = (uint32_t)nid; pl_pos[oi] = np; pl_rc[oi] = nrc;
            }
            for (int c = 0; c < MC; ++c) if (absorbed[c]) ++nabs;
            _mmark("containment pass");
            size_t before = contigs.size(); contigs.swap(keep);
            fprintf(stderr, "[MERGE] contain pass %d: %zu → %zu (absorbed %d)\n", cpass, before, contigs.size(), nabs);
            if (nabs == 0) break;
        }
    }

    if (ASM_TIMING) { auto t=_tnow(); fprintf(stderr,"[ASM] Phase1.6 merge: %.2fs\n", std::chrono::duration<double>(t-_tp0).count()); _tp0=t; }

    // ── Phase 1.5: consensus polish ───────────────────────────────────────────
    // Each contig base was set by the first read to cover it (carrying that
    // read's errors). Replace it with the majority vote of all reads covering
    // the position → mismatches drop toward the true per-read error rate,
    // shrinking the aux stream. Positions are unchanged; phase 2 recomputes all
    // mismatches vs the polished contigs, so this stays lossless.
    {
        std::vector<std::vector<std::array<uint16_t,4>>> votes(contigs.size());
        for (size_t c = 0; c < contigs.size(); ++c)
            votes[c].assign(contigs[c].size(), std::array<uint16_t,4>{0,0,0,0});
        for (int oi = 0; oi < n; ++oi) {
            uint32_t cid = pl_cid[oi], pos = pl_pos[oi];
            const std::string tt = pl_rc[oi]
                ? reverse_complement(reads[(size_t)oi].seq) : reads[(size_t)oi].seq;
            auto& vc = votes[cid];
            const int rl = (int)tt.size();               // per-read length
            for (int j = 0; j < rl; ++j) {
                uint8_t b = encode_base(tt[(size_t)j]);
                uint32_t p = pos + (uint32_t)j;
                if (b < 4 && p < vc.size() && vc[p][b] < 60000) ++vc[p][b];
            }
        }
        for (size_t c = 0; c < contigs.size(); ++c) {
            auto& vc = votes[c];
            std::string& cs = contigs[c];
            for (size_t p = 0; p < cs.size(); ++p) {
                int best = 0, bv = vc[p][0];
                for (int b = 1; b < 4; ++b) if (vc[p][b] > bv) { bv = vc[p][b]; best = b; }
                if (bv > 0) cs[p] = BASE_TO_CHAR[best];
            }
        }
    }

    if (ASM_TIMING) { auto t=_tnow(); fprintf(stderr,"[ASM] Phase1.5 polish: %.2fs\n", std::chrono::duration<double>(t-_tp0).count()); _tp0=t; }

    // ── Phase 2: concatenate contigs → pg, emit per-read streams ──────────────
    std::vector<uint32_t> offset(contigs.size());
    size_t run = 0;
    for (size_t c = 0; c < contigs.size(); ++c) { offset[c] = (uint32_t)run; run += contigs[c].size(); }
    r.pg.reserve(run);
    for (auto& c : contigs) r.pg += c;

    // Global pg position of each read.
    std::vector<uint32_t> gpos(n);
    for (int oi = 0; oi < n; ++oi) gpos[oi] = offset[pl_cid[oi]] + pl_pos[oi];

    // Emit reads sorted by pg position so pg_pos is monotonic → the position
    // stream delta-codes to almost nothing. Names/quality follow chain_order,
    // so the output is a valid reordering of the read set (lossless).
    std::vector<uint32_t> ord(n);
    std::iota(ord.begin(), ord.end(), 0u);
    std::sort(ord.begin(), ord.end(),
              [&](uint32_t a, uint32_t b){ return gpos[a] < gpos[b]; });

    for (uint32_t oi : ord) {
        int rc = pl_rc[oi];
        std::string target = rc ? reverse_complement(reads[(size_t)oi].seq) : reads[(size_t)oi].seq;
        record_mapped(r, oi, gpos[oi], rc, target, L);
        if (pl_new[oi]) ++r.n_chain_starts; else ++r.n_deltas;
    }

    // ── Instrumentation: export polished contigs as FASTA (reference-free align) ──
    if (const char* cp = getenv("ARCS_CONTIGS_DUMP")) {
        FILE* cf = fopen(cp, "w");
        if (cf) {
            for (size_t c = 0; c < contigs.size(); ++c)
                fprintf(cf, ">contig_%zu\n%s\n", c, contigs[c].c_str());
            fclose(cf);
            fprintf(stderr, "[CONTIGS-DUMP] wrote %zu contigs to %s\n", contigs.size(), cp);
        }
    }

    // ── Expose placements for in-process reference-free calling (arcs call) ──────
    // Same (post-merge) data ARCS_PILEUP_SAM emits, handed to the C++ caller so it
    // needs no external aligner. Indexed by original read index.
    if (call_out) {
        call_out->contigs = contigs;
        call_out->read_cid.assign((size_t)n, 0);
        call_out->read_pos.assign((size_t)n, 0);
        call_out->read_rc.assign((size_t)n, 0);
        for (int oi = 0; oi < n; ++oi) {
            call_out->read_cid[(size_t)oi] = pl_cid[oi];
            call_out->read_pos[(size_t)oi] = pl_pos[oi];
            call_out->read_rc[(size_t)oi]  = pl_rc[oi];
        }
        call_out->valid = true;
    }

    // ── Self-contained pileup SAM: read→consensus alignments from the placements ──
    // ARCS already computes each read's placement (contig, position, strand) during
    // compression. Emitting it as a SAM gives the reference-free variant caller its
    // read-vs-consensus pileup WITHOUT any external aligner (bwa) — the "zero-cost
    // byproduct" made literal. Contig names match ARCS_CONTIGS_DUMP (contig_N) so the
    // same FASTA serves for downstream coordinate mapping. Placement is ungapped, so
    // CIGAR is {len}M; SEQ/QUAL are emitted in contig (forward) frame. MAPQ=60 (placed).
    // Off unless ARCS_PILEUP_SAM=<path>. Recomputes nothing that affects the archive.
    if (const char* sp = getenv("ARCS_PILEUP_SAM")) {
        FILE* sf = fopen(sp, "w");
        if (sf) {
            for (size_t c = 0; c < contigs.size(); ++c)
                fprintf(sf, "@SQ\tSN:contig_%zu\tLN:%zu\n", c, contigs[c].size());
            for (int oi = 0; oi < n; ++oi) {
                const std::string& os = reads[(size_t)oi].seq;
                const std::string& qs = reads[(size_t)oi].qual;
                const int rl = (int)os.size();
                const bool rcf = pl_rc[oi] != 0;
                std::string sseq = rcf ? reverse_complement(os) : os;       // contig frame
                std::string squal = qs;
                if (rcf) std::reverse(squal.begin(), squal.end());          // quality follows read
                if ((int)squal.size() != rl) squal.assign((size_t)rl, 'I'); // length guard
                // MAPQ from placement mismatch count: ARCS placement is ungapped and can
                // force-place reads a gapped aligner would soft-clip or reject, injecting
                // spurious pileup columns. Mismatch count IS the placement confidence
                // (bwa-mapq analog); mapping it to MAPQ lets the caller's existing mapq<20
                // filter drop the noisy tail. Linear/principled (NOT tuned per region).
                const std::string& cc = contigs[pl_cid[oi]];
                uint32_t cs = pl_pos[oi]; int mm = 0;
                for (int j = 0; j < rl; ++j) {
                    uint32_t p = cs + (uint32_t)j;
                    if (p >= cc.size()) { mm = rl; break; }        // ran off contig end → untrusted
                    char a = sseq[(size_t)j];
                    if (encode_base(a) < 4 && cc[p] != a) ++mm;
                }
                int mapq = 60 - 6 * mm; if (mapq < 0) mapq = 0;    // mm>6 → mapq<20 (filtered)
                fprintf(sf, "r%d\t%d\tcontig_%u\t%u\t%d\t%dM\t*\t0\t0\t%s\t%s\n",
                        oi, rcf ? 16 : 0, pl_cid[oi], pl_pos[oi] + 1u, mapq, rl,
                        sseq.c_str(), squal.c_str());
            }
            fclose(sf);
            fprintf(stderr, "[PILEUP-SAM] wrote %d read alignments to %s\n", n, sp);
        }
    }

    // ── Instrumentation: pileup / consensus-confidence dump (measurement only) ──
    // Writes per-(read,position) TSV: quality, column-agreement, depth, mismatch,
    // pos-in-read. Lets us measure whether consensus confidence predicts quality
    // beyond local context, and whether it separates errors from variants.
    // Off unless ARCS_PILEUP_DUMP=<path>. Recomputes votes; no effect on output.
    if (const char* dp = getenv("ARCS_PILEUP_DUMP")) {
        std::vector<std::vector<std::array<uint32_t,4>>> V(contigs.size());
        for (size_t c = 0; c < contigs.size(); ++c)
            V[c].assign(contigs[c].size(), std::array<uint32_t,4>{0,0,0,0});
        for (int oi = 0; oi < n; ++oi) {
            const std::string tt = pl_rc[oi] ? reverse_complement(reads[(size_t)oi].seq)
                                             : reads[(size_t)oi].seq;
            auto& vc = V[pl_cid[oi]];
            for (int j = 0; j < (int)tt.size(); ++j) {
                uint8_t b = encode_base(tt[(size_t)j]); uint32_t p = pl_pos[oi] + (uint32_t)j;
                if (b < 4 && p < vc.size()) ++vc[p][b];
            }
        }
        FILE* f = fopen(dp, "w");
        if (f) {
            // One row per orig-read position. Adds global column id + consensus base
            // + read base (target frame) so mismatches can be grouped by column to
            // test the error-vs-variant signature (recurring consistent alt = variant).
            fprintf(f, "rid\tq\tprevq\tb0\tb1\tb2\tagree\tdepth\tmm\tpos\tcol\tcons\trb\n");
            for (int oi = 0; oi < n; ++oi) {
                const std::string& qs = reads[(size_t)oi].qual;
                const std::string& os = reads[(size_t)oi].seq;    // orig frame
                uint32_t cid = pl_cid[oi], pos = pl_pos[oi];
                int prevq = 0;
                for (int k = 0; k < L; ++k) {                     // k = orig position
                    int jt = pl_rc[oi] ? (L - 1 - k) : k;         // target-frame index
                    uint32_t p = pos + (uint32_t)jt;
                    int q = (k < (int)qs.size() && qs[(size_t)k] >= 33) ? qs[(size_t)k]-33 : 0;
                    int b0 = (k>0)   ? encode_base(os[(size_t)k-1]) : 4;
                    int b1 = encode_base(os[(size_t)k]);
                    int b2 = (k+1<L) ? encode_base(os[(size_t)k+1]) : 4;
                    int agree=0; uint32_t depth=0; int mm=0; int cb=4, rb=4;
                    long col = -1;
                    if (p < V[cid].size()) {
                        auto& v = V[cid][p];
                        depth = v[0]+v[1]+v[2]+v[3];
                        uint32_t best = std::max(std::max(v[0],v[1]), std::max(v[2],v[3]));
                        agree = depth ? (int)(1000ull*best/depth) : 0;
                        cb = encode_base(contigs[cid][p]);
                        rb = encode_base(pl_rc[oi] ? "TGCAN"[b1<4?b1:4] : os[(size_t)k]);
                        mm = (rb < 4 && cb < 4 && rb != cb) ? 1 : 0;
                        col = (long)offset[cid] + (long)p;           // global column id
                    }
                    fprintf(f, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%u\t%d\t%d\t%ld\t%d\t%d\t%d\n",
                            oi, q, prevq, b0, b1, b2, agree, depth, mm, k, col, cb, rb, pl_rc[oi]);
                    prevq = q;
                }
            }
            fclose(f);
            fprintf(stderr, "[PILEUP-DUMP] wrote %s\n", dp);
        }
    }

    fprintf(stderr,
        "[CHAIN-PG-MC] reads=%d contigs=%zu placed_mapped=%zu (%.1f%%) new=%zu (%.1f%%) "
        "pg_len=%zu B (%.2fx smaller than raw reads)\n",
        n, contigs.size(), r.n_deltas, 100.0 * (double)r.n_deltas / n,
        r.n_chain_starts, 100.0 * (double)r.n_chain_starts / n,
        r.pg.size(), (double)((size_t)n * L) / (double)std::max<size_t>(1, r.pg.size()));

    return r;
}
