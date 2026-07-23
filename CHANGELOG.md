# Changelog

All notable changes to ARCS are documented here.

---

## v2.1.0 — 2026-07-23

### Dual-mode sequence codec

- **`ARCS_PG_FAST_DECODE=1`**: LZMA-ASCII fast-decode mode (format 0x07, codec_id 0x01).
  Plain LZMA-9 on ASCII pseudogenome — no 2-bit packing. GIAB: 2.8× faster decompress
  (1.0 s vs 2.8 s), +0.23% archive cost, lossless roundtrip verified.
  Single block always — multi-block LZMA regresses on fragmented pgs (same cross-block
  context loss as multi-block FCM). Speed comes from LZMA streaming replay vs FCM
  adaptive table re-training from scratch.
- Default (FCM adaptive) unchanged: 4,310,976 B GIAB, byte-exact baseline.
- Genozip gap (0.45 s) is architectural: storing FCM state ~150 MB vs 4.3 MB archive.
  Fast-decode is the practical ceiling without trading ratio crown.

### Parallel decode pipeline

- **Step 1**: pg + names async parallel stream decode — GIAB decompress 2.15 s → 1.93 s (−10%).
- **Step 2**: Seeded multi-block format 0x06 — 26-byte context seed from preceding block
  tail warms FCM shift-registers before each block. Auto-activates at pg ≥ 100 MB (large WGS).
  DS7-scale: 4-block parallel decode, ~3–4× pg-decode speedup. Small datasets (GIAB 743 KB)
  stay single-block, zero ratio regression.
- DS1 name bug fixed: `build_columnar_names()` restored to Illumina fast path.

---

## v2.0.0 — 2026-07-22

### New default: chain-pg pseudogenome encoder

- **chain-pg is now the default encoder** (no flags required). Previous default was the
  MST encoder which produced archives 19% larger. All no-flag benchmarks now reflect
  the true best-ratio path.
- Amplicon data (AKC score < 0.15) is still correctly routed to the dedicated
  cluster-sort encoder — AKC check takes priority over the chain-pg default.
- Amplicon header flag (0x10) fixed — was incorrectly gated on `!use_chain_pg`.

### Variant calling: `arcs call`

- `arcs call reads.fq out.vcf` — self-contained reference-free heterozygous variant caller.
  Uses ARCS's own placements (CallData), internal 31-mer counting, no external aligner.
- Het-SNV F1 = 0.957 averaged across 5 GIAB individuals / 3 ancestries (HG001–HG005),
  validated by rtg vcfeval (GA4GH engine). Beats DiscoSNP++ (0.912) and Kmer2SNP (0.533).
- Het-indel calling via contig-bubble detection. Real GIAB HG002 chr20 F1 = 0.505
  (DiscoSNP++ 0.513 — zero-cost byproduct ties the dedicated tool).
- Polyploid support (`ARCS_PLOIDY=k`): triploid synthetic F1 = 0.994–1.000 (4 seeds).
- `arcs compress --call` — fused single-pass compress + call (byte-identical archive).

### Compression improvements

- **8/8 dataset wins vs SPRING** (was 7/8 before auto-chunk threshold fix).
  Auto-chunk threshold raised 200 MB → 2000 MB (prevents fragmented assembly on large files).
- **Decode speedup**: buffered FASTQWriter — write phase 0.88 s → 0.05 s, GIAB decode
  2.5 s → 1.59 s. Byte-exact, no ratio cost.
- **Order-preserving default**: byte-exact output by default; `ARCS_ORDER_FREE=1` opts
  into smaller multiset mode. CRLF preserved (`\r\n` round-trips correctly).
- **Own range-coded name coder** (format 0x05): order-0 range coder for Illumina name
  coordinates; wins DS5 over SPRING on names.
- Wins 8/8 smallest lossless archive vs SPRING:
  DS1 −3.8%, DS2 −1.9%, DS4 −1.9%, DS5 −8.6%, DS6 −17%, DS7 −28%, GIAB −23%, ECOLI30x −0.9%.

### Speed and RAM

- O(1) reverse-complement IR + prefetch (`dna_coder.cpp`): DS1 pg-decode 13.3 s → 8.0 s.
- `unordered_map` → chunked FreqMap (`qual_cm.cpp`): DS7 decode 29.8 s → 16.8 s (−44%).
- Flat CSR KmerIndex for contig MERGE (`chain_encoder.cpp`): DS1 assembly −27%.
- uint8 Row4 in `dna_coder.cpp`: DS2 decode RAM 661 MB → 341 MB (bit-identical).
- Block-parallel quality CM + names LZMA auto-scale with core count.

### Edge-case hardening

- 12 adversarial input classes now round-trip correctly: empty reads, 1 bp reads, sub-seed
  reads, 0-length reads, mixed-length, lowercase/IUPAC bases, full Phred 0–93,
  control-char names, name traps. >65535 bp reads: clean error refusal.
- SILENT BUG FIXED: quality clamped at Phred 42 — legal FASTQ is 0–93. Now encodes
  full range losslessly.

---

## v1.0.0 — 2026-06-23

Initial release.

- Reference-free lossless FASTQ compression (single-end).
- AKC pilot scoring for automatic WGS / amplicon regime detection.
- WGS path: LSH-accelerated k-NN graph → Kruskal MST → DFS ordering →
  context-rANS quality codec.
- Amplicon path: PCR dedup → cluster-order sort → LZMA-9 quality codec.
- Platforms: Windows 11 (MinGW-w64), Linux Ubuntu 22.04.
