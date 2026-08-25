# ARCS Benchmark Dataset — FINAL LOCKED 2026-08-25

**DO NOT CHANGE THIS FILE.** All accessions verified via NCBI SRA/DDBJ 2026-08-25.

## The 10 Datasets

| # | Accession | Organism | Kingdom | Read len | _1 size (uncompressed) | Auto-chunk |
|---|-----------|----------|---------|---------|------------------------|-----------|
| 1 | SRR2584863 | E. coli B REL606 WGS | Bacteria | ~108 bp | ~576 MB | No |
| 2 | ERR552797 | M. tuberculosis H37Rv WGS | Bacteria | ~301 bp | ~217 MB | No |
| 3 | SRR554369 | P. aeruginosa PAO1 WGS | Bacteria | 100 bp | ~334 MB | No |
| 4 | ERR5181310 | SARS-CoV-2 amplicon WGS | Virus | ~150 bp | ~30 MB | No |
| 5 | (GIAB S3) | H. sapiens HG002 chr20 | Animalia | 150 bp | ~37 MB | No |
| 6 | SRR16357346 | C. elegans N2 WGS (CeNDR) | Animalia | ~135 bp | ~1.42 GB | No |
| 7 | DRR976266 | S. cerevisiae WGS | Fungi | ~150 bp | ~1.67 GB | No |
| 8 | SRR1945765 | Arabidopsis thaliana WGS | Plantae | 102 bp | ~1.95 GB | No |
| 9 | SRR36741279 | Leishmania major WGS | Protista | ~75 bp | ~1.70 GB | No |
| 10 | SRR37283774 | P. falciparum WGS | Protista | ~100 bp | ~669 MB | No |

**Diversity:** 3 bacteria · 1 virus · 2 animals · 1 fungus · 1 plant · 2 protists — 6 kingdoms

**All files are under 2 GB. Auto-chunk does NOT fire for any dataset. ARCS runs single-pass assembly on all 10.**

**Sizes are uncompressed FASTQ _1 file sizes** (SRA download is 3-5× smaller).

## Spot counts and coverage

| # | Accession | Spots | Read len | Genome | Coverage |
|---|-----------|-------|----------|--------|----------|
| 1 | SRR2584863 | ~2.5M | 108 bp | 4.6 Mb | ~59× |
| 2 | ERR552797 | ~720K | 301 bp | 4.4 Mb | ~49× |
| 3 | SRR554369 | ~1.5M | 100 bp | 6.3 Mb | ~24× |
| 4 | ERR5181310 | ~100K | 150 bp | 30 kb | ~500× |
| 5 | GIAB HG002 chr20 | ~240K | 150 bp | 50 Mb | ~0.7× subset |
| 6 | SRR16357346 | 4.75M | ~135 bp | 100 Mb | ~6.4× |
| 7 | DRR976266 | 5.08M | ~150 bp | 12 Mb | ~63× |
| 8 | SRR1945765 | ~8.6M | 102 bp | 135 Mb | ~6.5× |
| 9 | SRR36741279 | 9.49M | ~75 bp | 32 Mb | ~22× |
| 10 | SRR37283774 | 2.91M | ~100 bp | 23 Mb | ~12.6× |

## Download Command

```bash
bash benchmark/download.sh /data/fastq
```

## ARCS Compression Command (NO FLAGS — default only)

```bash
arcs compress INPUT.fq OUTPUT.arcs
```

No flags. No auto-chunk on any dataset. Single-pass assembly on all 10.

## Disk space estimate (250 GB SSD)

| Files | Size |
|-------|------|
| SRA downloads (compressed) | ~12 GB |
| Uncompressed FASTQ _1 files | ~9 GB |
| ARCS archives + intermediates | ~5 GB |
| SPRING/Genozip archives | ~8 GB |
| **Total peak** | **~34 GB** |

250 GB SSD has ~216 GB headroom.

## Banned accessions — never use again

- SRR390728 — H. sapiens RNA-Seq (NOT WGS)
- SRR1296601 — Glycine max ncRNA-Seq (NOT M. tuberculosis)
- ERR015526 — H. sapiens WGS (NOT E. coli 36bp)
- SRR1294122 — H. sapiens RNA-Seq stem cells (NOT WGS)
- ERR174310 — H. sapiens WGS NA12877 26.6 GB (NOT C. elegans, too large)
- SRR988075 — D. melanogaster 8.7 GB _1 (too large, replaced)
- SRR065390 — C. elegans 7.8 GB _1 (too large, replaced by SRR16357346)
- SRR327342 — S. cerevisiae 3.4 GB _1 (too large, replaced by DRR976266)
- SRR1663585 — D. melanogaster 2.33 GB _1 (auto-chunk fires, replaced by SRR36741279)
- SRR870667 — T. cacao 17 GB _1 (too large, replaced by SRR37283774)

## Published SPRING bpb (PgRC2 2025, for cross-validation)

| Accession | SPRING bpb | Notes |
|-----------|-----------|-------|
| SRR554369 | 0.2416 | P. aeruginosa — use to verify SPRING is working correctly |

Cross-validate: server SPRING bpb for SRR554369 must be ~0.2416 ±2%. T. cacao gold bpb no longer applies (accession replaced).
