# ARCS Benchmark Dataset — FINAL LOCKED 2026-08-23

**DO NOT CHANGE THIS FILE.** All accessions verified via NCBI SRA 2026-08-22/23.

## The 10 Datasets

| # | Accession | Download cmd | Organism | Kingdom | Read len | _1 size (uncompressed FASTQ) | Auto-chunk |
|---|-----------|-------------|----------|---------|---------|------------------------------|-----------|
| 1 | SRR2584863 | fasterq-dump | E. coli B REL606 WGS | Bacteria | ~108 bp | ~576 MB | No |
| 2 | ERR552797 | fasterq-dump | M. tuberculosis H37Rv WGS | Bacteria | ~301 bp | ~217 MB | No |
| 3 | SRR554369 | fasterq-dump | P. aeruginosa PAO1 WGS | Bacteria | 100 bp | ~334 MB | No |
| 4 | ERR5181310 | fasterq-dump | SARS-CoV-2 amplicon WGS | Virus | ~150 bp | ~30 MB | No |
| 5 | (GIAB S3) | aws s3 cp | H. sapiens HG002 chr20 | Animalia | 150 bp | ~37 MB | No |
| 6 | SRR065390 | fasterq-dump | C. elegans N2 WGS | Animalia | 101 bp | ~7.8 GB | Yes |
| 7 | SRR327342 | fasterq-dump | S. cerevisiae I14 WGS | Fungi | ~100 bp | ~3.4 GB | Yes |
| 8 | SRR1945765 | fasterq-dump | Arabidopsis thaliana WGS | Plantae | 102 bp | ~1.95 GB | No |
| 9 | SRR1663585 | fasterq-dump | D. melanogaster WGS | Animalia | 101 bp | ~2.33 GB | No |
| 10 | SRR870667 | fasterq-dump | T. cacao Matina1-6 WGS | Plantae | 109 bp | ~17 GB | Yes |

**Diversity:** 3 bacteria · 1 virus · 3 animals · 1 fungus · 2 plants

**Sizes are uncompressed FASTQ _1 file sizes** (SRA download is 3-5× smaller).
Auto-chunk fires automatically in ARCS for files >2 GB. This is expected default behaviour.

## Download Command

```bash
bash benchmark/download.sh /data/fastq
```

## ARCS Compression Command (NO FLAGS — default only)

```bash
arcs compress INPUT.fq OUTPUT.arcs
```

Do NOT add any flags. Auto-chunk is ARCS default.

## Disk space estimate (250 GB SSD)

| Files | Size |
|-------|------|
| SRA downloads (compressed) | ~30 GB |
| Uncompressed FASTQ _1 files | ~33 GB |
| ARCS archives + intermediates | ~15 GB |
| SPRING/Genozip archives | ~20 GB |
| **Total peak** | **~98 GB** |

250 GB SSD is sufficient with ~150 GB headroom.

## Datasets NEVER to use again (all wrong — RNA-Seq or mis-annotated)

- SRR390728 — H. sapiens RNA-Seq B-cell lymphoma (NOT WGS)
- SRR1296601 — Glycine max ncRNA-Seq (NOT M. tuberculosis)
- ERR015526 — H. sapiens WGS 109bp (NOT E. coli 36bp)
- SRR1294122 — H. sapiens RNA-Seq stem cells (NOT WGS)
- ERR174310 — H. sapiens WGS NA12877 26.6 GB (NOT C. elegans, too large)
- SRR988075 — D. melanogaster, 8.7 GB _1 (replaced by SRR1663585, 2.33 GB)
- SRR065390 was previously mis-noted as "E. coli" — it IS C. elegans ✅ (confirmed NCBI)

## Published SPRING bpb (PgRC2 2025, for cross-validation)

| Accession | SPRING bpb | Notes |
|-----------|-----------|-------|
| SRR554369 | 0.2416 | P. aeruginosa |
| SRR870667 | 1.2621 | T. cacao — SPRING very poor on complex plant |

Cross-validate: server SPRING measurement for these two should match ±2%.
