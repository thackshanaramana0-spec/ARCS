# ARCS Benchmark Dataset — FINAL LOCKED 2026-08-25

**DO NOT CHANGE THIS FILE.** All accessions verified via NCBI SRA/DDBJ 2026-08-25.

---

## The 10 Full Datasets (Claim 1 primary + Claim 3)

All are complete, unsubsampled SRA runs as deposited by sequencing labs.

| # | Accession | Organism | Kingdom | Read len | _1 size (uncompressed) | Peak RAM | Auto-chunk |
|---|-----------|----------|---------|---------|------------------------|---------|-----------|
| 1 | SRR2584863 | E. coli B REL606 WGS | Bacteria | ~108 bp | ~576 MB | ~5 GB | No |
| 2 | ERR552797 | M. tuberculosis H37Rv WGS | Bacteria | ~301 bp | ~217 MB | ~3 GB | No |
| 3 | SRR554369 | P. aeruginosa PAO1 WGS | Bacteria | 100 bp | ~334 MB | ~4 GB | No |
| 4 | ERR5181310 | SARS-CoV-2 amplicon WGS | Virus | ~150 bp | ~30 MB | ~1 GB | No |
| 5 | ERR17740259 | S. aureus WGS (Firmicutes, 33% GC) | Bacteria | ~148 bp | ~970 MB | ~6 GB | No |
| 6 | SRR065390 | C. elegans N2 WGS | Animalia | 100 bp | ~11 GB | ~18 GB | Suppressed* |
| 7 | DRR976266 | S. cerevisiae WGS | Fungi | ~150 bp | ~1.67 GB | ~8 GB | No |
| 8 | SRR870667 | Theobroma cacao WGS | Plantae | 69 bp SE | ~15 GB | ~28 GB | Suppressed* |
| 9 | SRR36741279 | Leishmania major WGS | Protista | ~75 bp | ~1.70 GB | ~7 GB | No |
| 10 | SRR37283774 | P. falciparum WGS | Protista | ~100 bp | ~669 MB | ~5 GB | No |

**\* Auto-chunk suppressed via `ARCS_AUTOCHUNK_MB=25000` env var in run_block1.sh. Server has 90 GB RAM; both datasets fit single-pass assembly comfortably.**

**Diversity:** 4 bacteria · 1 virus · 1 animal · 1 fungus · 1 plant · 2 protists — 6 kingdoms  
**Note:** S. aureus (slot 5) adds Firmicutes / Gram-positive / 33% GC — distinct from slots 1–3 (all Gram-negative, 51–67% GC).  
**Note:** Slots 6 and 8 were previously banned due to RAM constraints on 16 GB servers. On 90 GB server, peak RAM is ~18 GB and ~28 GB respectively — both safe.

---

## The 5 Human Subset Datasets (Claim 1 supplementary + Claim 2)

Full human WGS (~150 GB) cannot be single-pass assembled on current hardware; these chr20 subsets demonstrate ARCS on human data. Limitation stated in paper.

| # | File | Individual | Ancestry | Size | Source |
|---|------|-----------|---------|------|--------|
| S1 | HG001_pooled.fq | NA12878 (HG001) | European CEU | ~37 MB | GIAB S3 |
| S2 | HG002_pooled.fq | NA24385 (HG002) | Ashkenazi Jewish son | ~37 MB | GIAB S3 |
| S3 | HG003_pooled.fq | NA24149 (HG003) | Ashkenazi Jewish father | ~37 MB | GIAB S3 |
| S4 | HG004_pooled.fq | NA24143 (HG004) | Ashkenazi Jewish mother | ~37 MB | GIAB S3 |
| S5 | HG005_pooled.fq | NA24631 (HG005) | Han Chinese son | ~37 MB | GIAB S3 |

**All 5 are under 2 GB. Auto-chunk does NOT fire. ARCS runs single-pass assembly on all 5.**  
**Ancestry diversity:** European, Ashkenazi Jewish (trio), Han Chinese — demonstrates reference-free advantage on non-European genomes.

---

## Spot counts and coverage

| # | Accession | Spots | Read len | Genome | Coverage |
|---|-----------|-------|----------|--------|----------|
| 1 | SRR2584863 | ~2.5M | 108 bp | 4.6 Mb | ~59× |
| 2 | ERR552797 | ~720K | 301 bp | 4.4 Mb | ~49× |
| 3 | SRR554369 | ~1.5M | 100 bp | 6.3 Mb | ~24× |
| 4 | ERR5181310 | ~100K | 150 bp | 30 kb | ~500× |
| 5 | ERR17740259 | ~3.0M | ~148 bp | 2.8 Mb | ~158× |
| 6 | SRR065390 | ~33.8M PE | 100 bp | 100 Mb | ~34× (_1 only) |
| 7 | DRR976266 | 5.08M | ~150 bp | 12 Mb | ~63× |
| 8 | SRR870667 | ~74M SE | 69 bp | 430 Mb | ~12× |
| 9 | SRR36741279 | 9.49M | ~75 bp | 32 Mb | ~22× |
| 10 | SRR37283774 | 2.91M | ~100 bp | 23 Mb | ~12.6× |
| S1-S5 | GIAB HG001-HG005 | ~240K each | 150 bp | chr20 50 Mb | ~0.7× (subset) |

## Download Command

```bash
bash benchmark/download.sh /data/fastq
```

## ARCS Compression Command

```bash
ARCS_AUTOCHUNK_MB=25000 arcs compress INPUT.fq OUTPUT.arcs
```

No command-line flags. `ARCS_AUTOCHUNK_MB=25000` is a server-side env var (set in run_block1.sh) that raises the auto-chunk threshold to 25 GB, ensuring single-pass assembly on all 10 datasets including SRR065390 (~11 GB) and SRR870667 (~15 GB).

## Disk space estimate (250 GB SSD)

| Files | Size |
|-------|------|
| SRA downloads (compressed) | ~30 GB |
| Uncompressed FASTQ _1 files | ~32 GB |
| ARCS archives + intermediates | ~8 GB |
| SPRING/Genozip archives | ~12 GB |
| **Total peak** | **~82 GB** |

250 GB SSD has ~168 GB headroom.

## Banned accessions — never use again

- SRR390728 — H. sapiens RNA-Seq (NOT WGS)
- SRR1296601 — Glycine max ncRNA-Seq (NOT M. tuberculosis)
- ERR015526 — H. sapiens WGS (NOT E. coli 36bp)
- SRR1294122 — H. sapiens RNA-Seq stem cells (NOT WGS)
- ERR174310 — H. sapiens WGS NA12877 26.6 GB (too large, human excluded from scope)
- SRR988075 — D. melanogaster 8.7 GB _1 (no published gold bpb, replaced)
- SRR16357346 — C. elegans CeNDR 6.4× coverage (too low coverage, replaced by SRR065390)
- SRR1945765 — Arabidopsis thaliana 6.5× coverage (too low coverage, replaced by SRR870667)
- SRR327342 — S. cerevisiae 3.4 GB _1 (replaced by DRR976266 higher coverage)
- SRR1663585 — D. melanogaster 2.33 GB _1 (no published gold bpb, replaced by SRR36741279)

## Published SPRING bpb (PgRC2 2025, for cross-validation)

| Accession | SPRING bpb | Notes |
|-----------|-----------|-------|
| SRR554369 | 0.2416 | P. aeruginosa — use to verify SPRING is working correctly |
| SRR870667 | 1.2621 | T. cacao — SPRING's worst dataset; ARCS expected biggest win here |

Cross-validate: server SPRING bpb must match ±2%. If SRR870667 deviates, SPRING benchmark setup is wrong.
