# ARCS full benchmark — 8 datasets × 5 tools

Every value links to the raw log it was measured from. Native ext4, best-of-2 wall time, peak RAM via kernel VmHWM. Lossless verified by byte comparison.

Datasets: **DS1** (E.coli 151bp), **DS2** (Human 127bp), **DS4** (M.tb 51bp), **DS5** (SARS 221bp), **DS6** (bacteria 37bp), **DS7** (Human chr22 ~1M reads), **GIAB** (HG002 150bp), **NA12878** (HG001 148bp)

## Table 1 — Archive size (compression ratio; smaller = better)

| Dataset | raw | ARCS | SPRING | Genozip | fqzcomp | PgRC2 |
|---|---|---|---|---|---|---|
| DS1 | 37.04 MB | [6.25 MB](logs/DS1__ARCS.log) | [6.50 MB](logs/DS1__SPRING.log) | [8.37 MB](logs/DS1__Genozip.log) | [7.82 MB](logs/DS1__fqzcomp.log) | — |
| DS2 | 28.58 MB | [4.82 MB](logs/DS2__ARCS.log) | [4.92 MB](logs/DS2__SPRING.log) | [5.08 MB](logs/DS2__Genozip.log) | [4.91 MB](logs/DS2__fqzcomp.log) | — |
| DS4 | 16.98 MB | [3.06 MB](logs/DS4__ARCS.log) | [3.13 MB](logs/DS4__SPRING.log) | [3.23 MB](logs/DS4__Genozip.log) | [3.25 MB](logs/DS4__fqzcomp.log) | — |
| DS5 | 48.84 MB | [1.37 MB](logs/DS5__ARCS.log) | [1.39 MB](logs/DS5__SPRING.log) | [1.35 MB](logs/DS5__Genozip.log) | [1.39 MB](logs/DS5__fqzcomp.log) | — |
| DS6 | 10.48 MB | [0.60 MB](logs/DS6__ARCS.log) | [0.73 MB](logs/DS6__SPRING.log) | [0.75 MB](logs/DS6__Genozip.log) | [0.82 MB](logs/DS6__fqzcomp.log) | — |
| DS7 | 343.71 MB | [48.84 MB](logs/DS7__ARCS.log) | [67.53 MB](logs/DS7__SPRING.log) | [58.97 MB](logs/DS7__Genozip.log) | [86.49 MB](logs/DS7__fqzcomp.log) | — |
| GIAB | 39.38 MB | [4.24 MB](logs/GIAB__ARCS.log) | [5.51 MB](logs/GIAB__SPRING.log) | [5.02 MB](logs/GIAB__Genozip.log) | [7.26 MB](logs/GIAB__fqzcomp.log) | — |
| NA12878 | 29.02 MB | [4.08 MB](logs/NA12878__ARCS.log) | [5.62 MB](logs/NA12878__SPRING.log) | [4.83 MB](logs/NA12878__Genozip.log) | [6.86 MB](logs/NA12878__fqzcomp.log) | — |

## Table 2a — Compression speed (wall time; lower = better)

| Dataset | ARCS | SPRING | Genozip | fqzcomp | PgRC2 |
|---|---|---|---|---|---|
| DS1 | [59.23s](logs/DS1__ARCS.log) | [5.28s](logs/DS1__SPRING.log) | [1.19s](logs/DS1__Genozip.log) | [0.46s](logs/DS1__fqzcomp.log) | — |
| DS2 | [96.27s](logs/DS2__ARCS.log) | [5.73s](logs/DS2__SPRING.log) | [1.16s](logs/DS2__Genozip.log) | [0.26s](logs/DS2__fqzcomp.log) | — |
| DS4 | [23.98s](logs/DS4__ARCS.log) | [3.26s](logs/DS4__SPRING.log) | [0.74s](logs/DS4__Genozip.log) | [0.25s](logs/DS4__fqzcomp.log) | — |
| DS5 | [10.47s](logs/DS5__ARCS.log) | [2.70s](logs/DS5__SPRING.log) | [2.93s](logs/DS5__Genozip.log) | [0.26s](logs/DS5__fqzcomp.log) | — |
| DS6 | [11.79s](logs/DS6__ARCS.log) | [2.34s](logs/DS6__SPRING.log) | [0.53s](logs/DS6__Genozip.log) | [0.10s](logs/DS6__fqzcomp.log) | — |
| DS7 | [180.49s](logs/DS7__ARCS.log) | [37.12s](logs/DS7__SPRING.log) | [87.71s](logs/DS7__Genozip.log) | [7.60s](logs/DS7__fqzcomp.log) | — |
| GIAB | [12.13s](logs/GIAB__ARCS.log) | [5.54s](logs/GIAB__SPRING.log) | [7.18s](logs/GIAB__Genozip.log) | [0.74s](logs/GIAB__fqzcomp.log) | — |
| NA12878 | [10.56s](logs/NA12878__ARCS.log) | [5.12s](logs/NA12878__SPRING.log) | [6.07s](logs/NA12878__Genozip.log) | [0.55s](logs/NA12878__fqzcomp.log) | — |

## Table 2b — Decompression speed (wall time; lower = better)

| Dataset | ARCS | SPRING | Genozip | fqzcomp | PgRC2 |
|---|---|---|---|---|---|
| DS1 | [15.27s](logs/DS1__ARCS.log) | [2.17s](logs/DS1__SPRING.log) | [0.41s](logs/DS1__Genozip.log) | [0.60s](logs/DS1__fqzcomp.log) | — |
| DS2 | [20.08s](logs/DS2__ARCS.log) | [2.31s](logs/DS2__SPRING.log) | [0.38s](logs/DS2__Genozip.log) | [0.41s](logs/DS2__fqzcomp.log) | — |
| DS4 | [5.94s](logs/DS4__ARCS.log) | [1.02s](logs/DS4__SPRING.log) | [0.44s](logs/DS4__Genozip.log) | [0.26s](logs/DS4__fqzcomp.log) | — |
| DS5 | [1.51s](logs/DS5__ARCS.log) | [1.98s](logs/DS5__SPRING.log) | [0.33s](logs/DS5__Genozip.log) | [0.36s](logs/DS5__fqzcomp.log) | — |
| DS6 | [3.15s](logs/DS6__ARCS.log) | [0.55s](logs/DS6__SPRING.log) | [0.21s](logs/DS6__Genozip.log) | [0.12s](logs/DS6__fqzcomp.log) | — |
| DS7 | [21.47s](logs/DS7__ARCS.log) | [13.52s](logs/DS7__SPRING.log) | [8.19s](logs/DS7__Genozip.log) | [9.06s](logs/DS7__fqzcomp.log) | — |
| GIAB | [2.21s](logs/GIAB__ARCS.log) | [3.45s](logs/GIAB__SPRING.log) | [0.52s](logs/GIAB__Genozip.log) | [0.86s](logs/GIAB__fqzcomp.log) | — |
| NA12878 | [3.66s](logs/NA12878__ARCS.log) | [2.48s](logs/NA12878__SPRING.log) | [0.78s](logs/NA12878__Genozip.log) | [0.61s](logs/NA12878__fqzcomp.log) | — |

## Table 3 — Peak RAM (compress / decompress, VmHWM)

| Dataset | ARCS | SPRING | Genozip | fqzcomp | PgRC2 |
|---|---|---|---|---|---|
| DS1 | [2002 MB](logs/DS1__ARCS.log) / 654 MB | [378 MB](logs/DS1__SPRING.log) / 120 MB | [275 MB](logs/DS1__Genozip.log) / 139 MB | [61 MB](logs/DS1__fqzcomp.log) / 59 MB | — |
| DS2 | [2096 MB](logs/DS2__ARCS.log) / 661 MB | [350 MB](logs/DS2__SPRING.log) / 97 MB | [221 MB](logs/DS2__Genozip.log) / 120 MB | [60 MB](logs/DS2__fqzcomp.log) / 58 MB | — |
| DS4 | [1107 MB](logs/DS4__ARCS.log) / 329 MB | [307 MB](logs/DS4__SPRING.log) / 55 MB | [168 MB](logs/DS4__Genozip.log) / 80 MB | [60 MB](logs/DS4__fqzcomp.log) / 58 MB | — |
| DS5 | [383 MB](logs/DS5__ARCS.log) / 161 MB | [421 MB](logs/DS5__SPRING.log) / 145 MB | [444 MB](logs/DS5__Genozip.log) / 117 MB | [57 MB](logs/DS5__fqzcomp.log) / 57 MB | — |
| DS6 | [521 MB](logs/DS6__ARCS.log) / 167 MB | [289 MB](logs/DS6__SPRING.log) / 43 MB | [105 MB](logs/DS6__Genozip.log) / 50 MB | [57 MB](logs/DS6__fqzcomp.log) / 57 MB | — |
| DS7 | [3912 MB](logs/DS7__ARCS.log) / 1468 MB | [1306 MB](logs/DS7__SPRING.log) / 1004 MB | [1448 MB](logs/DS7__Genozip.log) / 610 MB | [62 MB](logs/DS7__fqzcomp.log) / 60 MB | — |
| GIAB | [423 MB](logs/GIAB__ARCS.log) / 120 MB | [382 MB](logs/GIAB__SPRING.log) / 126 MB | [413 MB](logs/GIAB__Genozip.log) / 120 MB | [60 MB](logs/GIAB__fqzcomp.log) / 59 MB | — |
| NA12878 | [474 MB](logs/NA12878__ARCS.log) / 215 MB | [349 MB](logs/NA12878__SPRING.log) / 94 MB | [341 MB](logs/NA12878__Genozip.log) / 130 MB | [61 MB](logs/NA12878__fqzcomp.log) / 59 MB | — |

## Table 5 — Losslessness (byte-identical roundtrip)

| Dataset | ARCS | SPRING | Genozip | fqzcomp | PgRC2 |
|---|---|---|---|---|---|
| DS1 | [LOSSLESS](logs/DS1__ARCS.log) | [LOSSLESS](logs/DS1__SPRING.log) | [LOSSLESS](logs/DS1__Genozip.log) | [LOSSY](logs/DS1__fqzcomp.log) | — |
| DS2 | [LOSSLESS](logs/DS2__ARCS.log) | [LOSSLESS](logs/DS2__SPRING.log) | [LOSSLESS](logs/DS2__Genozip.log) | [LOSSY](logs/DS2__fqzcomp.log) | — |
| DS4 | [LOSSLESS](logs/DS4__ARCS.log) | [LOSSLESS](logs/DS4__SPRING.log) | [LOSSLESS](logs/DS4__Genozip.log) | [LOSSY](logs/DS4__fqzcomp.log) | — |
| DS5 | [LOSSLESS](logs/DS5__ARCS.log) | [LOSSLESS](logs/DS5__SPRING.log) | [LOSSLESS](logs/DS5__Genozip.log) | [LOSSLESS](logs/DS5__fqzcomp.log) | — |
| DS6 | [LOSSLESS](logs/DS6__ARCS.log) | [LOSSLESS](logs/DS6__SPRING.log) | [LOSSLESS](logs/DS6__Genozip.log) | [LOSSY](logs/DS6__fqzcomp.log) | — |
| DS7 | [LOSSLESS](logs/DS7__ARCS.log) | [LOSSLESS](logs/DS7__SPRING.log) | [LOSSLESS](logs/DS7__Genozip.log) | [LOSSY](logs/DS7__fqzcomp.log) | — |
| GIAB | [LOSSLESS](logs/GIAB__ARCS.log) | [LOSSLESS](logs/GIAB__SPRING.log) | [LOSSLESS](logs/GIAB__Genozip.log) | [LOSSY](logs/GIAB__fqzcomp.log) | — |
| NA12878 | [LOSSLESS](logs/NA12878__ARCS.log) | [LOSSLESS](logs/NA12878__SPRING.log) | [LOSSLESS](logs/NA12878__Genozip.log) | [LOSSY](logs/NA12878__fqzcomp.log) | — |

## Table 4 — Reference-free variant calling (analysis-free byproduct)

Only datasets with a truth VCF (GIAB, NA12878, E.coli-sim). Values from the ref-free eval docs; see linked result files.

| Dataset | ARCS F1 | Kmer2SNP F1 | DiscoSNP++ | source |
|---|---|---|---|---|
| NA12878 (held-out) | 0.816 | 0.833 | 0.71 | [REFFREE](../docs/REFFREE_NA12878_COMPARISON.md) |
| NA12878 (avg 3-region) | 0.840 | 0.818 | 0.71 | [REFFREE](../docs/REFFREE_NA12878_COMPARISON.md) |
| E.coli sim (haploid) | >0.99 P/R | — | — | [SHARED_LATENT](../docs/SHARED_LATENT_ERROR_STATE.md) |
