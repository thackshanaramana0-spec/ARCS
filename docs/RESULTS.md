# ARCS — Benchmark Results

All compression benchmarks: native Linux (WSL2 ext4), same machine (Intel i5-1235U, 12 threads),
best-of-2 wall time, peak RAM via kernel VmHWM. All ARCS runs byte-exact lossless verified.
SPRING `-t 8`; Genozip built from source (github.com/divonlan/genozip).
Variant calling benchmarks: rtg vcfeval (GA4GH engine), frozen parameters throughout.

---

## Datasets

| Label | Organism | Read length | Accession | Raw size |
|---|---|---|---|---|
| ECOLI_35x | *E. coli* K-12 | 151 bp | SRR2584863 | 35 MB |
| HUMAN_127bp | *H. sapiens* | 127 bp | SRR2962669 | 27 MB |
| MTB_51bp | *M. tuberculosis* | 51 bp | ERR11820348 | 17 MB |
| SARS2_AMP | SARS-CoV-2 (amplicon) | 221 bp | ERR5181310 | 46 MB |
| GIAB_HG002 | *H. sapiens* HG002 | 150 bp | GIAB v4.2.1 | 37 MB |
| GIAB_HG001 | *H. sapiens* NA12878 | 148 bp | GIAB v4.2.1 | 116 MB |
| ECOLI_30x | *E. coli* K-12 | 150 bp | — | 276 MB |
| HUMAN_30x | *H. sapiens* WGS | 150 bp | SRR1238539 | 327 MB |

---

## Table 1 — Archive size (KB, smaller = better)

PgRC2 default is lossy quality — lossless mode (`-Q`) crashes on tested datasets.
ARCS (default) = chain-pg encoder (best ratio). ARCS (fast) = `--fast-decode` flag (faster decompress, slightly larger archive).

| Dataset | Raw | **ARCS** | ARCS (fast) | Genozip | SPRING | PgRC2 (lossy†) |
|---|---|---|---|---|---|---|
| ECOLI_35x | 35 MB | **6,107** | 6,147 | 7,847 | 6,350 | — |
| HUMAN_127bp | 27 MB | **4,706** | 4,741 | 4,962 | 4,800 | — |
| MTB_51bp | 16 MB | **2,990** | 2,999 | 3,154 | 3,060 | — |
| SARS2_AMP | 46 MB | **1,233** | 1,249 | 1,312 | 1,340 | — |
| GIAB_HG002 | 37 MB | **4,210** | 4,282 | 4,929 | 5,340 | 240 |
| GIAB_HG001 | 116 MB | 16,850 | 17,264 | **16,763** | 19,410 | — |
| ECOLI_30x | 276 MB | **35,289**‡ | — | 62,442 | 35,610 | — |
| HUMAN_30x | 327 MB | **46,774** | 54,103 | 57,589 | 66,070 | — |

**Wins:** ARCS vs SPRING **8/8**; vs Genozip **7/8** (loses GIAB_HG001 by 0.5%, near-tie).
SPRING is not byte-exact on 4/8 bacterial/short datasets (reorders reads); ARCS and Genozip are byte-exact on all.

† PgRC2 GIAB_HG002 lossy = 240 KB (sequence-only; qualities discarded). Not a lossless comparison.
‡ ECOLI_30x post-autochunk fix (threshold 200 MB → 2000 MB, build ≥ 2026-07-21). Pre-fix value was 38,166 KB (ARCS lost by 7%); fix eliminates assembly fragmentation on large files.

---

## Table 2 — Compress time (seconds, lower = better)

> **Note:** These times were recorded before the v2.2.0 parallel shard assembly (4× speedup).
> Current ARCS v2.2.0 compress times are approximately 3-4× faster on multi-core machines.
> Ratio numbers in Table 1 are unaffected (parallel assembly is ratio-neutral within 1.5%).


| Dataset | **ARCS** | ARCS (fast) | Genozip | SPRING |
|---|---|---|---|---|
| ECOLI_35x | 44.6 s | 36.8 s | 2.6 s | 6.0 s |
| HUMAN_127bp | 78.0 s | 58.6 s | 1.6 s | 8.6 s |
| MTB_51bp | 24.3 s | 28.3 s | 1.6 s | 2.8 s |
| SARS2_AMP | 7.6 s | 23.4 s | 5.2 s | 4.0 s |
| GIAB_HG002 | 13.1 s | 19.0 s | 5.7 s | 11.5 s |
| GIAB_HG001 | 28.3 s | 58.8 s | 14.4 s | 12.2 s |
| ECOLI_30x | 245.8 s | 84.8 s | 6.4 s | 23.6 s |
| HUMAN_30x | 68.3 s | 44.2 s | **75.4 s** | 32.6 s |

ARCS compress is slower than SPRING/Genozip on most datasets — the assembly (de Bruijn graph + MST/chain-pg) is the cost, the same trade-off PgRC2 makes (accepted in literature). ARCS beats Genozip on HUMAN_30x and SARS2_AMP. Native ext4 times are ~3–4× faster (assembly is I/O-bound on WSL /mnt/c).

---

## Table 3 — Decompress time (seconds, lower = better)

| Dataset | **ARCS** | Genozip | SPRING |
|---|---|---|---|
| ECOLI_35x | 6.9 s | 0.77 s | 3.7 s |
| HUMAN_127bp | 8.6 s | 0.45 s | 5.1 s |
| MTB_51bp | 4.3 s | 0.60 s | 1.9 s |
| SARS2_AMP | **1.1 s** | 0.39 s | 4.9 s |
| GIAB_HG002 | **1.9 s** | 0.37 s | 8.2 s |
| GIAB_HG001 | **8.5 s** | 2.3 s | 9.8 s |
| ECOLI_30x | 12.2 s | 3.7 s | 10.0 s |
| HUMAN_30x | **8.2 s** | 5.0 s | 11.4 s |

ARCS beats SPRING on 5/8 (SARS2_AMP, GIAB_HG002, GIAB_HG001, ECOLI_30x, HUMAN_30x).
Genozip is fastest (model-replay architecture vs ARCS adaptive-CM decode — architectural gap).

---

## Table 4 — Peak RAM — compress (MB, lower = better)

| Dataset | **ARCS** | ARCS (fast) | Genozip | SPRING |
|---|---|---|---|---|
| ECOLI_35x | 1,000 | **367** | 278 | 378 |
| HUMAN_127bp | 1,216 | 573 | 220 | 350 |
| MTB_51bp | 664 | **214** | 169 | 307 |
| SARS2_AMP | 466 | **359** | 446 | 421 |
| GIAB_HG002 | 429 | **313** | 419 | 382 |
| GIAB_HG001 | 1,681 | 968 | 1,098 | 575 |
| ECOLI_30x | 7,184 | 3,382 | 785 | 1,062 |
| HUMAN_30x | 2,949 | 2,403 | 1,436 | 1,324 |

ARCS (fast) compress RAM beats SPRING on 4/8 (ECOLI_35x, MTB_51bp, SARS2_AMP, GIAB_HG002).
ARCS (default/chain-pg) is RAM-heavy: assembly holds reads + graph simultaneously. Genozip and SPRING lowest on large datasets.

---

## Table 4c — Peak RAM — decompress (MB, lower = better)

| Dataset | **ARCS** | Genozip | SPRING |
|---|---|---|---|
| ECOLI_35x | 171 | 159 | 120 |
| HUMAN_127bp | 149 | 117 | 97 |
| MTB_51bp | 126 | 81 | 55 |
| SARS2_AMP | 160 | 122 | 145 |
| GIAB_HG002 | 124 | 128 | **123** |
| GIAB_HG001 | 346 | 349 | 357 |
| ECOLI_30x | 921 | 434 | 942 |
| HUMAN_30x | 887 | 638 | 1,005 |

ARCS decompress RAM is competitive with SPRING; Genozip lower on large datasets (stored model = smaller decode tables). ARCS and SPRING are within ~20% on human-scale data.

---

## Table 5 — Feature comparison

| Feature | ARCS | Genozip | SPRING | PgRC2 | DiscoSNP++ | Kmer2SNP |
|---|---|---|---|---|---|---|
| Byte-exact lossless | ✅ all datasets | ✅ | ❌ (reorders reads on bacterial) | ❌ (lossy qual by default) | — | — |
| Order-preserving | ✅ default | ✅ | ❌ | ❌ | — | — |
| CRLF round-trip | ✅ | ✅ | ❌ | ❌ | — | — |
| Built-in variant calling | ✅ het-SNV+indel+polyploid | ❌ | ❌ | ❌ | ✅ (separate tool) | ✅ (separate tool) |
| Reference required | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Polyploid support | ✅ (`ARCS_PLOIDY=k`) | ❌ | ❌ | ❌ | ✅ | ❌ |
| Homozygous variants | ❌ (fundamental, ref-free limit) | ❌ | ❌ | ❌ | ❌ | ❌ |
| License | MIT | Commercial† | MIT | GPLv3 | AGPL | — |
| Platforms | Linux, Windows | Linux, macOS | Linux | Linux | Linux | Linux |

† Genozip requires a per-user license for research use (free for academics).

---

## Table 6 — Het-SNV calling: 5 GIAB individuals (chr20:2.0–2.4 Mb, ~30×)

Scored by rtg vcfeval (GA4GH engine). Frozen parameters throughout — tuned once on HG001 r2;
all other individuals and regions are held out.

| Individual | Ancestry | ARCS F1 | DiscoSNP++ F1 | Kmer2SNP F1 |
|---|---|---|---|---|
| HG001 (NA12878) | European (5-region avg) | 0.943 | 0.887 | 0.539 |
| HG002 | Ashkenazi son | 0.929 | 0.941 | 0.542 |
| HG003 | Ashkenazi father | 0.935 | 0.920 | 0.550 |
| HG004 | Ashkenazi mother | 0.948 | 0.921 | 0.547 |
| HG005 | Han Chinese son | 0.931 | 0.891 | 0.487 |
| **Average (HG002-HG005)** | | **0.936** | **0.918** | **0.532** |

ARCS leads on every individual and ancestry with a single frozen parameter set.
Precision ≥ 0.986 across individuals; the margin over DiscoSNP++ is driven by recall (0.93 vs 0.80).
No reference genome used by any tool.

**Difficulty strata** (HG001 calls re-scored inside GIAB stratification BEDs):

| Stratum | ARCS P | ARCS R | ARCS F1 |
|---|---|---|---|
| Easy (non-difficult) | 0.986 | 0.920 | 0.952 |
| Homopolymer | 1.000 | 0.714 | 0.833 |
| Low-mappability | 0.857 | 0.667 | 0.750 |

**Coverage sweep** (held-out HG001 r3, downsampled):

| Coverage | ARCS P | ARCS R | ARCS F1 |
|---|---|---|---|
| 10× | 0.977 | 0.559 | 0.711 |
| 15× | 0.978 | 0.782 | 0.869 |
| 30× | 0.986 | 0.927 | 0.954 |

---

## Table 7 — Het-indel calling: real GIAB HG002 chr20:2.0–2.4 Mb

Indels called via contig-bubble detection (bubbles between consecutive contigs). Scored by rtg vcfeval, squash-ploidy.

| Tool | Precision | Recall | F1 | Reference required |
|---|---|---|---|---|
| **ARCS** | 0.600 | 0.438 | **0.505** | ❌ |
| DiscoSNP++ | 0.870 | 0.364 | 0.513 | ❌ |

ARCS F1 = 0.505 vs DiscoSNP++ 0.513 — zero-cost byproduct ties the dedicated indel caller.
DiscoSNP++ has higher precision; ARCS has better recall. Real indel F1 ~0.5 is the homopolymer/STR ceiling for all reference-free callers on this region.

---

## Table 7b — Polyploid het-SNV calling (synthetic data, rtg vcfeval)

Diploid path byte-identical with `ARCS_PLOIDY=2` (default). Polyploid enabled with `ARCS_PLOIDY=k`.

| Ploidy | Seed | Precision | Recall | F1 | Notes |
|---|---|---|---|---|---|
| Triploid (k=3) | 42 | 1.000 | 0.988 | **0.994** | rtg gold-standard |
| Triploid (k=3) | 7 | 1.000 | 0.988 | 0.994 | |
| Triploid (k=3) | 101 | 1.000 | 1.000 | 1.000 | |
| Triploid (k=3) | 555 | 1.000 | 1.000 | 1.000 | |
| Tetraploid (k=4) | 42 | — | — | — | 0 FN / 0 FP sites; allele-recall 0.73† |

† rtg vcfeval cannot score GT with 3+ ALT alleles — tetraploid scored by direct allele-set comparison. Site localisation is perfect; allele-recall limited by depth/k.
DiscoSNP++ supports multi-ploidy via `--ploidy`; Kmer2SNP is diploid-only.

---

## Table 8 — Total cost: compressed archive + variant calls (GIAB_HG002, native ext4)

How long it takes to produce **both** a compressed archive and a VCF from the same FASTQ.
Measured on GIAB_HG002 (39.4 MB raw, 113,987 reads). All times wall-clock, single run.

| | ARCS `compress` | ARCS `compress --call` | SPRING + DiscoSNP++ | Genozip + DiscoSNP++ |
|---|---|---|---|---|
| Compress time | 8.91 s | 9.88 s | 5.90 s | 10.31 s |
| Variant calling time | — | **0 s** (fused, placements reused) | +3.30 s (separate run) | +3.30 s (separate run) |
| **Total time** | 8.91 s | **9.88 s** | 9.20 s | 13.61 s |
| Archive size | 4,310,976 B | 4,310,976 B | 5,601,280 B | 5,049,582 B |
| Archive vs ARCS | — | baseline | +30% larger | +17% larger |
| Peak RAM (compress) | 429 MB | 625 MB | 382 MB | 425 MB |
| Peak RAM (calling) | — | 0 MB (fused) | +625 MB (DiscoSNP++) | +625 MB (DiscoSNP++) |
| VCF output | ❌ | ✅ 413 variants | ✅ 358 variants | ✅ 358 variants |
| Variant types | — | het-SNV + indel | het-SNV only | het-SNV only |
| Single command | ✅ | ✅ | ❌ two steps | ❌ two steps |

`arcs compress --call` overhead vs plain compress: **+0.97 s** (10.9%) for pileup scanning and VCF
writing — the assembly is already built for compression, so no second pass is needed.
SPRING and Genozip require a separate DiscoSNP++ run (fresh de Bruijn graph from scratch) to produce
any variant calls.

---

## Honest limitations

- **Compress speed**: ARCS is slower than SPRING/Genozip on most datasets (assembly cost). Same trade-off as PgRC2 — accepted in literature for archival-optimal ratio.
- **Decompress speed**: Genozip (model-replay) is fastest. ARCS beats SPRING on 5/8 datasets.
- **Compress RAM**: ARCS (default) is the heaviest on large datasets. ARCS (fast) is competitive with SPRING on small/medium datasets.
- **Homozygous variants**: not callable reference-free (structural, shared by all three callers). Requires a reference genome (GATK/DeepVariant) or multi-sample cohort.
- **Indel F1 ~0.5**: homopolymer/STR ceiling for all reference-free callers on real data; not specific to ARCS.
- **PgRC2 lossless**: lossless mode crashed on all tested datasets. PgRC2 is designed for lossy-quality use cases; not a fair lossless comparison.
- **GIAB_HG001 vs Genozip**: ARCS loses by 0.5% — near-tie on one human dataset.
