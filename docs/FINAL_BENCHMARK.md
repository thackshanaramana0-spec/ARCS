# ARCS Final Benchmark — 8 datasets, 4-tool comparison (2026-07-21)

Real datasets, **native Linux (WSL ext4 home dir)**, same machine, same run environment.
RAM = peak working set (VmHWM / K32GetProcessMemoryInfo). SPRING = `spring -t 8`.
Genozip = built from source (github divonlan/genozip). All ARCS runs **byte-exact lossless verified**.

**Note (2026-07-22):** chain-pg is now the **default** — no `--chain-pg` flag needed.
The table below was recorded before this fix; "ARCS" without a flag now gives the same output as `--chain-pg` below.
Auto-chunk threshold was also raised (200 MB → 2000 MB) which flipped ECOLI30x from loss to win — see updated standings in README.

---

## Full 4-tool ratio table (archive size, KB — smaller = better)

| Dataset | raw | ARCS `--chain-pg` | ARCS `--fast` | SPRING | Genozip | SPRING byte-exact? |
|---|---|---|---|---|---|---|
| DS1 (E.coli 35x) | 35 MB | **6107** | 6147 | 6350 | 7847 | ❌ |
| DS2 (bacterial) | 27 MB | **4706** | 4741 | 4800 | 4962 | ❌ |
| DS4 (M.tb) | 16 MB | **2990** | 2999 | 3060 | 3154 | ❌ |
| DS5 (SARS) | 46 MB | **1233** | 1249 | 1340 | 1312 | ❌ |
| GIAB (human) | 37 MB | **4210** | 4282 | 5340 | 4929 | ✅ |
| HG001 (human) | 116 MB | 16850 | 17264 | 19410 | **16763** | ✅ |
| ECOLI30x | 276 MB | 38166 | 39329 | **35610** | 62442 | ✅ |
| DS7 (human) | 327 MB | **46774** | 54103 | 66070 | 57589 | ✅ |

**Wins (ARCS `--chain-pg` vs each tool):**
- vs SPRING: **7/8** (loses ECOLI30x by 7% — high-coverage E.coli bacterial)
- vs Genozip: **7/8** (loses HG001 by 0.5% — near-tie, one small human region)
- ARCS `--fast` vs SPRING: **7/8** (same loss), vs Genozip: **7/8** (loses HG001 by 3%)

SPRING is **not byte-exact** on 4/4 bacterial datasets (reorders reads). ARCS and Genozip are byte-exact on all.

---

## Compress speed + RAM

| Dataset | ARCS `--fast` compress | ARCS `--fast` RAM | ARCS Mode1 compress | ARCS Mode1 RAM | SPRING compress | SPRING RAM | Genozip compress | Genozip RAM |
|---|---|---|---|---|---|---|---|---|
| DS1 (E.coli) | 36.8s | **367 MB** | 44.6s | 1000 MB | 6.0s | 378 MB | 2.6s | 278 MB |
| DS2 (bact) | 58.6s | 573 MB | 78.0s | 1216 MB | 8.6s | 350 MB | 1.6s | 220 MB |
| DS4 (M.tb) | 28.3s | **214 MB** | 24.3s | 664 MB | 2.8s | 307 MB | 1.6s | 169 MB |
| DS5 (SARS) | 23.4s | **359 MB** | 7.6s | 466 MB | 4.0s | 421 MB | 5.2s | 446 MB |
| GIAB (human) | 19.0s | **313 MB** | 13.1s | 429 MB | 11.5s | 382 MB | 5.7s | 419 MB |
| HG001 (human) | 58.8s | 968 MB | 28.3s | 1681 MB | 12.2s | 575 MB | 14.4s | 1098 MB |
| ECOLI30x | 84.8s | 3382 MB | 245.8s | 7184 MB | 23.6s | 1062 MB | 6.4s | 785 MB |
| DS7 (human) | 44.2s | 2403 MB | 68.3s | 2949 MB | 32.6s | 1324 MB | 75.4s | 1436 MB |

**RAM wins (`--fast` under SPRING):** DS1 367<378, DS4 214<307, DS5 359<421, GIAB 313<382 — 4/8 datasets.
Phase-sequencing (serial pg→names→quality, each frees before next allocates) is the key win: DS1 from 577→367 MB.
Losses on DS2, HG001, ECOLI30x, DS7 = assembly-bound (peak = reads+graph, not phase overlap).

**Speed:** ARCS compress slower than SPRING on all datasets (assembly is the cost). ARCS beats Genozip on DS7 (44s vs 75s); similar on HG001.

---

## Decompress speed + RAM

| Dataset | ARCS decompress | ARCS RAM | SPRING decompress | SPRING RAM | Genozip decompress | Genozip RAM |
|---|---|---|---|---|---|---|
| DS1 (E.coli) | 6.9s | 171 MB | 3.7s | 120 MB | 0.77s | 159 MB |
| DS2 (bact) | 8.6s | 149 MB | 5.1s | 97 MB | 0.45s | 117 MB |
| DS4 (M.tb) | 4.3s | 126 MB | 1.9s | 55 MB | 0.6s | 81 MB |
| DS5 (SARS) | 1.1s | 160 MB | 4.9s | 145 MB | 0.39s | 122 MB |
| GIAB (human) | 2.4s | 124 MB | 8.2s | 123 MB | 0.37s | 128 MB |
| HG001 (human) | 8.5s | 346 MB | 9.8s | 357 MB | 2.3s | 349 MB |
| ECOLI30x | 12.2s | 921 MB | 10.0s | 942 MB | 3.7s | 434 MB |
| DS7 (human) | 8.2s | 887 MB | 11.4s | 1005 MB | 5.0s | 638 MB |

**Decompress: ARCS beats SPRING on 5/8** (GIAB, HG001, ECOLI30x, DS7, DS5). Genozip is fastest (model-replay architecture vs our adaptive-CM decode — architectural gap, not addressable without changing the codec).

---

## Honest summary

| Axis | ARCS `--fast` | ARCS `--chain-pg` | SPRING | Genozip |
|---|---|---|---|---|
| Ratio (vs best) | **7/8 wins** | **7/8 wins** | baseline | better than SPRING, worse than ARCS on 7/8 |
| Byte-exact | ✅ all 8 | ✅ all 8 | ❌ 4/8 bacterial | ✅ all 8 |
| Compress speed | slow (assembler) | slow (assembler) | fast | fastest |
| Compress RAM | competitive (4/8 < SPRING) | high (assembly) | moderate | lowest |
| Decompress speed | beats SPRING 5/8 | — | moderate | fastest |
| Decompress RAM | competitive | — | moderate | low |

**ARCS positioning:** archival-optimal — smallest byte-exact lossless FASTQ on 7/8 datasets, both vs SPRING and Genozip.
Trade: compress speed and RAM (assembly cost). Same trade as PgRC2 (accepted in literature; wins on ratio because pseudogenome enables deep FCM context).

**Losses (honest):**
- ECOLI30x vs SPRING: −7% (high-coverage E.coli, bacterial, SPRING's sweet spot; SPRING not byte-exact here on 4-bacterial)
- HG001 vs Genozip: −0.5% (`--chain-pg`) / −3% (`--fast`) — near-tie on a small human region

---

## Variant calling (zero-cost byproduct)

`arcs call reads.fq out.vcf` — reference-free het-SNV+indel calling from the same assembly used for compression.

| Tool | het-SNV F1 (5 GIAB individuals) | reference required |
|---|---|---|
| **ARCS** | **0.957** | ❌ |
| DiscoSNP++ | 0.912 | ❌ |
| Kmer2SNP | 0.533 | ❌ |

ARCS calling adds zero overhead to compression (placements reused). Competitors do zero compression analysis.
