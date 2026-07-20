# Reference-free heterozygous-SNV calling — head-to-head vs DiscoSNP++ and Kmer2SNP

All three tools scored by the **same** gold-standard pipeline: `rtg vcfeval` 3.12.1
(the GA4GH engine `hap.py` uses internally), allele-aware, `--squash-ploidy`,
restricted to the GIAB confident regions, no arbitrary matching window. Each tool's
calls are placed in genome coordinates by the same means (BWA), then compared to the
GIAB HG001 truth set. Five chr20 regions (~400 kb each, ~30×), including one held out.

ARCS calls are a **byproduct of the lossless compression assembly** — no external
aligner or k-mer tool (`arcs call reads.fq out.vcf`). DiscoSNP++ and Kmer2SNP are
dedicated reference-free callers, each given its best reasonable configuration
(DiscoSNP++ default `-c 3`; Kmer2SNP coverage band chosen by sweep).

## Heterozygous-SNV F1 (rtg vcfeval)

| region            | ARCS  | DiscoSNP++ | Kmer2SNP |
|-------------------|:-----:|:----------:|:--------:|
| r2 (chr20:3.0–3.4M)| 0.923 |   0.826    |  0.550   |
| r3 (4.0–4.4M, held-out)| **0.954** | 0.886 | 0.544 |
| na (2.0–2.4M)     | 0.964 |   0.899    |  0.528   |
| r4 (5.0–5.4M)     | 0.939 |   0.884    |  0.503   |
| r5 (6.0–6.4M)     | 0.934 |   0.942    |  0.570   |
| **average**       |**0.943**|  0.887   |  0.539   |

ARCS leads on 4 of 5 regions (DiscoSNP++ edges r5 by 0.008). The advantage is in
**recall**: all three tools reach high precision (0.95–0.99), but ARCS recovers more
true heterozygous sites (recall ≈ 0.93 vs DiscoSNP++ ≈ 0.80 vs Kmer2SNP ≈ 0.38) —
the assembly-based consensus positions and supports more variants than short-k-mer
bubble detection.

## Cross-individual generalization (frozen parameters, held-out sample)

The five regions above are all HG001 (NA12878). To test whether the result is
sample-specific or overfit, we ran the **same caller with identical frozen
parameters** — no re-tuning — on a **different individual, HG002** (Ashkenazi son),
on a **different reference build (GRCh38)**: chr20:2.0–2.4 Mb, Illumina ~30×.

| tool        | HG001 (5-region avg) | HG002 (held-out individual) |
|-------------|:--------------------:|:---------------------------:|
| **ARCS**    | **0.943**            | **0.971** (P .988 R .955)   |
| DiscoSNP++  | 0.887                | 0.941                       |
| Kmer2SNP    | 0.539                | 0.542                       |

The ranking is **stable across individuals** (ARCS > DiscoSNP++ > Kmer2SNP), and
ARCS's HG002 F1 (0.971) is *higher* than its HG001 average — with **zero parameter
changes**. This is direct evidence the ~0.94 is **not hypertuned to one sample**: it
holds, and leads, on a held-out individual and a different genome build. (Tuning was
done once on HG001 r2; every other region and the entire HG002 sample are held out.)

## Generalization matrix (frozen parameters throughout)

Beyond the 5 HG001 regions, the frozen caller was tested across four additional
axes to check the ~0.94 is a robust property, not a single-condition artifact:

**Cross-individual** (held-out sample, GRCh38) — ARCS 0.971 vs DiscoSNP++ 0.941 vs
Kmer2SNP 0.542. Ranking stable, ARCS leads (see above).

**Difficulty strata** (HG001 calls re-scored inside GIAB stratification BEDs):

| stratum          | ARCS P | ARCS R | ARCS F1 |
|------------------|:------:|:------:|:-------:|
| easy (non-difficult) | 0.986 | 0.920 | 0.952 |
| homopolymer      | 1.000  | 0.714  | 0.833   |
| low-mappability  | 0.857  | 0.667  | 0.750   |

Precision holds high (0.86–1.0) even in hard regions; recall falls in
homopolymer/low-map — the expected profile (hard for every caller).

**Coverage sweep** (held-out r3, downsampled):

| coverage | ARCS P | ARCS R | ARCS F1 |
|----------|:------:|:------:|:-------:|
| 10×      | 0.977  | 0.559  | 0.711   |
| 15×      | 0.978  | 0.782  | 0.869   |
| 30×      | 0.986  | 0.927  | 0.954   |

Precision is coverage-stable (~0.98 at all depths); recall scales with depth —
graceful degradation, no breakage at low coverage.

**Other chromosome** (HG001 chr21:30.0–30.4 Mb): ARCS 0.915 vs DiscoSNP++ 0.930 vs
Kmer2SNP 0.613 — a region where ARCS and DiscoSNP++ are neck-and-neck (DiscoSNP++
edges by 0.015), reported honestly.

**Summary:** across individual, difficulty, coverage, and chromosome, ARCS holds
F1 ≈ 0.90–0.97 het-SNV with high stable precision, leading or tying the best
reference-free competitor on every axis except one chr21 region — all with a single
frozen parameter set tuned once on HG001 r2.

## Scope and honesty

- **This is the heterozygous-SNV task**, the established reference-free niche. On
  **all** SNVs (het + homozygous) ARCS averages F1 ≈ 0.74; the gap is homozygous
  variants, which are **structurally invisible** to reference-free bubble calling (a
  homozygous variant forms no bubble in a self-assembled consensus) — a limitation
  shared by all three tools, not specific to ARCS.
- Reference-free calling **complements, not replaces**, reference-based pipelines
  (GATK/DeepVariant, ~0.99 all-variant) that are required for homozygous variants.
- The BWA step maps calls to reference coordinates **for evaluation only**; in true
  reference-free use (novel organism) calls remain in contig coordinates.

## Reproduce

```bash
arcs call reads.fastq calls.vcf          # ARCS (byproduct of the archive)
scripts/lift_vcf.py ... | rtg vcfeval ... # score vs GIAB truth (see scripts/)
```
Validation drivers: `scripts/lift_vcf.py`, `scripts/kmer2snp_to_vcf.py`,
`scripts/eval_vcf.py`.
