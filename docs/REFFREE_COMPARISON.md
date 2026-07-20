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
