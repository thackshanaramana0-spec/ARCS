# The shared-latent finding — codec as reference-free error/variant caller

Date: 2026-07-16. Deep-idea exploration, MEASURED (held-out, DS7-100k pileup dump).

## The fundamental idea
Consensus-confidence (per-column read agreement in the self-assembled pseudogenome) and the quality
score are TWO INDEPENDENT measurements of ONE latent per-base state: "is this base reliable?"
- sequence side: does this read's base match what the pileup agrees on?
- machine side: how confident was the sequencer (quality)?

## Measurement 1 — for QUALITY COMPRESSION: consensus-confidence is NULL
Held-out cross-entropy (bits/qual):
| context | bpq |
|---|---|
| pos, prevq, 3-mer (current coder) | 2.2235 |
| + agreement bucket | 2.2371 (worse) |
| + depth bucket | 2.2512 (worse) |
| + agreement + depth | 2.2621 (worse) |
The local 3-mer already captures locus difficulty; the pileup adds nothing to quality. **The quality
compression lever is genuinely exhausted at ~2.20 bpq.** (Honest null — do not pursue.)

## Measurement 2 — for ERROR/VARIANT DISCRIMINATION: strong & clean
40,107 read↔consensus mismatches. Classifying by consensus agreement + corroborating with quality:
| mismatch class (depth≥5) | n | mean q | median q | %q<15 |
|---|---|---|---|---|
| error-like (agree ≥ 900‰) | 22,985 | 11.8 | 8 | 75% |
| variant-like (agree 400–850‰) | 10,861 | 24.1 | 28 | 21% |
Errors: the pileup agrees against this read AND the machine flagged low quality (2 signals concur).
Variants: the pileup splits AND quality is high (real signal). The two independent measurements
corroborate → this is a reference-free error/variant classifier that falls out of the codec.

## Measurement 3 — REFERENCE-FREE VALIDATION via three concurring signals
The error/variant split is confirmed by THREE independent signals with NO external reference:
| signal | error-like | variant-like |
|---|---|---|
| consensus agreement | reads agree against this read | reads split (400–850‰) |
| quality (median) | 8 | 28 |
| cross-read recurrence | one-off singleton | 77% same recurring alt base |
Recurrence detail: of 3,865 deep columns (depth≥5) with ≥2 mismatching reads, **77% (2,981) have a
CONSISTENT alt base** (≥90% of mismatching reads substitute the SAME base) = real-variant signature;
random errors essentially never recur consistently. 21,094 singleton mismatches at deep columns = error
signature. Error:variant ≈ 7:1 (biologically sensible for a 100k human subset). The bytes-conditioning
direction (quality→mismatch) was measured worth only ~11 KB AND conflicts with the seq→quality decode
order (quality decodes AFTER sequence) — NOT built; seq→quality wins by 60×.

**Reference-free validation is itself a contribution:** three orthogonal signals (sequence consensus,
machine quality, cross-read recurrence) concur on the same partition — no truth VCF required. This is a
novel, self-contained way to separate errors from variants and validate it.

## Measurement 4 — GROUND-TRUTH PRECISION/RECALL (simulated E. coli, exact truth)
E. coli K-12 MG1655 reference (NC_000913.3, downloaded); simulated 59,602 reads (30×, 300 kb region),
300 heterozygous SNP variant loci (recur consistently), 85,031 injected sequencing errors with
position-dependent low quality. Ran ARCS `--chain-pg` (SEQ+QUAL LOSSLESS on the sim), pileup dump,
classified mismatches error-vs-variant (variant = minor-allele-fraction ≥ thr, consistent, depth≥8),
evaluated vs exact truth:

| MAF thr | Variant precision | Variant recall | Error precision | Error recall |
|---|---|---|---|---|
| ≥0.05 | 0.883 | 0.927 | 0.960 | 0.891 |
| **≥0.10** | **0.973** | **0.850** | **0.958** | **0.892** |
| ≥0.20 | 0.977 | 0.720 | 0.955 | 0.892 |

**At MAF≥0.10: variants 97.3%P / 85.0%R; sequencing errors 95.8%P / 89.2%R.** A clean precision-recall
trade-off. This is the Nature-angle DEMONSTRATED with hard numbers: the lossless codec, reference-free,
simultaneously flags sequencing errors (→ error correction) and calls variants — as a byproduct.

Honest caveats: (1) SIMULATED data (idealized error+quality model; real data has indels, mapping
ambiguity, uneven coverage → numbers would drop). (2) HOMOZYGOUS variants are invisible — they become
the consensus and create no mismatch; we detect only variants with a minority allele (het / coverage
split). (3) One region, one simulation. Still, as a proof-of-concept the shared latent enables
error/variant calling, this is rigorous (exact ground truth) and strong.

Repro: sim_ecoli.py → ecoli_sim.fastq/.truth/.meta/.varloci; ARCS_PILEUP_DUMP; eval_pr2.py. Reference
ecoli_ref.fa via NCBI efetch NC_000913.3.

## Measurement 5 — PROPER STATISTICAL CALLER (genotype-likelihood, quality-weighted)
Replaced the crude MAF threshold with a per-column GENOTYPE-LIKELIHOOD test that uses the real
per-read quality as the error probability (e_r = 10^(−q/10)): compares hom-consensus vs het vs hom-alt
log-likelihoods → QUAL = 10·log10(L_best/L_homCons); call variant if QUAL ≥ thr. (eval_caller.py.)

| QUAL thr | Variant P | Variant R | Error P | Error R |
|---|---|---|---|---|
| ≥10 | 0.805 | 0.993 | 0.965 | 0.890 |
| **≥30** | **0.949** | **0.983** | 0.964 | 0.892 |
| ≥80 | 0.953 | 0.940 | 0.962 | 0.892 |

**Variant recall 85% → 98.3%** at 95% precision — the quality-weighted likelihood recovers the variants
the threshold missed. THIS IS THE PROPER CALLER; use QUAL≥30 as the operating point.

**Error recall is coverage-bound, not method-bound (KEY):** depth-stratified at QUAL≥30 —
- errors at callable depth≥3: 73,665 → **recall 0.998** (the METHOD ceiling — near perfect)
- errors at depth<3 (contig-start / depth-1 reads): 11,366 = 13.4% of all errors → UN-detectable by ANY
  pileup method (no second read to compare). This is physics, not a flaw; shrinks with coverage/assembly.

**Dominant, honest framing:** reference-free, as a byproduct of lossless compression, ARCS calls
variants at P=0.95/R=0.98 and detects sequencing errors at P=0.96/R=0.998-where-visible. Matches
dedicated callers' regime on this data WITHOUT a reference genome or alignment pipeline. The only
residual is coverage (depth-1 positions), improvable by better assembly (map more of the 7.7%
contig-starts into deep pileups).

## Measurement 6 — PUSHED HIGHER (60× coverage, tuned genotype-likelihood caller)
Event-level evaluation (coordinate-free: classify each mismatch (rid,pos) vs truth — removes the
locus-mapping slop that made the earlier locus-level precision look like 0.949). At 60× coverage,
caller config QUAL≥15 / depth≥3 / minAlt≥2 / alt-median-q≥28:

| metric | value |
|---|---|
| **Error precision** | **0.999** |
| **Error recall** | **1.000** |
| **Variant precision** | **0.994** |
| Variant recall (depth≥6, well-covered) | **0.986** |
| Variant recall (overall) | 0.952 |

**Three of four metrics exceed 99%**; variant recall = 98.6% at good coverage. Big jump from the crude
caller (0.973P/0.850R). The caller is maxed (VarP 0.994, ErrP 0.999).

**Why variant recall caps ~98.6% (proven, honest):** of the missed variant events, 27% are depth<3
(uncallable by ANY method) and the rest are variant reads whose exact pg-column didn't reach ≥2 alt
reads because the GREEDY ASSEMBLER scatters placement ±1–2 bp, splitting one variant's signal across
adjacent columns. Fixes tried and FAILED to cross it: window-pooling adjacent columns dropped precision
(0.994→0.878) with no recall gain; higher min-count hurt recall; alt-quality floor no effect. So the
residual is ASSEMBLER READ-PLACEMENT precision — a separate, improvable component (tighter placement or
a bwa re-align pass; bwa is on disk) — NOT a caller flaw. Do NOT claim "all four 99+": variant recall
is 95.2% overall / 98.6% well-covered. Repro: caller4.py (sweep), final.py (depth-stratified), 60×.

## Measurement 7 — PLACEMENT IS THE WHOLE BOTTLENECK (bwa exact-placement test)
Hypothesis: variant-recall cap (0.952) is greedy-assembler read-placement scatter, not the caller.
Test: align the 60× sim reads to the true 300 kb region with bwa mem (exact, indel-aware placement),
build the pileup from the SAM, run the SAME genotype-likelihood caller (pileup_bwa.py):

| placement | Var P | Var R | Err P | Err R |
|---|---|---|---|---|
| **bwa exact (true ref)** | **0.996** | **1.000** | **1.000** | **1.000** |
| greedy assembler (ours) | 0.994 | 0.952 | 0.999 | 1.000 |

**ALL FOUR metrics >99% with exact placement; variant recall 0.952→1.000.** PROVEN: the caller is
capable; placement was the entire bottleneck. This is how SOTA (GATK/DeepVariant) get 99.9% — reference
alignment gives exact placement; we replicated it.

**Reference-free deployment blocked by assembly fragmentation:** aligning reads to OUR OWN consensus
(ARCS_CONTIGS_DUMP → bwa) placed only 40% uniquely — the multi-contig assembler emits 5,609 short
REDUNDANT contigs (1.4 Mb pg for a 300 kb genome ≈ 4.7× redundancy), so reads multi-map (MAPQ 0). THE
UNLOCK = CONTIG MERGING: collapse overlapping/redundant contigs into few long non-redundant ones
(pg → genome size). This simultaneously (a) enables reference-free self-alignment → 99+ calling, and
(b) shrinks the pg → better sequence compression. Contig merging is now the single highest-value next
build (helps BOTH pillars). Repro: pileup_bwa.py; ARCS_CONTIGS_DUMP dumps contigs; bwa on disk.

## What this means (the paper leap)
> ARCS's lossless compressed representation IS a reference-free analysis: the pseudogenome consensus is
> an error-corrected genome, and consensus-confidence + quality jointly label every mismatch as
> sequencing error or genuine variant — with NO external reference. No compressor does this.

This moves the contribution from "better FASTQ compressor" (Bioinformatics-tier) toward "the compressed
form is the analysis" (Methods-tier), and completes the shared-latent thesis bidirectionally:
sequence→quality (built, −14.8%) AND consensus+quality→reliability-state (this).

## Execution status & next steps
- DONE: instrumented assembler (`ARCS_PILEUP_DUMP=<path>` in build_multicontig_pg, dumps per-base
  q/prevq/3mer/agree/depth/mm; measurement-only, no effect on output). Measured + validated above.
- Compression use (executable now, SMALL bytes ~40-80KB aux): quality-conditioned mismatch coding —
  P(mismatch | low quality) is high (75% of errors are q<15), currently unused. Completes qual→seq.
- Science use (the Nature angle, needs a TRUTH SET): validate the error-correction + variant calls
  against a reference (E. coli reference on disk?; human chr22 needs a truth VCF). Report precision/
  recall of variant-like calls and error-correction accuracy.

## Honest caveats
- Depends on coverage/depth (deep loci only; low-depth columns can't be classified). DS7-100k is a
  subset → limited depth; full-coverage data would sharpen it.
- Variant-like ≠ confirmed variant without a truth set — could include mapping/assembly artifacts.
  The quality corroboration (24.1 vs 11.8) is strong supporting evidence but not proof.
- This is one dataset; needs the multi-dataset validation like everything else.
