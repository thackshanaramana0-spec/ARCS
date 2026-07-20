# ARCS — paper thesis & contribution map

Date: 2026-07-16. The consolidation target after the quality-frontier exploration showed the 3-mer +
interpolation lever is exhausted (~2.20 bpq) and the easy next signals (read-level, tile, spatial) are
dead on the test data.

## The one big idea (thesis)
> A self-assembled **pseudogenome** is not merely a sequence-compression device (as in PgRC) — it is a
> **shared latent structure** that conditions BOTH streams of a FASTQ file. The same assembled
> consensus that lets us store reads as edits also predicts their quality scores. This **unifies** two
> problems every prior tool solves separately, and it is only possible reference-free because the reads
> are aligned to a structure built from the reads themselves.

One structure → two arrows: (1) sequence residual coding, (2) quality conditioning.

## Why this is novel (and where we are honest about borrowing)

**Borrowed / standard (state clearly — do not claim):**
- Greedy multi-contig assembly + majority-vote consensus (shape shared with NanoSpring, different
  regime: short Illumina vs nanopore).
- Context-mixing DNA entropy coding (GeCo3 / PAQ / cmix lineage).
- PPM-style interpolation / backoff (classic; Jelinek–Mercer / Witten–Bell).
- Adaptive arithmetic (range) coding.

**Genuinely novel — the paper's claims:**
1. **Unification.** One self-assembled pseudogenome conditions BOTH sequence and quality. No prior tool
   does this: PgRC uses the pseudogenome for sequence only; fqzcomp/ACO compress quality with no
   sequence model; CRAM/MPEG-G condition quality on an EXTERNAL reference (needs an aligner).
2. **Reference-free, sequence-motif-conditioned, LOSSLESS quality.** We predict Q_i from the local DNA
   motif — the physical cause of base-call difficulty — with zero stored side information (the sequence
   decodes before quality, so the decoder recomputes the context for free). fqzcomp/ACO do not use the
   base sequence; Quartz/GeneCodeq use it but LOSSILY; CRAM needs an external reference and uses the
   column MEAN, not the motif.
3. **Empirical demonstration** (held-out) that the local sequence motif carries ~10.6% of quality
   entropy beyond the quality-history model everyone uses, and that PPM interpolation realizes it in a
   real lossless codec.

## Headline results (DS7-100k, human chr22-class Illumina, lossless, dependency-free)
| track | ARCS | SOTA | delta |
|---|---|---|---|
| sequence (total_seq) | 353,065 B | PgRC2 360,025 B | **−1.9%** (5k: −14.2%) |
| quality | 2.213 bpq | fqzcomp context 2.419 bpq | **−8.5%** |
| quality vs general | 2.213 | xz 2.899 / bzip2 3.052 | −24% / −27% |
| total archive (quality lever) | 4,982,611 B | (was 5,695,934) | **−12.5%** |

All record-lossless (reorders reads, same contract as PgRC). Dependency-free: own DNA CM coder
(dna_coder), own quality coder (qual_cm), own assembler.

## Methodological rigor (a strength to foreground)
- **Held-out train/test** measurements throughout — not in-sample. We CORRECTED an in-sample "2.09
  floor" mirage to the true held-out floor 2.204, and hit it (2.213).
- **Negative results reported**: idea-B pg-surprise = null; read-level/tile/spatial frontiers = dead;
  binary-decomposition + logistic mixer = failed. Science, not cherry-picking.
- **keep-smaller gates**: every new coder is chosen only if it beats the incumbent → provably no
  regression on any dataset.

## Honest limitations (must state, and ideally close before submission)
1. **One instrument/organism tested.** MUST validate on: binned NovaSeq (4–8 quality levels — does
   sequence-conditioning still win when the alphabet is tiny?), other organisms, longer/shorter reads.
   This is the single most important missing experiment.
2. **fqzcomp not run head-to-head.** We compared against its context's measured entropy (fair, slightly
   generous to it). Should install and run the actual fqzcomp5 binary for a direct number.
3. **Reorders reads** (record-lossless). If original order is required, costs a permutation.
4. **pg-consensus quality lever unexplored** (predicted marginal <1–2%). Optional.

## Suggested paper structure
1. Intro: FASTQ = sequence + quality; prior tools solve them SEPARATELY; the missed opportunity.
2. The shared-pseudogenome thesis (the unifying figure).
3. Sequence: multi-contig consensus assembler → beats PgRC2.
4. Quality: sequence-motif conditioning + PPM interpolation → beats fqzcomp; held-out analysis;
   the −10.6% sequence-entropy result.
5. Unified dependency-free codec; losslessness argument (RC-involutive; deterministic integer coders).
6. Results (multi-dataset — TO COMPLETE), ablations, negative results.
7. Limitations, reproducibility.

## The single most valuable next action for the paper
Run the multi-dataset validation (esp. binned NovaSeq) + a real fqzcomp5 head-to-head. If
sequence-conditioning holds across datasets, this is a strong, honest, novel paper. If it only helps on
full-range quality, the claim narrows but stays valid — and we report that honestly.
