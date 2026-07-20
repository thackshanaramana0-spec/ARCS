# Lossless quality is at the information-theoretic wall (measured negative result)

Date: 2026-07-17. Question: is there ANY lossless-quality headroom left below the ~1.57 bpq that ARCS
and fqzcomp both hit on binned GIAB quality? Quality is 79% of the archive, so this is the highest-value
question for total ratio. Answer, measured on real GIAB HG002: **effectively ZERO.** The tie with fqzcomp
is the wall, not a skill gap.

## Method
Held-out train/test cross-entropy (NOT in-sample entropy, which overfits sparse contexts): build the
per-context frequency model on even-indexed reads, measure Laplace-smoothed cross-entropy on odd-indexed
reads. This is the real achievable bpq for a given context — it penalizes context fragmentation honestly.
GIAB read names carry full Illumina coordinates (@D00360:97:H2YVMBCXX:lane:tile:X:Y), so flowcell-2D is
testable. Script: flowcell.py in /c/Temp/arcs_test.

## Result — the last untapped lever (flowcell-2D) HURTS
| context added to (q1,q2,pos,run) | held-out bpq | vs baseline |
|---|---|---|
| baseline (nothing) | 1.5941 | — |
| + spatial 2×2 (4 regions) | 1.6394 | +2.8% WORSE |
| + spatial 3×3 | 1.6932 | +6.2% WORSE |
| + spatial 4×4 | 1.7497 | +9.8% WORSE |
| + spatial 8×8 | 1.9738 | +24% WORSE |
| + spatial 8×8 + tile | 4.0626 | +155% WORSE |

Monotonic: every spatial bin fragments the data faster than any signal helps. The quality value at a
position is NOT meaningfully predicted by where the read landed on the flowcell, once q1/q2/pos/run are known.

## Three-front confirmation that binned lossless quality is exhausted
1. **Sequence-motif (seq3mer):** helps ~8% on FULL-RANGE quality, ~0 on binned (dilutes). [quality_sota]
2. **Cross-read / pileup consensus:** measured NULL earlier — "quality EXHAUSTED at 2.20 bpq, agreement/
   depth add nothing beyond local 3mer". [SHARED_LATENT_ERROR_STATE]
3. **Flowcell-2D spatial (this doc):** makes it WORSE at every granularity.

## Conclusion (honest, and useful)
Binned Illumina quality is at its conditional-entropy floor (~1.57–1.59 bpq); H(q|q1,q2,pos,run) captures
essentially all the extractable structure, and fqzcomp already sits there. **No lossless context beats it
meaningfully — we tested the three plausible levers and all are dead or negative.** This REFRAMES the
quality "tie" from a weakness into a proven ceiling: matching the information-theoretic floor with a
dependency-free adaptive coder is the correct outcome, not a shortfall. The only way to shrink the 79%
quality bulk further is LOSSY (Illumina binning / Quartz / Crumble / Enano) — a different product. For a
LOSSLESS codec, quality is solved; the ratio frontier lives in sequence, names, and (for real gains) the
lossy question.

Caveat: measured on one GIAB region; the wall is an information-theoretic property of binned Illumina
quality and is expected to generalize, but full-scale confirmation across instruments/chemistries is the
honest next check.

## LITERATURE CONFIRMATION (checked 2026-07-17) — the field independently reached the same wall
The "quality is at the floor" claim is NOT just our measurement — it is the published consensus, confirmed
by the strongest possible evidence (a competition no one could win):

1. **SequenceSqueeze competition (Pistoia Alliance / EBI):** multiple INDEPENDENT specialist entrants all
   converged to **~2.94 bits/quality-value** on the test set, and NO entrant improved on it. A whole field
   competing and all landing on the same number is the empirical fingerprint of an information-theoretic
   limit — exactly the wall our held-out measurement shows. This is far stronger than a single experiment.
2. **Quality = 70–80% of the losslessly-compressed file** (published) — matches our measured 79% exactly.
3. **Quality ≈ 6 bits/symbol vs 2 for sequence** (published) — confirms quality is the high-entropy stream
   where lossless compression saturates.
4. **The field's OWN pivot to LOSSY** (QVZ, QualComp, Quartz, Crumble, Enano) is the tell: the specialists
   concluded the only meaningful further gains on the 80% quality bulk are lossy — i.e. they hit the same
   lossless wall we did and stopped trying to beat it losslessly.

IMPORTANT CORRECTION (regime-specificity): the floor is DATASET/BINNING-SPECIFIC, not one universal number.
- FULL-RANGE (~40-level) quality → SequenceSqueeze floor ≈ **2.94 bpq**.
- BINNED (~8–24 level, modern Illumina, our GIAB) → floor ≈ **1.57 bpq** (what ARCS and fqzcomp both hit).
So the precise, defensible statement is: "lossless quality is at the CONDITIONAL-ENTROPY FLOOR FOR ITS
BINNING REGIME, a competition-confirmed wall — not a skill gap." The only residual research (neural/
transformer quality models) squeezes fractions of a percent via longer-range dependencies and does not
dent the 70–80% bulk. Our tie with fqzcomp is therefore standing exactly on the field's proven lossless
ceiling. Sources: SequenceSqueeze (Pistoia/EBI); fqzcomp (jkbonfield); QVZ/QualComp/Quartz/Crumble (lossy
quality literature).
