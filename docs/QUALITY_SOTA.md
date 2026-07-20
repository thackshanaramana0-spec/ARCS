# SOTA of lossless quality-score compression + our design

Survey date: 2026-07-16. Goal: make ARCS quality genuinely best-in-class (quality ≈ 80% of archive).

## ★★★ FINAL: INTERPOLATED SEQUENCE-CONDITIONED QUALITY — 2.213 bpq, ~8% below fqzcomp (2026-07-16)

Two independent levers, both principled (measured, not tuned):
1. **Sequence-conditioning** — predict quality from the local DNA motif (the physical cause of
   base-call difficulty), which only we have via pseudogenome-aligned reads (seq decodes before
   quality → free, lossless). Detail in the section below.
2. **PPM-style interpolation** — code each symbol against λ·P_child + (1−λ)·P_parent where the child is
   {q1,q2,pos}⊗3-mer and the parent is the dense {q1,q2,pos}, λ = obs_child/(obs_child+K). Per-symbol
   backoff removes the sparse child's cold-start cost (integer, deterministic → lossless). Held-out
   measurement said this reaches 2.204 bpq vs 2.287 for a single 3-mer model.

**FINAL RESULT (default, lossless, 7/7 ctests):**
| coder (DS7-100k quality) | bytes | bpq | Δ vs static |
|---|---|---|---|
| static-rANS (old) | 4,807,480 | 2.599 | — |
| adaptive Q-ary, quality-history | 4,711,835 | 2.547 | −1.99% |
| + sequence context (hard child) | 4,460,085 | 2.411 | −7.23% |
| **+ interpolation (final)** | **4,094,157** | **2.213** | **−14.84%** |

- **vs fqzcomp-q2 floor 2.419 → we are at 2.213 = ~8.5% below fqzcomp**, real coded, lossless.
- Essentially at the held-out floor (2.204).
- **Total archive 5,695,934 → 4,982,611 = −12.5%** (quality is ~80% of it).
- Interpolation also FIXED the 5k case: quality 253,939 → 229,218 (−9.7%) — the backoff handles sparse
  data, so the coder now wins on BOTH 5k and 100k (the earlier hard-child version lost +17% on 5k).
  Default ON. Tunables in header: FREQ_INC (adapt speed), QCM_K (interpolation strength, flat 12..32),
  USE_SEQ. Disable seq via ARCS_QUAL_NOSEQ, CM via ARCS_QUAL_NOCM.

**Note on the 2.09 "floor":** that earlier number was IN-SAMPLE (overfit). The real HELD-OUT floor is
2.204, and wider-than-3-mer contexts do not beat it (measured). We are at 2.213 — at the floor.

**Next-frontier signals MEASURED (held-out, DS7-100k) — ALL DEAD on this data:**
| context added on top of pos,q1,q2,3mer (base 2.326) | bpq | verdict |
|---|---|---|
| + read-global mean quality | 2.341 | worse (autocorrelation already has it) |
| + tile id | 2.953 | much worse (no signal, just sparsity) |
| + x/y spatial bin | 2.949 | much worse |
Read-level, tile, and naive spatial add NOTHING (ACO's real 2D gain needs physical-neighbor quality via
read reordering — big architecture change, contested ~5%). The only un-measured lever is the
**pg-consensus motif** (condition on the assembled consensus instead of per-read bases) — likely
marginal (<1-2%: reads ≈ consensus except at ~1% error positions, already flagged by is_dev) but it is
the concrete expression of the paper thesis. VERDICT: this lever is exhausted at ~2.20; further quality
gains need a fundamentally different signal, not more tuning. Pivot to consolidation/paper.

---

## ★★ GAME-CHANGER: SEQUENCE-CONDITIONED QUALITY — beats fqzcomp's context (2026-07-16)

**The fundamental insight (not fine-tuning, not copying):** quality is not fundamental — it is the
sequencer's *estimate of how hard a base was to call*, and difficulty is driven by the **local DNA
sequence** (homopolymers, GC/GGC motifs, low-complexity). Everyone (fqzcomp/ACO/CRAM) predicts Q_i
from PAST QUALITY (q1,q2,delta,pos) — a lagging *proxy* for the real cause. We predict Q_i from the
**local sequence context itself**, which we uniquely have clean access to because our reads are aligned
to a self-assembled **pseudogenome** and the sequence stream decodes BEFORE quality (so the decoder has
it for free — lossless, zero stored bytes).

**Held-out measurement (50/50 train/test cross-entropy, so NOT overfitting):**
| quality context | test bpq | vs quality-history |
|---|---|---|
| pos,q1,q2 (fqzcomp's core idea) | 2.559 | — |
| + homopolymer run-length | 2.534 | −0.97% |
| + current base | 2.509 | −1.98% |
| **+ local 3-mer sequence context** | **2.287** | **−10.65%** |
| sequence-only (no quality history) | 2.584 | +0.97% |
3-mer is the single-model sweet spot; wider (5/7/9-mer) needs model mixing (future).

**REAL CODED RESULT (default, lossless, 7/7 ctests):**
| coder | DS7-100k quality | bpq | Δ vs static |
|---|---|---|---|
| static-rANS (old) | 4,807,480 | 2.599 | — |
| adaptive Q-ary (quality-history only) | 4,711,835 | 2.547 | −1.99% |
| **adaptive Q-ary + SEQUENCE context** | **4,460,085** | **2.411** | **−7.23%** |

**vs SOTA (same data): fqzcomp-q2 context floor = 2.419 bpq (in-sample, optimistic). We are at 2.411
REAL coded bits — at/below fqzcomp's best context**, using our own sequence-conditioning lever, fully
lossless. (fqzcomp binary not installed → compared against its context's measured entropy; our 2.287
held-out vs their 2.419 in-sample is an even stronger, stricter statement.) Total archive 5,695,934 →
5,348,539 (−6.1%).

**Why this is novel, not a copy:** classic fqzcomp/ACO do NOT condition quality on the base sequence;
CRAM/MPEG-G positional-quality needs an EXTERNAL reference + aligner and uses column *mean*;
Quartz/GeneCodeq use sequence but LOSSILY. Ours = reference-FREE (self-assembled pseudogenome) +
local sequence MOTIF context + fully LOSSLESS + unified codec. It is the quality-side twin of the
pseudogenome idea (exploit sequence structure), and it is what let us pass fqzcomp's context.

**Gate behaviour (honest):** on the 5k set the sequence context is too sparse (745K values over ~265K
contexts → cold-start), so it LOSES (+17%) and the keep-smaller gate correctly keeps static — no
regression. It wins only where there is enough data (100k). Both lossless. Default ON; disable seq via
ARCS_QUAL_NOSEQ, disable CM via ARCS_QUAL_NOCM.

**Next lever toward the 2.09 floor:** mix multiple sequence-context widths (3/5/7-mer) via deterministic
interpolation so the wide contexts contribute only where populated; and use the pg CONSENSUS motif
(cleaner than the per-read bases). Data-adaptive context width would also let 5k benefit.

---

## ⚑ ENTROPY SCOREBOARD + WIN (2026-07-16)

**Scoreboard — every model's score on OUR DS7-100k quality** (14.8M values, 35 levels, full-range
Illumina; entropy rows = bits/qual a perfect adaptive coder with that context reaches):

| model / tool | bits/qual | bytes | vs our old static |
|---|---|---|---|
| bzip2 -9 | 3.052 | 5,645,939 | +17% |
| xz -9 | 2.899 | 5,362,436 | +12% |
| q1 only (order-1) | 2.798 | 5,176,870 | — |
| q1+pos | 2.565 | 4,744,601 | ≈ our tier |
| **OUR OLD static-rANS** | **2.599** | **4,807,480** | baseline |
| fqzcomp-q1 (q1,max(q2,q3),pos) | 2.474 | 4,576,644 | −4.8% |
| fqzcomp-q2 (+delta) | 2.419 | 4,475,265 | −6.9% |
| 3rd-order + pos | 2.253 | 4,168,812 | −13.3% |
| 3rd-order + fine pos | 2.090 | 3,865,715 | −19.6% |

Our old static coder sat at the BASIC `q1+pos` tier; clear headroom below. We beat general tools
already but not the specialized context tiers.

**WIN — new adaptive Q-ary coder (`src/qual_cm.{h,cpp}`, quality mode 0x07, DEFAULT ON, lossless,
7/7 ctests):**

| coder | DS7-100k quality | Δ vs static | 5k | Δ |
|---|---|---|---|---|
| static-rANS (old) | 4,807,480 | — | 253,939 | — |
| **adaptive Q-ary (new, default)** | **4,711,835** | **−1.99%** | **251,565** | **−0.93%** |

Design = classic Subbotin range coder + **whole-symbol** per-context ADAPTIVE frequency table (no
transmitted model → rich context is free). Context is our own generative-model formulation:
{q1, q2, **decaying running-max ceiling**, **EWMA volatility**, pos8, is_dev} + **context inheritance**
(seed a fresh rich context from the warm low-order parent {q1,q2,pos}). Tunables stored in header:
FREQ_INC (adaptation speed, default 4; env ARCS_QCM_INC), INHERIT (default on; env ARCS_QCM_INHERIT).
Keep-smaller gate vs static means it can never enlarge the archive.

**What FAILED first (recorded honestly):** binary-decomposition + 4-model logistic mixer — tied static
at best, exploded to +490% at wrong learning rate; abandoned as finicky. The Q-ary whole-symbol
adaptive coder is the design that actually won.

**Gap to floor / next lever:** we reach 2.548 bpq vs the context's 2.448 floor — ~0.1 bpq is inherent
one-pass adaptation lag (inheritance bought only +0.03%, so it is NOT cold-start). Pure one-pass
adaptive cannot reach the −5.8% floor. Next: deterministic **2-model interpolation** (integer blend of
rich + low-order frequencies, PPM-style — stable, unlike logistic mixing) and richer context (q3,
finer pos). CAVEAT: `test5k.fastq` is CRLF — strip `\r` before paste|sort|cmp or it false-fails.

Everything below is the pre-build SOTA survey + design.
---

## 1. State of the art

### fqzcomp / fqzcomp5 (James Bonfield) — the reference SOTA
Power comes entirely from the **context** used to predict quality Q_i, coded with an **ADAPTIVE**
arithmetic coder — the frequency model updates online and is **never transmitted**, so adding context
dimensions is essentially free. Context for Q_i:
- Q_{i-1} (previous quality, exact)
- **max(Q_{i-2}, Q_{i-3})** (second-order history)
- boolean **Q_{i-2} == Q_{i-3}**
- a **running delta / variability** measure (a few bits; how bumpy the recent run is)
- **position** in read (fine-grained; param `qloc`)
- (aligned/CRAM mode) a match/selector context
Tunable in fqzcomp5 via qbits/qshift/qloc. Levels -q1..-q3 add these progressively.
Source: https://github.com/jkbonfield/fqzcomp5 , https://gensoft.pasteur.fr/docs/fqzcomp/4.6/

### ACO (BMC Bioinformatics 2022) — ~5% better than fqzcomp
- Adaptive coding ORDER: traverse quality symbols along the most correlated trajectory so similar
  symbols cluster.
- Compound context including the **read's global average quality**.
- 2D flowcell spatial correlation (cross-talk between adjacent clusters; needs tile/x/y from names).
Source: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04712-z

### Practical notes
- Modern NovaSeq quality is often binned to ~4-8 levels → tiny alphabet → the context model dominates.
- CRAM 3.1 (2022) uses fqzcomp-style adaptive quality with a match/reference selector context.

## 2. Where ARCS quality is today (the gap)
Current model: `ctx = q_prev_exact[0..42] × is_dev[0..1] × pos_bin[0..9]` = 860 contexts, static rANS,
and **the model is serialized into the archive**. Two gaps vs SOTA:
1. Missing context dims: no max(q2,q3), no equality bit, no running-delta, coarse position.
2. Architecture: static + serialized model → every extra context dim costs model-transmission bytes.
   (This is exactly why idea B's richer context lost: model overhead ate the gain.)
One thing we have that fqzcomp lacks reference-free: the **mismatch-vs-pg flag** (additive).

## 3. Design — adaptive CM quality coder (mode 0x07)
Switch to a fully ADAPTIVE model (no transmitted model), then enrich context to fqzcomp-parity +
our extra dims. Reuse dna_coder's carryless binary arithmetic coder + squash/stretch bit models.

- Code each quality value by binary decomposition (ceil(log2(Qmax+1)) bits); each bit predicted by an
  adaptive 12-bit probability selected by (context, bit-index, partial-bits). Encoder and decoder
  update the same bit-model after each bit → symmetric, lossless, **zero model bytes**.
- Context (hashed): {q1 exact, max(q2,q3) quantized, (q2==q3) bit, running-delta quantized,
  position bin (finer, ~16), is_dev mismatch flag, [phase 2: read-average-quality bin]}.
- New quality mode byte 0x07; encoder keeps whichever is smaller (current vs adaptive) → no regression.
- Losslessness: decoder rebuilds identical context (q-history from decoded symbols; is_dev from
  dev_sets; read-average needs 1 pass — computable since quality decoded per read).

### Phases (each gated on measured lossless gain)
- Q1: adaptive coder + {q1, max(q2,q3), q2==q3, delta, pos16, is_dev}. Expect the big win (fqzcomp parity).
- Q2: add read global-average-quality context (ACO idea).
- Q3 (optional, needs names): flowcell tile/x/y 2D context (ACO idea).

### Risks / rules
- Adaptive bit-model must update identically on both sides (integer only) — reuse dna_coder's proven
  update. Gate every phase on measured bytes; keep-smaller guarantees no regression. Verify
  losslessness via paste|sort|cmp on DS7-100k each phase.

## 4. Sources
- fqzcomp5 — https://github.com/jkbonfield/fqzcomp5
- fqzcomp 4.6 docs — https://gensoft.pasteur.fr/docs/fqzcomp/4.6/
- ACO — https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04712-z
- Human short-read benchmark 2025 — https://www.nature.com/articles/s41598-025-00491-8
- fastqz/fqzcomp context description — http://www.mattmahoney.net/dc/fastqz/fastqz.pdf
