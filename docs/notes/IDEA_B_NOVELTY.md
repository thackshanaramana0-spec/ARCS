# Idea B — does it actually make our novelty stronger? (honest ceiling)

Date: 2026-07-16. Direct answer to: "will idea B make it stronger / make novelty?"

## ⚑ BUILT AND MEASURED (2026-07-16) — NULL RESULT
Idea B was implemented end-to-end (lossless, 7/7 tests pass) and measured on DS7-100k. **It does
not help in practice.** Order sweep of the pg-surprise quality context (total quality bytes vs the
plain mismatch-flag baseline of 4,807,480 B):

| surprise FCM order | quality bytes | Δ vs plain |
|---|---|---|
| 2 | 4,807,158 | **−0.01%** (322 B — noise) |
| 3 | 4,822,088 | +0.30% |
| 4 | 4,825,037 | +0.37% |
| 6 | 4,839,642 | +0.67% |
| 8 | 4,847,498 | +0.83% |
| 11 | 4,847,865 | +0.84% |
| 14 | 4,845,224 | +0.79% |

Best case is −322 bytes on a 4.8 MB stream; everything above order 2 is worse (added model contexts
with no signal). **Interpretation:** the mismatch-vs-pg flag (already a quality context since the
chain-pg work) already captures the useful cross-stream signal; the pseudogenome's *predictability*
is essentially uncorrelated with Illumina quality. So idea B, as built, does NOT strengthen us
empirically. It is correct, lossless, and gated to never regress — but **OFF by default** (opt-in via
`ARCS_QUAL_SURP`) so it doesn't cost 2× quality-encode time for ~0 gain. This is a valid, honest
measured negative — recorded, not hidden.

Everything below is the pre-build analysis, kept for context.
---

## Short answer
**Yes, it strengthens our position — but it does NOT make us foundationally novel.** The general idea
"use sequence predictability to help compress quality" is heavily prior-arted. Our defensible slice is
narrow and specific. Build it because it is (a) our most-differentiated angle and (b) aimed at ~80% of
the archive — not because it is a new invention.

## What idea B is
While the sequence coder walks each read against the self-assembled pseudogenome it knows, per base:
(1) a **mismatch-vs-pg flag**, and (2) the CM model's **predicted probability / surprise** for the
called base. Feed both as CONTEXT into the LOSSLESS quality coder, in one unified integer codec.

## Prior art that limits the claim (found 2026-07-16)
- **Predict quality from a model, code the residual:** patent US 12080384 — "a prediction quality
  score value is computed at the current position, and the difference between predicted and actual is
  passed to an entropy coder… count-based or neural-network prediction-based arithmetic coding." This
  is the *general* prediction-based quality-coding idea. **Prior art.**
- **Match/mismatch-vs-reference as quality context:** patent US 12125562 (MPEG-G), for aligned data.
  **Prior art.**
- **Base predictability → drop/smooth quality:** Quartz, GeneCodeq, sequence-based smoothing
  (BMC 2019) — but LOSSY. **Prior art (lossy).**
- **Cluster + Markov-mixture quality models:** PMC5649045 — context = cluster id + position +
  previous quantized value. **Prior art.**
- **CM on quality generally:** fqzcomp, ACO, ZPAQ. **Prior art.**

## The narrow slice that is (so far) NOT covered by what we found
All four together:
1. **Reference-FREE** — the "reference" is a self-assembled short-read pseudogenome, no external genome
   and no aligner (the MPEG-G/CRAM context work assumes aligned data against a real reference).
2. **Real-valued model SURPRISE as the context**, not just a binary match/mismatch flag — the quality
   coder is conditioned on *how confident the DNA CM model was*, quantized into buckets.
3. **Fully LOSSLESS** (the strong prior-art predictability work is lossy smoothing).
4. **ONE unified integer codec** — the same engine's internal state is exported directly as quality
   context; encode and decode recompute it bit-identically (dna_coder is deterministic all-integer).

I did not find a tool combining all four. That is the honest ceiling: a **mechanism + combination
novelty**, in the same category as PgRC's own contribution — not a new primitive.

## So how much stronger does it make us?
- **Without idea B:** our assembler ≈ NanoSpring's idea; our CM coder ≈ GeCo3's; our pseudogenome ≈
  PgRC/Leon lineage. The novelty is "a new *combination* for short reads" — real but thin, and easy
  for a reviewer to pattern-match to existing tools.
- **With idea B (if it yields measurable lossless gains):** we add one mechanism (surprise→quality)
  that none of the cited tools implement in this reference-free lossless unified form. That moves us
  from "a re-combination of known parts" to "a re-combination **plus** a specific coupling mechanism
  we can point to as ours." Still a combination/mechanism claim, but a materially stronger one.

## Honest caveats before building
- It only counts as novelty **if it actually compresses** — a mechanism that doesn't beat the
  independent quality coder is not a contribution, it's a footnote. Gate on measured lossless gain.
- Frame it precisely as: *reference-free, surprise-conditioned (not flag-conditioned), lossless,
  unified.* Drop any one of those four qualifiers and it collapses into existing prior art.
- Do NOT claim "first to use sequence info for quality" — patents US 12080384 / 12125562 and Quartz
  own that framing. Claim only the specific mechanism + combination.

## Verdict
Build idea B next. It is our strongest genuinely-differentiated angle and hits the 80% of the archive
that sequence work can't. But log it as a **mechanism/combination novelty with a measured-gain gate**,
not as a foundational invention. Design + milestones: `PLAN_ideaB_crossstream.md`.

## Sources
- Prediction-residual quality patent US 12080384 — https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12080384
- MPEG-G mismatch-context patent US 12125562 — https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12125562
- Cluster + Markov quality models — https://pmc.ncbi.nlm.nih.gov/articles/PMC5649045/
- Sequence-based quality smoothing (lossy) — https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2883-5
- GeneCodeq — https://academic.oup.com/bioinformatics/article/32/20/3124/2196578
- ACO quality CM — https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04712-z
