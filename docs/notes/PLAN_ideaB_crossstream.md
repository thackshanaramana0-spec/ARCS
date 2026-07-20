# PLAN — idea B: cross-stream quality conditioning

Status: BUILT + MEASURED 2026-07-16 → NULL RESULT (best −0.01%, usually worse). Shipped lossless,
7/7 tests, OFF by default (opt-in `ARCS_QUAL_SURP`, order via `ARCS_QUAL_SURP_ORDER`, stored in
meta[8]). See `IDEA_B_NOVELTY.md` for the measured order sweep and interpretation. The mismatch-flag
context already captures the useful cross-stream signal; pg predictability is ~uncorrelated with
Illumina quality. Sections below are the original design, kept for reference.

## Goal
Cut the quality stream (~80% of the archive) by conditioning the quality coder on what the
sequence side already knows per base, in one lossless unified integer codec.

## Signals the sequence side already produces per base (free, at encode time)
1. **mismatch flag** — was this base a mismatch vs the pseudogenome? (from `pg_mm_pos_flat`)
2. **model surprise** — ARCS-DNA's predicted probability for the emitted base at that position
   (quantize −log2(p) into a few buckets). Requires exposing the per-base prob from `dna_coder`.
3. **position in read** (already used by the current quality model).
4. **previous quality value** (already used).

## Design
- Extend the quality context model in `encoder.cpp` (`encode_quality_rans`) to mix in
  (mismatch_flag, surprise_bucket) alongside (position, prev_qual).
- Decoder reconstructs the SAME context: at quality-decode time the sequence is already
  decoded, so mismatch_flag and surprise are recomputable → context is symmetric → lossless.
- Keep it a rANS/CM context table; no change to the sequence path.

## Milestones
1. **A-first baseline** (cheap): run current independent quality coder, record DS7-100k quality bytes.
2. Expose per-base predicted prob from `dna_coder` (encode side) + recompute on decode side.
3. Add mismatch_flag context only → measure quality delta.
4. Add surprise_bucket context → measure quality delta.
5. Gate: keep only if quality shrinks with zero losslessness regression (project rule:
   "quality or compression no acceptable" to regress).

## Risks / notes
- Surprise recomputation on decode must be bit-identical to encode (same model state, same
  integer arithmetic) or losslessness breaks. dna_coder is already all-integer/deterministic — good.
- Expect the mismatch-flag context to help most (mismatch bases carry distinct quality profiles);
  surprise is the incremental, more-novel part.
- Novelty framing per PRIOR_ART.md: reference-free + surprise + lossless + unified. Don't overclaim.
