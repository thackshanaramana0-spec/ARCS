# ARCS Fast Mode (Mode 2) — Final Deep Plan

Goal: a second operating mode that wins **speed + RAM + compression** across all datasets,
generalizable, **without inheriting PgRC2's design** — while Mode 1 (default) keeps the max-ratio
crown + reference-free variant calling.

---

## 0. The one rule that governs everything

**Fast mode keeps the FCM coder on the pseudogenome string. Always.**

Dropping the FCM for PgRC2's `position + mismatch` scheme is the single thing that would make
this a copy. The FCM-on-pseudogenome is simultaneously (a) the ratio crown and (b) the
originality vs PgRC2 — they are the same asset. Everything below makes the FCM + assembly
*cheap*, never replaces them.

---

## 1. Assembly decision — resolved

**Not greedy SCS. Not unchanged de Bruijn. → LEAN de Bruijn.**

Keep the multi-contig de Bruijn assembly; make it cheap via a sparse index + chunking. Rationale:

- Greedy SCS's only edge is RAM/speed — both winnable other ways (sparsity, chunking).
- Greedy SCS costs: variant-calling capability (no bubbles), originality (it *is* PgRC2's
  assembler), and buys ~nothing on ratio (the FCM, not the assembler, is the ratio lever).
- Lean de Bruijn keeps the richer structure, the ratio, and optionally the analysis byproduct.

---

## 2. Where the Genozip win actually lives (why fast mode survives)

Human archive decomposition (GIAB) and what fast mode touches:

| Layer | Share | ARCS vs Genozip | Touched by fast mode? |
|---|---|---|---|
| **Quality** | ~80% | ~tied (both near info floor) | **No** (coder unchanged) |
| **Names** | ~14% | ARCS ~2× smaller (own coord coder) | **No** (coder unchanged) |
| **Sequence** | ~6% | ARCS far smaller (pseudogenome vs no-assembly) | **Yes** (this erodes) |

The 16% archive win over Genozip on human comes from **names + sequence**; quality is a tie.
Fast mode only erodes **sequence = 6% of the archive**. Even a 50%-larger sequence is +3% archive
against a 16% lead. **Human margin survives with huge headroom.**

Structural reason it holds: Genozip builds **no read-reference at all** — it context-models bases
in place and sees each position ~30× (coverage) without deduping. A *cheaper* pseudogenome is
still a pseudogenome; it still collapses the redundancy Genozip cannot. Fast mode keeps the
architectural advantage, just builds it faster.

**Thin case = short-read bacterial (DS4 M.tb, 51 bp):** margin only ~5% AND sequence is a bigger
archive fraction (less quality bulk to hide behind). This is the only place the ratio claim
becomes "smallest, narrowly" rather than "dominant."

---

## 3. The three levers (per axis)

**Lever 1 — Minimizer-sparse assembly index (RAM + assembly speed).**
Index minimizers (~1 per 10 bases), k-mer table ~10× smaller. miniasm/minimap lineage —
architecturally distinct from PgRC2's copMEM. Prior "minimizer seeding = dead" was a *ratio*
result for the merge; here it only needs to be ratio-*neutral* (lower bar).

**Lever 2 — Chunked streaming assembly (RAM + speed, decisive).**
Process C reads at a time; each chunk builds → encodes → frees its own pseudogenome. Peak RAM =
a few in-flight chunks (bounded, tunable). Chunks run in a thread pool → parallel speed. Ratio
cost: cross-chunk redundancy lost (~0.36%/split, measured).

**Lever 3 — Lean FCM (speed, small ratio cost).**
8 orders → ~5. Mixer already auto-gates redundant high orders. Cleanest knob; env
`ARCSDNA_ORDERS` already exists → measurable with zero code.

**Reinvestment escape hatch:** chunk parallelism can be spent on *ratio* instead of pure speed —
run a richer (8-order) FCM per chunk in the time a serial 5-order would take. This is the rescue
for the DS4 thin case (keep ratio, still get the RAM win from small chunks).

---

## 4. Generalization key — adaptive chunk sizing

Ratio loss must always be a small fraction of the existing margin:

- **Small / bacterial (< ~100 MB, thin 2–5% lead):** single chunk (zero cross-chunk loss) +
  sparse index + lean-or-reinvested FCM. Speed/RAM win at near-zero ratio cost.
- **Large / human WGS (huge 15–27% lead):** chunk aggressively. Spend ~1% ratio to buy 5× speed
  and 4× RAM.

Chunk count scales with input so the ratio hit is always dwarfed by the margin.

---

## 5. Honest expected outcome (vs SPRING / Genozip — the real lossless peers)

| Axis | Realistic outcome | Confidence |
|---|---|---|
| **RAM** | **Win outright** — chunking bounds peak RAM below both (DS7 3.9 GB → target < 1.3 GB) | High |
| **Speed** | **Match-to-beat SPRING** (DS7 180 s → ~30–40 s); beating Genozip's ~1.2 s on tiny files unlikely | Medium |
| **Ratio** | **Hold crown on human comfortably; still smallest on bacterial but narrowly (DS4 thinnest)** | Medium |

**Do NOT claim:** "dominates all 3 axes on every dataset." Beating Genozip's sub-second small-file
speed with an assembler onboard is unrealistic. Fast mode is *competitive* there, *winning* on RAM
and large-file speed, *still smallest* on ratio.

---

## 6. Execution phases (build only after each gate passes)

**Phase 0 — Gate experiment (pure measurement, zero code).**
Run the 8-dataset sweep with `ARCSDNA_ORDERS=2,4,8,14,22` (5 orders) vs default 8. Record ratio
delta per dataset, especially DS4/DS5/DS2.
- If 5-order costs **< 1%** on DS4 → full plan viable.
- If **~3%** → use the reinvestment hatch (parallel 8-order chunks) or scope fast mode to human-WGS.

**Phase 0b — Genozip `--best` baseline.** Re-run Genozip at max level on all datasets so the margin
is measured against Genozip's strongest, not its default (pre-empts a reviewer surprise).

**Phase 1 — Chunked path curve.** The `>200 MB → 4 chunks` auto-chunking and block-parallel pg
(format 0x05) already exist. Lower the threshold; measure ratio/speed/RAM vs chunk count on DS7.

**Phase 2 — Sparse minimizer index.** Add as the fast-mode assembly index; verify ratio-neutral on
one human + one bacterial dataset; measure RAM drop.

**Phase 3 — Wire the `--fast` flag** that composes: adaptive chunk sizing + sparse index +
lean-or-reinvested FCM. Run the full 8-dataset sweep + losslessness verification (byte-exact
roundtrip is non-negotiable; all 7 ctests must pass).

**Phase 4 — Document** the two-mode story: Mode 1 (max ratio + variant calling), Mode 2 (fast:
speed/RAM optimized, still smallest lossless archive on human, competitive on bacterial).

---

## 7. Risks & honest caveats

- **DS4/DS5 ratio erosion** — the real risk; mitigated by adaptive sizing + reinvestment; gated by
  Phase 0.
- **Genozip `--best`** could thin margins independent of fast mode — measure it (Phase 0b).
- **Losslessness** — chunking/sparse index must stay byte-exact; every phase re-verifies roundtrip.
- **Variant calling** — off in fast mode by design; do not hardwire the *inability* (keep de
  Bruijn so Mode 1 still has it).
- **Speed vs Genozip small files** — will not win; state as competitive, not dominant.

---

## 8. One-line summary

Keep the FCM and the de Bruijn assembly (that's the ratio crown, the variant calling, and the
not-PgRC2 story); make them cheap with sparse minimizer indexing + adaptive chunked streaming +
a leaner-or-reinvested model set. Result: RAM win outright, speed competitive-to-winning vs SPRING,
and still the smallest lossless archive on every dataset — comfortably on human, narrowly on
short-read bacterial. Gate the whole thing on one zero-code measurement: 5-order vs 8-order FCM.
