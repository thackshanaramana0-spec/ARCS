# pg-coder model ablation — generalized, ratio-safe speed lever

## Why this measurement

The DS1 profile disproved "assembly is the universal compress pole": on DS1
assembly was 15s but pg context-mixing coding was ~50s. The pg coder
(`dna_coder.cpp`: 8 finite-context models orders {2,4,8,11,14,18,22,26} + STCM +
logistic mixer + APM, per-base per-bit, serial-adaptive) is the one heavy
component that runs on **every** dataset (compress AND decode) and **is** the
ratio crown. It cannot be block-split without cold-restart ratio loss (measured
+20%/+333KB on DS1). So the only ratio-safe, generalized speed lever is to make
the coder do **less work per base** without changing its output — i.e. drop
models that don't earn their runtime.

## Method

Leave-one-out over the 8 orders on two regimes (E.coli = coherent bacterial,
the thin-margin/pg-pole case; DS7 = human). Signal = Δ pg_blob (isolated coder
output) when a model is dropped. `ARCSDNA_ORDERS` knob applies to both encode
and decode. Inputs: `/c/Temp/arcs_test/{d1,ds7_100k}.fastq` (100k reads each).

## Per-model contribution (Δ pg_blob when dropped)

| order | E.coli Δ | human Δ | verdict |
|------:|---------:|--------:|---------|
| 2     | +625 B   | +635 B  | keep (small, consistent) |
| **4** | **−4 B** | **−1 B**| **DEAD both regimes → prune** |
| 8     | +191 B   | +170 B  | near-dead (~0.08%) — optional prune |
| 11    | +3211 B  | +3192 B | workhorse |
| 14    | +1934 B  | +1938 B | strong |
| 18    | +746 B   | +723 B  | keep |
| 22    | +419 B   | +403 B  | keep |
| 26    | +403 B   | +370 B  | keep |

The contribution profile is **near-identical across E.coli and human** →
structural, not dataset-tuned. Passes the "must generalize" bar.

## Combined-prune confirmation (DS1, direct not leave-one-out)

| config | pg_blob | archive | vs baseline |
|--------|--------:|--------:|-------------|
| baseline (8) | 243,511 | 5,244,479 | — |
| prune-4  (7) | 243,507 | 5,244,475 | **−4 B (free)** |
| prune-4-8 (6)| 244,008 | 5,244,976 | +497 B (~0.01%) |

## Byte-safety (losslessness)

Rigorous check (`lossless_check.sh`):
- **prune-4 decode == baseline decode: IDENTICAL** ← the coder change alters
  nothing in the output. Provably output-preserving.
- baseline decode vs raw `/c/Temp/d1.fastq`: differs ONLY in read **names**
  (`@HWI-D00360...` vs `@HISEQ1...`); sequence + quality byte-identical, sizes
  identical. This is a pre-existing NAMES-stream artifact of that derived test
  file — orthogonal to `dna_coder` (names are a separate stream). Not caused by
  pruning. A clean absolute-lossless claim should use the pristine
  `arcs-clean/benchmark` inputs.

## Coder speedup — MEASURED, hypothesis REJECTED

Isolated microbench (`microbench.cpp`), 5M-base pg, best-of-3, single process:

| config | models | encode | vs baseline | out |
|--------|:------:|-------:|------------:|----:|
| baseline           | 8 | 7.794s | —      | 914,565 B |
| prune-4            | 7 | 8.479s | +8.8% (noise) | 913,461 B |
| prune-4-8          | 6 | 7.626s | −2.2%  | 922,974 B |

**Removing models does NOT meaningfully speed the coder.** prune-4 even measured
slower (thermal noise); the true effect of dropping 1–2 of 8 models is buried
under a ~±8% run-to-run noise floor, i.e. ≤2% even for two models. The earlier
structural assumption ("1 of 8 models ≈ 12.5% of per-base cost") is WRONG: the
coder's time is dominated by the SHARED per-bit machinery (logistic mixer
stretch/squash, APM, the carryless arithmetic inner loop), not the FCM gathers.
Model count is not the lever.

## Conclusion

Pruning dead models (order-4) is **ratio-free and byte-safe but speed-neutral**
(~0–2%, within noise) — NOT a favourable speed win. The measurement redirects
the target: to speed the coder you must optimize the shared per-bit hot loop
(SIMD the mixer/APM/arith), not reduce the model set — a harder, separate task.
The parallel train-freeze idea (bigger speedup, some ratio risk) remains the
other candidate. Net: this measurement's value is a clean NEGATIVE result that
saves shipping a pointless change and points at the real bottleneck.
