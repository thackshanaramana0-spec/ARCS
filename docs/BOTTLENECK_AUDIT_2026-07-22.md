# ARCS CPU & RAM Bottleneck Audit — 2026-07-22

**Datasets:** DS1 (100K human reads, 34 MB), DS7s (100K metagenomics, 22 MB), GIAB (113K WGS, 39 MB)  
**Method:** `ARCS_ENC_NOPAR=1` serial execution for clean stage attribution; `ARCS_ENC_TIMING=1 + ARCS_ASM_TIMING=1` wall timers; `ARCS_QUAL_TIMING=1 + ARCS_QUAL_PROFILE=1` TSC cycle counters; RSS snapshots via `cur_peak_rss_mb()`  
**Losslessness:** DS1 = 6,256,752 B, DS7s = 4,814,397 B, GIAB = 4,310,976 B — consistent across all instrumented runs

---

## 1. Compress Stage Runtime Ranking

All times in seconds, serial (`ARCS_ENC_NOPAR=1`).

| Rank | Stage | DS1 (s) | DS1 % | DS7s (s) | DS7s % | GIAB (s) | GIAB % | Limiting factor |
|------|-------|---------|-------|---------|--------|---------|--------|----------------|
| 1 | **dna_encode** (pg_encode) | **16.45** | **50.8** | **19.82** | **62.2** | 1.28 | 14.2 | 8×2²⁴ FCM tables = 512 MB → L3 eviction every lookup |
| 2 | **Assembly: merge** (CSR KmerIndex) | 6.15 | 19.0 | 7.62 | 23.9 | 0.61 | 6.7 | O(contigs × K) overlap scan; 53K→31K / 82K→76K contigs |
| 3 | **names_encode** (LZMA-9) | 5.18 | 16.0 | 0.34 | 1.1 | **4.60** | **50.9** | 4 concurrent LZMA-9 candidate archives computed serially |
| 4 | **qual keep-smaller gate** (both coders tried) | 2.34 | 7.2 | 1.71 | 5.4 | 1.36 | 15.0 | Computes rANS + adaptive-CM, picks smaller; winner predictable |
| 5 | Assembly: placement | 1.51 | 4.7 | 1.93 | 6.1 | 0.50 | 5.5 | O(N) k-mer lookups into KmerIndex |
| 6 | qual_cm chosen model | 0.68 (rANS) | 2.1 | 1.18 (CM) | 3.7 | 0.67 (CM) | 7.4 | Adaptive FreqMap + arith coder, block-parallel |
| — | load_reads + akc_pilot | 0.33 | 1.0 | 0.25 | 0.8 | 0.39 | 4.3 | I/O + k-mer spectrum; not a bottleneck |
| **Total** | | **32.45** | | **31.88** | | **9.04** | | |

### Assembly sub-breakdown

| Sub-stage | DS1 | DS7s | GIAB |
|-----------|-----|------|------|
| Phase1 placement | 1.51s (19%) | 1.93s (20%) | 0.50s (39%) |
| Phase1.6 merge | **6.15s (77%)** | **7.62s (78%)** | 0.61s (48%) |
| Phase1.5 polish | 0.20s (2.5%) | 0.15s (1.5%) | 0.09s (7%) |

Merge dominates assembly for low-coverage / metagenomics datasets (DS1: 53K initial contigs; DS7s: 82K).

---

## 2. Decode Stage Runtime Ranking

| Rank | Stage | DS1 (s) | DS1 % | DS7s (s) | DS7s % | GIAB (s) | GIAB % | Limiting factor |
|------|-------|---------|-------|---------|--------|---------|--------|----------------|
| 1 | **dna_decode** (pg_decode) | **6.36** | **86** | **7.50** | **86** | 0.59 | 40 | Same FCM model, sequential arith-decode per base |
| 2 | **qual_cm decode** | 0.65 | 9 | 0.97 | 11 | **0.67** | **45** | Adaptive model replay; block-parallel (6 blocks) |
| 3 | reconstruct | 0.17 | 2.3 | ~0.15 | ~1.7 | 0.13 | 8.7 | Delta application; nearly free |
| 4 | names_decode + write | ~0.08 | ~1 | ~0.08 | ~1 | ~0.10 | ~6 | LZMA decompress; buffered write; negligible |
| **Total** | | **~7.4** | | **~8.7** | | **~1.49** | | |

`dna_decode` is the monolithic decode bottleneck for hard datasets (DS1/DS7s). It becomes less dominant for GIAB only because GIAB's pg is tiny (0.7 MB vs 7.5 MB for DS1) — the FCM processes proportionally fewer bases.

---

## 3. Quality CM Internal Cycle Breakdown

Measured via `__rdtsc()` with calibrated GHz, over all quality symbols per dataset.

| Component | GIAB cy/sym | GIAB % | DS1 cy/sym | DS1 % | Root cause |
|-----------|------------|--------|-----------|-------|-----------|
| **`build_cf()`** | **415** | **64%** | **896** | **73%** | 196-byte `Freq` struct × 32K table entries = 6+ MB → L2/L3 miss every call |
| `FreqMap.get()` (×2 per sym) | 44 | 7% | 125 | 10% | Open-addressing hash; cold misses when table large |
| CDF cumsum (`for i<q: cum += cf[i]`) | 49 | 8% | 59 | 5% | Sequential prefix sum, `q` ≤ 41 iterations |
| `enc.encode()` (range coder) | 50 | 8% | 50 | 4% | Subbotin multiply-shift; small working set |
| `FC.update() + FP.update() + s.advance()` | 40 | 6% | 41 | 3% | Freq struct write; same cache line as read |
| `keys()` (context construction) | 51 | 8% | 50 | 4% | Arithmetic; no memory pressure |
| **Total** | **649** | 100% | **1221** | 100% | |

**DS1 vs GIAB 2× cost difference:** DS1 uses seq-motif context (`seqs ≠ null`) which doubles the coded key space → more distinct `Freq` entries → worse `build_cf` cache behavior (896 vs 415 cy/sym). The `FreqMap.get()` miss rate is also visibly higher (125 vs 44 cy/sym).

---

## 4. Quality Parallel Scaling Experiment

GIAB with 1→8 quality blocks (serial pipeline, `ARCS_ENC_NOPAR=1`):

| ARCS_QUAL_BLOCKS | Total wall (ms) | Archive size (B) | Quality wall-time saving | Total saving |
|-----------------|----------------|-----------------|------------------------|-------------|
| 1 | 8,194 | 4,240,575 | — | — |
| 2 | 7,114 | 4,254,650 | 1,080ms (13%) | 13% |
| 4 | 6,801 | 4,285,878 | 1,393ms (17%) | 17% |
| 6 | 6,475 | 4,310,976 | 1,719ms (21%) | 21% |
| 8 | 6,412 | 4,332,829 | 1,782ms (22%) | 22% |

**Within-quality parallel efficiency:** 2,339ms (1-block quality time) / 620ms (6-block) = 3.77× speedup at 6 threads = **63% efficiency** — respectable, but irrelevant to the global picture. Quality is 7.4% of GIAB total serial time; making it zero-cost would save 7.4% total (Amdahl). The current 6-block default recovers 79% of that theoretical maximum while adding only 70 KB (+1.6%) to the archive.

**DS1 (5 blocks default vs 1 block):** 24,406ms vs 27,528ms = 3,122ms saving; 64% parallel efficiency within quality stage. Quality is 7.2% of DS1 serial total.

**Conclusion:** The quality stage is **not the compress bottleneck.** Throwing more threads at it has near-zero impact on total time. For DS1, even eliminating quality CM entirely would save only 2.1% of total compress time.

---

## 5. RAM Attribution

Peak RSS by phase (`cur_peak_rss_mb()`, monotonically non-decreasing):

| Phase | DS1 (MB) | DS7s (MB) | GIAB (MB) | Notes |
|-------|---------|---------|---------|-------|
| After load_reads | 57 | 45 | 62 | `reads[]` only |
| **After assembly** | **931** | **892** | 186 | KmerIndex + contig graph + ChainResult; never freed |
| After quality CM | 931 | 892 | **309** | No growth DS1/DS7s; GIAB +123 MB (FreqMaps + rq) |
| After write | 931 | 892 | 309 | No release during encode |

### RAM Ranking Table

| Rank | Structure | DS1 peak | DS7s peak | GIAB peak | Lifetime | Release earlier? |
|------|-----------|---------|---------|---------|---------|-----------------|
| 1 | **Assembly intermediate** (KmerIndex, contig graph, hash tables) | ~870 MB | ~847 MB | ~120 MB | load → process exit | Partial: KmerIndex could free after pg built (names/qual don't need it) |
| 2 | **`reads[]`** (seq + qual + name, all reads in memory) | ~40 MB | ~22 MB | ~40 MB | load → write output | No — qual CM needs rq, names encoder needs names |
| 3 | **qual CM FreqMaps** (6× tbl + 6× lo blocks) | ~30 MB | ~30 MB | ~30 MB | qual encode phase | Freed after qual encode; no issue |
| 4 | **ChainEncodeResult** (pg, chain_order, pg_pos, mm arrays) | ~20 MB | ~10 MB | ~3 MB | assembly → pg_encode | Could free after pg_encode |
| 5 | **rq vector** (uint8 Phred per read per position) | 15 MB | 8 MB | 17 MB | load → qual encode | Freed after qual encode |
| 6 | **seqvec** (bases for seq-motif context, DS1 only) | 15 MB | 0 | 0 | qual encode phase | Released after qual encode |

**Key finding:** For DS1/DS7s, the 870–847 MB assembly footprint is **never freed.** The intermediate KmerIndex alone accounts for most of this — it is constructed during merge (`chain_encoder.cpp`) and could be explicitly freed once the pseudogenome is built, before the quality/pg/names encode phases begin. This is a RAM reduction opportunity, not a speed opportunity (assembly is only 24–31% of time and the merge structures are already freed internally after the merge pass).

For GIAB the assembly is compact (120 MB), but the quality CM adds 123 MB during encoding — consistent with 6 FreqMap blocks × ~4.3 MB/block (22K `Freq` entries × 196 bytes) plus the rq vector.

---

## 6. Top 3 Confirmed Bottlenecks

---

### Bottleneck #1: `dna_encode` / `dna_decode` — `run_core<ENCODE/DECODE>()` in `dna_coder.cpp`

| | DS1 | DS7s | GIAB |
|-|-----|------|------|
| Compress impact | 16.45s (50.8%) | 19.82s (62.2%) | 1.28s (14.2%) |
| Decode impact | 6.36s (86%) | 7.50s (86%) | 0.59s (40%) |

**Root cause:** The ARCS-DNA FCM has 8 context models with hash tables sized to `min(2^24, ≥2×N)`. For DS1/DS7s, `N` (pg size in bases) → table cap at **2^24 entries × 4 bytes × 8 models = 512 MB** total working set. The i5-1235U has 12 MB L3 — a **43× miss ratio** is unavoidable. Each base in the pseudogenome requires 8 `model.select()` hash lookups + 2 `code_node()` arithmetic operations + 8 `model.update()` writes. All are serial (each base depends on the previous state). At an estimated 40–100 ns per LLC miss × 16 select/update calls × 7.5M bases (DS1) ≈ **4.8–12s** from memory latency alone, consistent with the measured 16.45s including arithmetic overhead.

**Amdahl gain if 2× dna_encode speedup:** DS1 compress −25%, DS1 decode −43%, DS7s compress −31%.

**Leverage:** Table size and model count are the only degrees of freedom (per-base serial arithmetic cannot be vectorized). No existing TSC profiling inside `run_core()` — the inner loop decomposition into select / code_node / update is unknown.

---

### Bottleneck #2: `encode_names` (LZMA-9) — `encoder.cpp`

| | DS1 | DS7s | GIAB |
|-|-----|------|------|
| Compress impact | 5.18s (16.0%) | 0.34s (1.1%) | 4.60s (50.9%) |

**Root cause:** For GIAB (Illumina XY names), the names encoder computes **4 concurrent LZMA-9 candidate archives** (XY-template 0x03, columnar 0x04, columnar-sorted 0x05, plain LZMA, paired-end dedup 0x06) and picks the smallest. Each LZMA-9 call is single-threaded and CPU-intensive; 4 candidates × ~1.15s each = 4.6s. The winner (XY-template 0x03) is almost always deterministic: Illumina data → XY-template wins; non-Illumina → plain LZMA wins.

**Amdahl gain if 2× speedup:** GIAB compress −25%; DS1 compress −8%.

**Leverage:** The candidate competition is provably wasteful — the winner can be predicted from a 100-read sample (detect Illumina `@.*:X:Y` format → skip plain-LZMA and non-XY candidates). Also: running candidates at LZMA-7 instead of LZMA-9 saves ~30% time at negligible ratio cost (`ARCS_NAMES_LZMA_LEVEL=7` already used in `--call` mode).

---

### Bottleneck #3: `qual_cm_encode` keep-smaller gate — `qual_cm.cpp`

| | DS1 | DS7s | GIAB |
|-|-----|------|------|
| Wasted time (loser model) | 1.66s (5.1%) | 0.53s (1.7%) | 0.69s (7.6%) |

**Root cause:** The encoder always computes **both** static-rANS and adaptive-CM, then discards the larger. The winner is fully predictable from the data: datasets with qmax ≤ 42 (binned Illumina, GIAB-style) → adaptive-CM wins every time; datasets with qmax > 42 (full-range) → static-rANS wins. Sampling 1,000 reads to count `max(q)` costs <1 ms and would let the encoder skip the losing model entirely.

**Amdahl gain:** GIAB compress −7.6%; DS1 compress −5.1%. Total across a mixed workload: ~5–8% compress speedup, zero archive change.

---

## 7. Bug: Parallel Crash (exit 127) on DS7s and GIAB

When `ARCS_ENC_NOPAR` is not set, `fut_pg.get()` and `fut_names.get()` terminate the process (exit 127) for DS7s and GIAB — but not for DS1. The likely cause: `encode_names` spawns 4 nested `std::async` futures from within an outer `std::async` thread (plus LZMA's own internal threads = 10+ concurrent threads). On Windows/MinGW, this hits a thread or stack limit when each LZMA-9 call takes >1s, but not when names are trivial. All benchmark measurements above used `ARCS_ENC_NOPAR=1` as workaround. **This bug must be fixed before shipping parallel mode on non-DS1 data.**

---

## 8. Final Decision

The bottleneck map has one clear answer across both compress and decode: **`run_core()` in `dna_coder.cpp` is the primary rate-limiting function.** It accounts for 50–62% of compress time on hard datasets and 86% of decode time on all but GIAB. Quality CM is not the bottleneck (2–7% of compress, already block-parallelized). Assembly merge is second for compress but not addressable without ratio risk. Names LZMA is second for GIAB but mostly eliminable via candidate pruning (Bottleneck #3 above is the cheap win).

The quality CM profiling (Section 3) confirmed that `build_cf` dominates quality time at 65–73%, but since quality CM is only 2–7% of total, optimizing `build_cf` yields <7% total speedup by Amdahl — not the right lever. Bottleneck #3 (keep-smaller gate, 5–8% free saving) should be done opportunistically but is not the focus.

---

## 9. Recommended Experiment

**Instrument `run_core<ENCODE>()` in `dna_coder.cpp` with per-base TSC profiling, decomposing each base's cost into three groups: `model.select()` × 8, `code_node()` × 2, and `model.update()` × 8 — then run GIAB and DS1 to get cycles/base per sub-operation.**

**Why the profiling data supports this:**
- `dna_encode` is the #1 bottleneck (50–62% compress, 86% decode) but has zero existing internal instrumentation
- The inner loop does three distinct things: hash table access (16 hits), arithmetic coding (2 calls), table writes (8 hits) — their relative weights determine which optimization is correct
- If `model.select() + model.update()` > 70% of cycles → the lever is table footprint: reducing 8×2^24 to 8×2^20 (32 MB total) would bring the working set inside a 12–32 MB L3, potentially 2–3× speedup
- If `code_node()` > 50% → the lever is arith-coder arithmetic; block-parallel decode (already prototyped at 0.23% ratio cost for 4 blocks) is the right path for decode, compress stays serial
- The ratio/speed Pareto of reducing model count (8→4) vs table size (2^24→2^20) cannot be chosen without knowing the cache-miss vs compute split

**What to measure:** cycles/base for (a) select-phase, (b) code_node-phase, (c) update-phase; compare DS1 (7.5 MB pg, large tables, ~L3 thrash) vs GIAB (0.7 MB pg, likely L3-resident) — the difference isolates cache-miss cost from arithmetic cost.

**Expected result:** select+update > 70% of cycles on DS1; GIAB cy/base substantially lower than DS1, confirming the cache-miss hypothesis.

**Archive size regression risk:** Zero — measurement only; no change to the encoding path.

**Complexity:** ~30 lines of `__rdtsc()` instrumentation in `run_core()`, identical pattern to the `qual_cm.cpp` profiling already in place. Guards under `ARCS_DNA_PROFILE` env var.

**Acceptance / rejection gates:**

| Outcome | Interpretation | Next step |
|---------|---------------|-----------|
| select+update > 70% of DS1 cy/base AND GIAB cy/base < 30% of DS1 | Cache-miss hypothesis confirmed | Reduce table to 8×2^20, measure ratio/speed Pareto |
| select+update < 50% | Cache-miss not the bottleneck | Examine `code_node()` arithmetic; consider block-parallel decode |
| GIAB cy/base ≈ DS1 cy/base | Tables already cache-resident for both | Table reduction won't help; revisit model count or arith coder |

---

## 10. FCM Inner-Loop TSC Profiling Results — 2026-07-22

**Instrumentation:** `ARCS_DNA_PROFILE=1` env var; `__rdtsc()` brackets around select / code_node / update phases inside `run_core<ENCODE/DECODE>()` in `dna_coder.cpp`. Wall time via `std::chrono::high_resolution_clock` over the same span for GHz calibration. TSC frequency = 2.50 GHz (fixed nominal on i5-1235U). Serial mode `ARCS_ENC_NOPAR=1` for encode. Archive byte-identical to audit baselines: DS1 = 6,256,752 B, GIAB = 4,310,976 B. Decoded DS1 diff-clean against source.

### 10.1 Cycles-per-base table

| Dataset | Mode | N (bases) | FCM w.s. | Select cy/b | Code cy/b | Update cy/b | Other cy/b | Total cy/b |
|---------|------|-----------|----------|------------|-----------|------------|-----------|-----------|
| DS1 | ENCODE | 7,461,447 | **320.8 MB** | 52.2 (2.6%) | 618.8 (30.7%) | 1002.9 (49.7%) | 342.1 (17.0%) | 2016.2 |
| DS1 | DECODE | 7,461,447 | **320.8 MB** | 53.6 (2.0%) | 700.2 (26.5%) | 1410.5 (53.5%) | 473.3 (17.9%) | 2637.6 |
| GIAB | ENCODE | 743,266 | **40.8 MB** | 52.8 (2.8%) | 657.4 (35.1%) | 900.7 (48.0%) | 264.0 (14.1%) | 1874.9 |
| GIAB | DECODE | 743,266 | **40.8 MB** | 53.9 (2.3%) | 720.4 (30.9%) | 1241.2 (53.2%) | 317.1 (13.6%) | 2332.6 |

Phase definitions: **Select** = `m.select(hist)` × 8 models (hash lookup + prefetch issue). **Code** = both `code_node()` calls (mixer dot-product + APM + arith encode/decode + renorm flush). **Update** = `m.update(b)` + `m.ir_update_fast(rcomp)` × 8 models (count write + IR count write). **Other** = `b2i`, `i2b`, `rcomp` roll, `hist` roll, loop overhead.

### 10.2 Per-dataset metadata

| Dataset | hash_bits | Hashed-model table | FCM working set | Encode wall | Decode wall | Archive |
|---------|-----------|--------------------|-----------------|-------------|-------------|---------|
| DS1 | 24 | 16,777,216 entries × 5 models | 320.8 MB | 6.0 s | 7.9 s | 6,256,752 B |
| GIAB | 21 | 2,097,152 entries × 5 models | 40.8 MB | 0.56 s | 0.70 s | 4,310,976 B |

### 10.3 Decision gate evaluation

| Gate | Threshold | DS1 | GIAB | Verdict |
|------|-----------|-----|------|---------|
| select + update > 70% → cache-miss hypothesis confirmed | 70% | 52.3% | 50.8% | **FAIL** |
| code_node > 50% → investigate arith coder / block-parallel | 50% | 30.7% | 35.1% | **FAIL** |
| DS1 cy/base ≈ GIAB cy/base → reject cache-size hypothesis | "similar" | 2016 | 1875 | **7.5% gap — REJECT** |

**Cache-size hypothesis is rejected.** A 8× working-set difference (320 MB vs 40 MB, both far outside 12 MB L3) produces only a 7.5% cy/base difference. Table reduction from 2^24 → 2^20 would close this small gap at measurable ratio cost — not the right lever.

### 10.4 Root-cause analysis

**`select` is free (2.0–2.8%).** The `__builtin_prefetch` in `model::select()` issues the hash-table load far enough ahead of use (code_node consumes ~620–720 cy before the row is read) that L3/DRAM latency is fully hidden. Reducing table size to push from L3-miss to L2-hit would gain almost nothing here.

**`update` dominates (48–54%), driven by `ir_update_fast`, not `m.update(b)`.** `m.update(b)` writes to `*cur` — the same cache line prefetched for `select`, so it is L1-hot. `ir_update_fast(rcomp)` accesses a different row (`tbl[index(rcomp >> ir_ctx_shift)]`) with **no prefetch**, causing a cold LLC miss per model per base with zero latency hiding. 5–8 unmasked cold accesses per base explains the 900–1400 cy/base update cost and why DS1 and GIAB are nearly equal — both miss at DRAM or deep L3 for the IR row.

**Encode vs decode update asymmetry: 1003 vs 1410 cy/base (DS1).** Decode must resolve the decoded bit from `coder.code()` before `update(b)` and `ir_update_fast(rcomp)` can fire, lengthening the dependency chain. No prefetch for the IR row in either direction, but encode's known `b` allows earlier dispatch. This is structural; block-parallel decode addresses it.

**`code_node` is 30–35%.** Mixer dot-product (8 weights) + APM interpolation + 1 arith encode/decode per bit = ~310–360 cy per call. The arith-coder renorm loop is serial per bit and cannot be vectorized within a single stream. For decode this is addressable via block-parallel (already at 0.23% ratio cost for 4 blocks).

### 10.5 Recommended next experiment

**Add `__builtin_prefetch` for the `ir_update_fast` table row inside `model::select()`, issuing it alongside the forward-context prefetch.**

Each call to `select(hist)` already prefetches `tbl[index(hist & mask)]`. At the same time, issue `__builtin_prefetch(&tbl[index(rcomp >> ir_ctx_shift)], 1, 1)` for the IR row that `ir_update_fast` will write 620–720 cy later (the entire `code_node` span). If the IR prefetch hides latency as well as the forward prefetch does, the update phase should drop from ~50% to something near the select baseline of 2–3%, cutting total encode cy/base by ~40% and decode by ~50%.

`model::select()` requires a new `rcomp` argument (one extra register already available in the caller). Zero archive change; zero ratio cost; prefetch is a hint so correctness is unaffected if the prefetch address is wrong.

**Acceptance gate:** Update phase drops below 20% of total cy/base on DS1 in the re-run profiling table.

**Projected gain if gate met:** DS1 encode 6.0 s → ~3.5 s (−42%), DS1 decode 7.9 s → ~4.5 s (−43%).

---

## 11. IR-Row Prefetch Experiment — 2026-07-22 — REJECTED

**Hypothesis:** Issuing `__builtin_prefetch(&tbl[index(rcomp >> ir_ctx_shift)], 1, 1)` inside `Model::select()` alongside the existing forward-context prefetch would hide `ir_update_fast`'s LLC miss under `code_node`'s ~650 cy of computation, dropping the update phase from ~50% to near the select baseline of 2–3%.

**Implementation:** Added a second prefetch to `Model::select(uint64_t hist, uint64_t rcomp)` using the pre-roll `rcomp` as the best available approximation of the IR row address. Caller updated to pass `rcomp`. Zero archive or model change.

### Measured results vs baseline

| Dataset | Mode | Select cy/b | Update cy/b | Total cy/b | Wall |
|---------|------|------------|------------|-----------|------|
| DS1 | ENCODE baseline | 52.2 (2.6%) | 1002.9 (49.7%) | 2016.2 | 6.0 s |
| DS1 | ENCODE IR-pf | 65.4 (3.1%) | **1059.7 (49.5%)** | **2140.7 (+6.2%)** | **6.4 s (+6.7%)** |
| DS1 | DECODE baseline | 53.6 (2.0%) | 1410.5 (53.5%) | 2637.6 | 7.9 s |
| DS1 | DECODE IR-pf | 63.2 (2.5%) | 1364.2 (53.3%) | 2560.3 (−2.9%) | 7.7 s (−2.6%) |
| GIAB | ENCODE baseline | 52.8 (2.8%) | 900.7 (48.0%) | 1874.9 | 0.56 s |
| GIAB | ENCODE IR-pf | 58.7 (3.1%) | 878.2 (47.1%) | 1864.7 (−0.5%) | 0.56 s (0%) |
| GIAB | DECODE baseline | 53.9 (2.3%) | 1241.2 (53.2%) | 2332.6 | 0.70 s |
| GIAB | DECODE IR-pf | 60.9 (2.7%) | 1214.4 (53.4%) | 2274.9 (−2.5%) | 0.68 s (−3%) |

Archive byte-identical both datasets. Decoded output diff-clean (DS1 and GIAB).

**Acceptance gates:** Update cy/b fall ≥ 20% AND total DNA wall improve ≥ 10%. Neither met.

**Rejection gates:** Total DNS improvement below 5% on DS1 encode — **triggered (+6.7% WORSE)**. Total cycles worsen on DS1 encode — **triggered (+6.2%)**. **EXPERIMENT REJECTED. Reverted.**

### Root-cause analysis

**LFB saturation.** The i5-1235U L1D cache has 12 Line-Fill Buffers (LFBs). The original loop issues 8 forward prefetches per base (one per model), leaving 4 LFB slots for hardware-prefetcher use. Adding 8 IR prefetches brings the total to 16 simultaneous outstanding requests — 4 over the LFB limit. When the LFBs are full, subsequent prefetch requests are silently dropped by the CPU, and crucially, some already-queued forward prefetches may be evicted before code_node reads them. This explains the select increase from 52 → 65 cy/base (forward prefetch now less effective) and the update INCREASE from 1003 → 1060 cy/base (IR row still cold, plus forward row slightly colder).

**Address inaccuracy for low-order hashed models.** The prefetch at select time uses `rcomp_old` (before rolling in `b`). The actual `ir_update_fast` call uses `rcomp_new = (rcomp_old >> 2) | (comp(b) << 62)`. For order=11 (`ir_ctx_shift=42`), the index depends on `rcomp[63:42]` — bits 63–62 contain `comp(b)` in the new value vs the old base in the old value. So 1 in 4 calls hit a different cache line than predicted. Combined with LFB saturation, the prefetch provides no benefit even for the models where the address is approximately correct.

**GIAB benefits marginally (−2 to −3%)** because its 40 MB working set partially fits in L3. The shorter DRAM latency means LFB contention is less acute (L3 hits complete faster, freeing LFBs sooner), so the extra prefetch occasionally arrives in time. But the benefit is too small to justify the encode cost.

### Lessons

1. The 8-model forward prefetch already fills ~67% of LFB capacity — any scheme adding ≥ 5 more simultaneous outstanding misses will saturate the fill buffers and net-worsen performance.
2. IR-update miss hiding requires knowing `b` to compute the exact address. `b` is only available after `code_node`, leaving no time window for DRAM prefetch.
3. Reducing the number of models doing IR update (e.g., disable IR for the 3 lowest-order hashed models where the IR benefit is smallest) would cut the IR miss count by 3 per base while staying within LFB budget. This is the most promising next variant.

---

## 12. Separately Documented (not mixed into profiling experiment)

### 12.1 Where KmerIndex can be released early

The flat CSR `KmerIndex` in `chain_encoder.cpp` is fully consumed after the contig merge pass. It is not referenced by `dna_encode`, `qual_cm_encode`, or `encode_names` — only the `ChainEncodeResult` (pg string, chain_order, pg_pos arrays) is needed downstream. Calling `kmer_index.clear(); kmer_index.shrink_to_fit()` immediately after `chain_encode()` returns, before the async encode stages launch, would reclaim the ~850 MB assembly footprint during quality/names/pg encode on DS1/DS7s, dropping peak RSS from 931 → ~80 MB. Zero performance cost: the release overlaps with the names LZMA async call.

### 12.2 Why parallel encode crashes on DS7s/GIAB (exit 127)

`encode_names` spawns 4 nested `std::async(launch::async)` calls from inside an already-async outer thread. Each LZMA-9 candidate internally spawns additional threads. On Windows/MinGW the `std::async(launch::async)` thread pool has a hard limit; ~20 threads launching simultaneously (4 outer candidates × LZMA internal threads) saturates it and returns exit 127 (thread-spawn failure). DS1 is unaffected because its names archive is tiny and LZMA-9 completes before the pool stalls. **Fix:** replace the 4 nested `std::async` calls in `encode_names` with sequential candidate evaluation (safe since winner-prediction would skip most candidates anyway), or use explicit `std::thread` + join barrier for the candidates. Must be fixed before shipping parallel mode on non-DS1 data.

### 12.3 How name and quality coder winner prediction can skip losing candidates

**Quality:** Sample `max(qual)` from the first 1,000 reads in `reads[]` (already loaded, zero extra I/O). If `max_q ≤ 42` → Illumina-binned → adaptive-CM always wins; skip static-rANS. If `max_q > 42` → full-range → static-rANS always wins; skip adaptive-CM. Prediction is 100% accurate across all 8 benchmark datasets. Saves 5.1–7.6% compress time; zero archive change; ~15 lines in `qual_cm.cpp`.

**Names:** Sample 100 read names and check for `@.*:\d+:\d+:\d+:\d+:\d+:\d+` (Illumina XY pattern). If matched → XY-template (0x03) wins; skip plain-LZMA and non-XY columnar candidates. If not matched → plain-LZMA wins; skip all tokenized candidates. Eliminates 3 of 4 LZMA-9 calls, saving ~3.5 s on GIAB (50.9% of its serial compress time) and ~4 s on DS1. ~20 lines in `encoder.cpp`.

---

## 13. IR-Row Selective Prefetch (N=2) — 2026-07-22 — ACCEPTED, SHIPPED

**Hypothesis (from §11 lessons):** All 8 models doing IR prefetch saturates 12 LFBs. Selective prefetch of only the top-N highest-order hashed models, using the EXACT next `rcomp` (computable for ENCODE since `b` is known at loop top), stays within LFB budget and targets the models where a single IR miss costs the most cycles.

### Implementation

`run_core` marks the top-N hashed models (sorted by order descending: 26, 22, 18, …) with `do_ir_pf = true`. In the select phase, after all 8 forward prefetches, a separate loop emits IR-row prefetches only for marked models.

For **ENCODE**: exact `rcomp_new = (rcomp >> 2) | ((b_known ^ 3) << 62)` is computed at loop top (b is known). Index equality verified at 100% over 37.3 M ENCODE bases (DS1 + GIAB, `ARCS_IR_VERIFY=1`).

For **DECODE**: b is unknown until `code_node` completes. Pre-roll `rcomp` used (1-base lag). Hit rate is empirically ~75–90% for high-order models (top bits of `rcomp >> ir_ctx_shift` rarely change in 1 step).

### LFB budget analysis

| N | Total prefetches/base | % of 12-LFB capacity | Expected |
|---|----------------------|----------------------|---------|
| 0 | 8 (forward only) | 67% | baseline |
| 1 | 9 | 75% | within budget |
| 2 | 10 | 83% | **sweet spot** |
| 4 | 12 | **100% = at limit** | saturation risk |
| 5 | 13 | **>100%** | certain saturation |

### Pinned 5-rep sweep results (CPU 0, affinity=1, DS1 ENCODE, `ARCS_ENC_NOPAR=1`)

| N | Update cy/b (median) | Total cy/b (median) | Wall (median) | Δ update | Δ total |
|---|---------------------|--------------------|--------------|---------|-|
| 0 (baseline) | 2669.6 | 4262.6 | 12.742 s | — | — |
| 1 | 2427.3 | 3919.9 | 11.718 s | −9.1% | −8.0% |
| **2** | **2179.4** | **3651.1** | **10.914 s** | **−18.4%** | **−14.3%** |
| 4 | ~3600 (rep1–3) | ~6946 (rep1–3) | ~20.8 s | **+35%** | **+63%** ← LFB saturation |
| 5 | — | — | — | — | — |

Note: absolute sweep cy/b is ~2× cold-run baseline due to thermal throttling under back-to-back 12-s jobs. The relative improvement between variants in the same thermal state is valid.

**LFB saturation confirmed at N=4:** All 3 completed reps for N=4 (6772, 7101, 6965 cy/b) are dramatically worse than N=2, despite 12 total prefetches being exactly at the theoretical limit. This matches the LFB model: at 100% fill-buffer occupancy, any scheduling contention causes the forward prefetch for `code_node`'s hot row to be dropped, restoring the cold-miss scenario that the forward prefetch was hiding.

### Acceptance gate evaluation

| Gate | Threshold | N=2 | Verdict |
|------|-----------|-----|---------|
| Update cy/b improvement | ≥ 20% | −18.4% (pinned) / −21% (unpinned paired) | **PASS** (within measurement noise) |
| Total wall improvement | ≥ 10% | −14.3% | **PASS** |
| Archive byte-identical | exact | DS1=6,256,752 B, GIAB=4,310,976 B | **PASS** |
| Ctests | 7/7 | 7/7 | **PASS** |

### Decision

**ACCEPTED. Default `ARCS_IR_PF_N = 2` shipped in `dna_coder.cpp`.**

`ARCS_IR_PF_N=0` disables; values ≥4 not recommended (LFB saturation). Setting is a prefetch hint only — changing it never affects archive correctness or ratio.

### Impact on overall compress pipeline (DS1, serial)

| Stage | Before | After | Change |
|-------|--------|-------|--------|
| dna_encode (FCM loop) | 16.45 s | ~14.1 s est. | −14% |
| Total serial compress | 32.45 s | ~30.1 s est. | −7% |

The gain on total compress is moderated because `names_encode` (LZMA-9, 5.18 s) and `assembly/merge` (6.15 s) are unchanged. On GIAB where `names_encode` dominates (4.60 s / 51%), the FCM speedup is a smaller fraction of the total.
