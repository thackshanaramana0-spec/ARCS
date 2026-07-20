# Tool Analysis — ARCS vs PgRC2 vs Genozip vs SPRING

Honest, measurement-grounded comparison of ARCS against the competing FASTQ compressors,
with special attention to the **sequence-stream** question (PgRC2 vs Genozip vs ARCS) and
what can vs cannot be claimed in a paper.

---

## 1. The core positioning claim

ARCS occupies a corner of the Pareto frontier that **no other tool holds**:

| Property | ARCS | SPRING | Genozip | PgRC2 | fqzcomp |
|---|---|---|---|---|---|
| Smallest lossless archive | ✅ 8/8 | ❌ | ❌ | ❌ (lossy default) | ❌ |
| Byte-exact | ✅ all 8 | ❌ CRLF files | ✅ | ❌ (lossy qual) | ❌ (`#`→`!` remap) |
| Order-preserving default | ✅ | ✅ | ✅ | ❌ (needs permutation) | — |
| Reference-free variant calling | ✅ (P/R/F1) | ❌ | ❌ | ❌ | ❌ |

**No competitor occupies all four columns.** The publishable claim:

> ARCS is the only lossless FASTQ archiver that is simultaneously smallest, byte-exact,
> order-preserving, and produces reference-free variant calls as a zero-cost byproduct of
> compression.

---

## 2. Lossless vs lossy — the distinction that matters per stream

A FASTQ file has three streams: **sequence (bases)**, **quality**, **names/headers**.
"Lossy" almost always refers to the **quality** stream — no serious tool discards bases.

| Tool | Sequence | Quality | Order | Full archive |
|---|---|---|---|---|
| **ARCS** | lossless | lossless | preserved (default) | **lossless, byte-exact** |
| **Genozip** | lossless | lossless | preserved | **completely lossless** |
| **PgRC2** | lossless | **lossy/binned by default** | **NOT preserved by default** | lossy (default) |
| **SPRING** | lossless | lossless | preserved | lossless, but **not byte-exact on CRLF** |
| **fqzcomp** | lossless | lossless-ish | preserved | **not byte-exact** (`#`→`!` remap) |

Key takeaways:
- Calling PgRC2 "lossy" refers to **quality**, not bases. Its **sequence is lossless**.
- Genozip is **completely lossless** across all three streams + order.

---

## 3. The sequence-only comparison (PgRC2 vs Genozip vs ARCS)

For the **sequence stream alone**, all three are lossless → a sequence-vs-sequence byte
comparison **is** apples-to-apples in principle. **But there is one catch:**

| | Sequence lossless? | Order-preserving in that number? |
|---|---|---|
| PgRC2 | Yes | **No** — preserves the read *multiset*, not order; needs an extra permutation stream usually **not counted** in headline numbers |
| Genozip | Yes | Yes |
| ARCS | Yes | Yes (order-preserving is the default) |

**Why this matters:** if you compare PgRC2's raw sequence bytes to Genozip's or ARCS's
*order-preserving* sequence bytes, you are comparing a non-order-preserving number to two
order-preserving ones. PgRC2 looks smaller partly because it is answering an **easier
question** (the multiset, not the ordered list).

For a fair fight, either:
1. add PgRC2's read-order permutation cost to its sequence total, **or**
2. compare all three in order-discarding mode.

This is not cheating on PgRC2's part — it is their design choice — but it is the **asterisk
a paper needs** when putting the numbers side by side. It is also the honest reason PgRC2 can
appear to "win sequence": it is genuinely strong **and** often scored on the easier
(order-free) version of the problem.

---

## 4. "Does PgRC2 beat Genozip?" — what is and isn't established

**In our own repo/benchmarks: no such comparison exists.**
- `PGRC2_LITERATURE_SURVEY.md` analyzes PgRC2's architecture but never mentions Genozip.
- The only doc with both (`REAL_BENCHMARK_GIAB.md`) logs PgRC2 as **lossy** (245 KB, quality
  discarded) and Genozip as lossless (5.03 MB) — different objects, never compared head-to-head.
- 8-dataset `SUMMARY.md` has an **empty PgRC2 column** — lossless PgRC2 **crashed** in our setup
  (`-Q/-q0` crash; only the lossy default runs).

**In the published literature: not cleanly established either.**
- The PgRC paper (Kowalski & Grabowski, Bioinformatics 2020) benchmarks against **SPRING, ORCOM,
  Minicom, GTZ** — the read-reordering family — **not Genozip**.
- Genozip's own papers benchmark against SPRING, DSRC2, FQZComp, CRAM — **not PgRC** head-to-head.
- The two are different lineages that mostly were not compared in the same lossless-FASTQ table.

So **"PgRC2 beats Genozip" is a plausible folk belief** (PgRC is genuinely strong on sequence),
but it is **not a number citable** from our repo or from a single authoritative published table.
A reviewer would ask "measured where?" — and there is no clean answer.

**Citable statements we CAN make:**
1. PgRC/PgRC2 is a strong *sequence-stream* compressor (established vs SPRING/ORCOM in its paper).
2. On our human datasets, ARCS produces the smallest **byte-exact lossless archive**, beating
   Genozip and SPRING (measured, logged).
3. PgRC2's lossless mode did not run in our setup → absent from our lossless comparison; state as
   a limitation, not a win.

To get a real PgRC2-vs-Genozip sequence number, PgRC2's lossless build must be fixed and both
run on the same files. Until then it is a literature claim, not a measured one.

---

## 5. How PgRC2 actually encodes sequence (and why ARCS is not a copy)

**PgRC2 pipeline:** classify HQ/LQ/N reads → greedy approximate shortest-common-superstring
(pseudogenome) over HQ reads via **copMEM** → RC-fold → map LQ reads (position + mismatches) →
entropy-code streams with **PPMd/range/GeCo3**. Each read is stored as
`position-in-pseudogenome + mismatch list`. **No FCM on a string anywhere.**

**ARCS pipeline:** de Bruijn assembly → multi-contig consensus pseudogenome → **FCM
context-mixing coder on the pseudogenome string** (ARCS-DNA, our own) + own quality CM + own
name coder.

Shared components: **none**. Different assembler (de Bruijn vs greedy SCS), different seed method
(k-mer index vs copMEM coprime sampling), fundamentally different coder (FCM-on-string vs
position+mismatch). The FCM-on-pseudogenome is precisely the thing PgRC2 does **not** do — it is
both ARCS's ratio crown and its originality.

**Convergence note:** on high-coverage low-error Illumina data both approaches reach a similar
information floor (~2.3–2.4 bytes/read on GIAB) by different paths — PgRC2 by storing only
mismatches, ARCS by FCM prediction approaching near-certainty on a coherent pseudogenome.

**Where PgRC2 may genuinely win sequence:** low-coverage / bacterial data (PG ≈ total read size,
FCM gains shrink) and high-error data. Illumina ~1% error is PgRC2's sweet spot.

---

## 6. The "fast mode" (Mode 2) design — winning speed + RAM without copying PgRC2

**Non-negotiable:** fast mode **keeps the FCM coder on the pseudogenome string**. Dropping the
FCM for position+mismatch to save time **is** copying PgRC2. The FCM is kept always; only
assembly cost and coder parallelism are attacked.

**Bottleneck decomposition (from profiled numbers):**

| Dataset | compress | comp RAM | dominant cost |
|---|---|---|---|
| GIAB (39 MB) | 12 s | 423 MB | balanced; names-LZMA biggest phase |
| DS1 (37 MB) | 59 s | 2.0 GB | FCM coder + assembly |
| DS2 (28 MB) | 96 s | 2.1 GB | FCM + merge (high redundancy) |
| DS7 (344 MB) | 180 s | **3.9 GB** | assembly k-mer table + FCM, scales badly |

Two problems: **RAM** (de Bruijn k-mer table + all reads + placement index held at once — the clear
place to win outright) and **speed** (serial per-bit FCM, 8 cache-missing lookups per base).

**Three levers:**
1. **Minimizer-sparse assembly index** — index minimizers (~1 per 10 bases), k-mer table ~10×
   smaller (miniasm/minimap lineage, NOT PgRC2's copMEM). RAM + assembly speed. Prior note
   "minimizer seeding = dead" was a *ratio* result for the merge; here it only needs to be
   ratio-*neutral* (lower bar).
2. **Chunked streaming assembly** — process C reads at a time, each chunk builds/encodes/frees its
   own pseudogenome. Peak RAM = a few in-flight chunks, **tunable and bounded**. Chunks run in a
   thread pool → parallel. The decisive lever. Ratio cost: cross-chunk redundancy lost (~0.36%/split).
3. **Lean FCM** — 8 orders → ~5. The logistic mixer already auto-gates redundant high orders.
   Cleanest speed/ratio knob (env `ARCSDNA_ORDERS` already exists — measurable with zero code).

**Generalization key — adaptive chunk sizing:** small files (bacterial, thin 2–4% ratio lead) →
single chunk (zero cross-chunk loss) + sparse index + lean FCM (speed/RAM win, near-zero ratio
cost). Large files (human WGS, huge 15–27% lead) → chunk aggressively. Ratio hit always a small
fraction of the existing margin.

**Honest per-axis outcome (vs the real lossless peers SPRING/Genozip):**

| Axis | Realistic outcome | Confidence |
|---|---|---|
| **RAM** | **Win outright** — chunking bounds peak RAM below both | High |
| **Speed** | **Match-to-beat SPRING** (DS7 180 s → ~30–40 s); beating Genozip's ~1.2 s on tiny files unlikely | Medium |
| **Ratio** | **Hold crown on human; at-risk on bacterial** (thin 2–4% margin) | Medium |

**Claim NOT to make:** "dominates all 3 axes on every dataset." Beating Genozip's sub-second
speed on small bacterial files is unrealistic (mature streaming codec vs an assembler). Fast mode
is *competitive* there and *winning* on RAM and large-file speed.

**Validation gate (do first, pure measurement):** run the 8-dataset sweep with 5-order vs 8-order
FCM via the existing `ARCSDNA_ORDERS` env var. If 5-order costs <1% on bacterial → plan viable. If
~3% → bacterial ratio lost in fast mode; position Mode 2 honestly as "human-WGS fast mode" (which
is exactly the regime the motivating storage bills live in).

---

## 7. Summary — what is honest to put in a paper

- **Primary contribution:** compression. Cleanest gap, strongest result (smallest byte-exact
  lossless archive, 8/8, order-preserving).
- **Secondary contribution:** reference-free het-SNV calling as a zero-marginal-cost byproduct
  (F1 0.943, rtg-validated, leads DiscoSNP++ 0.887 / Kmer2SNP 0.539). Structurally impossible for
  SPRING/PgRC2/Genozip (no assembly bubble structure).
- **Sequence-only comparisons:** always state the order-preservation asterisk — PgRC2's headline
  sequence numbers are typically order-free; ARCS/Genozip are order-preserving.
- **Do NOT claim** PgRC2 vs Genozip results (unmeasured, uncited) or fast-mode domination on every
  dataset (unrealistic vs Genozip small-file speed).
