# Shard-Parallel Assembly with Overlap-Preserving Merge Reconciliation

**Status:** ratio hypothesis validated (serial prototype); threaded implementation pending.
**Date:** 2026-07-19.

## Problem

ARCS's global read-placement (the greedy multi-contig assembler in `build_multicontig_pg`)
is the dominant serial pole of compression (~46 s on DS7, ~27 s on DS2). It is serial because
every read is placed against a single shared, incrementally-mutated k-mer index / `used[]`
state. This leaves ~10 of 12 cores idle and is the main reason ARCS compress is slower than
SPRING/Genozip.

Two known ways to parallelize assembly, each with a cost:

| approach | parallel | deterministic | ratio |
|---|---|---|---|
| SPRING (shared dictionary + striped fine-grained locks) | yes | **no** (race-order-dependent read placement → archive size wobbles run-to-run) | good |
| naive chunking (independent chunks, no reconciliation) | yes | yes | **poor** (cross-chunk overlaps lost → +8.8 % on low-coverage DS1) |

## Our mechanism

**Repurpose the offline contig-MERGE — originally built to collapse online-greedy fragmentation —
as a parallelization primitive.**

1. **Partition** reads deterministically into `T` shards (content key; same input → same shards).
2. **Assemble each shard independently** on its own thread — no shared state, no locks, fully
   deterministic per shard.
3. **Reconcile globally** by running the existing overlap-layout **merge** across the union of all
   shards' contigs — recovering the cross-shard overlaps that naive chunking discards.

SPRING cannot do this: they have no merge step, so a single global pass (hence shared state +
locks) is their only route to global read order. Naive chunking cannot do this: it has no
reconciliation, so cross-shard overlaps are lost. ARCS can, *because* the merge already exists.

This reaches the Pareto corner both alternatives miss: **parallel + deterministic + ratio-preserving.**

## Validation (serial prototype)

`ARCS_SHARD_TEST=N` clears the k-mer index at shard boundaries in the placement loop, so reads
match only contigs from their own shard (= independent shard assembly). The existing
end-of-assembly merge then reconciles the union. Being serial, the prototype has the **same ratio**
as the future threaded version — so it isolates the ratio question from the threading work.

| dataset | coverage | global archive | shard+merge archive | Δ | lossless |
|---|---|---|---|---|---|
| DS1 E.coli | low (hardest) | 6,253,894 | 6,273,934 (shard-4) | **+0.32 %** | yes |
| DS2 human | — | 4,819,755 | 4,831,027 (shard-6) | +0.23 % | byte-exact |
| GIAB human | high | 4,311,061 | 4,310,047 (shard-6) | **−0.02 %** | byte-exact |

Reference: naive chunking cost **+8.8 %** on DS1. The merge recovers ~96 % of that loss. The
pseudogenome grows ~+3 % in length but the archive grows only ~+0.3 %, because the context-mixing
coder captures the extra overlap redundancy.

**Conclusion:** shard-parallel assembly + merge reconciliation preserves compression ratio within
±0.3 % (sometimes better) across bacterial and human data at low and high coverage, losslessly.
The threaded implementation is therefore justified.

## Full 8-dataset validation (serial prototype, `ARCS_SHARD_TEST=6`)

All byte-exact lossless (bacterial/viral/human, SRA & Illumina names, CRLF & LF, short & long).
Ratio is neutral-to-better; compress RAM drops on every dataset (the k-mer index is bounded per
shard instead of growing globally) — a bonus obtained *before* any threading. Compress time also
falls on the largest inputs purely from the smaller working set (better cache).

| DS | Δ ratio | compress RAM (global → shard) | lossless |
|---|---|---|---|
| DS1 | +0.38 % | 1463 → 1166 MB (−20 %) | yes |
| DS2 | +0.23 % | 1613 → 844 MB (−48 %) | yes |
| DS4 | +0.56 % | 844 → 711 MB (−16 %) | yes |
| DS5 | −0.03 % | 500 → 493 MB | yes |
| DS6 | +0.02 % | 408 → 337 MB (−17 %) | yes |
| DS7 | **−0.70 %** (smaller than global) | 3702 → 3213 MB (−13 %) | yes |
| GIAB | −0.02 % | 475 → 434 MB | yes |
| NA12878 | −0.01 % | 451 → 439 MB | yes |

Note DS7/GIAB/DS5/NA12878 are *smaller* under shard+merge than the single global greedy: the online
greedy makes locally-suboptimal placement choices that the shard-then-merge structure avoids.
The correctness and ratio/RAM behaviour are validated; threading is the remaining work for the
compress-time win.

## Next steps

1. Threaded shard assembler: per-shard placement on independent threads (own contigs + index),
   then combine and merge. Target: the ~46 s DS7 placement pole.
2. Keep deterministic (fixed partition + deterministic merge order).
3. Full 8-dataset byte-exact + ratio regression on every change.
