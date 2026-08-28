# ARCS vs PgRC2, SPRING, Genozip — scope-matched comparison

Supplementary to Claim 1. **Does not replace it**: Claim 1 is whole-FASTQ
lossless compression of the 10 locked accessions (`DATASET_LOCKED.md`). This
document exists because PgRC2 does not compress whole FASTQ at all, so a
side-by-side needs a stated scope before any number means anything.

Measured 2026-08-28 on the SDC3 server (12 vCPU, ~90 GB). Every timed run was
run alone — no two concurrent.

## What each tool actually stores

| tool | sequences | qualities | names | original order |
|------|-----------|-----------|-------|----------------|
| ARCS | yes | yes | yes | yes |
| PgRC2 | yes | no (used only to split HQ/LQ reads, never stored) | no | only with `-o` |
| SPRING `--no-quality --no-ids` | yes | no | no | yes |
| Genozip | yes | yes | yes | yes |

PgRC2's `-q` is a quality-based read *division* threshold, not a quality
stream; its decompressor emits sequences one per line and nothing else.

So the only scope all four share is **DNA sequences**. Compare ARCS's
`total_seq` (pg_blob + pos_blob + aux_blob, printed by the encoder) against
PgRC2's whole archive. Comparing whole archive to whole archive measures
ARCS storing three streams the others do not have.

## Inputs

`scope/make_scope_inputs.sh` rewrites a real FASTQ keeping sequences but
replacing names with `@1..@N` and qualities with all-`I`, so names and
qualities cost ~nothing for every tool (ARCS quality stream: 128 B).

**Caveat, and it matters.** This transformation is not neutral for ARCS.
Assembly reorders reads, and on real FASTQ ARCS recovers the original order
for free from the numeric field of the read name — the encoder prints
`order-preserving: derive from name index col N (free)`. Flattening names
destroys that, forcing ARCS to store an explicit permutation: **2,568,957 B**
on 1M yeast reads, which is close to the information-theoretic floor for a
1M-element permutation (log2(1000000!) ~ 2,311,138 B) and is simply not paid
on real data. Any whole-archive number taken from these inputs therefore
charges ARCS ~2.5 MB it does not pay in practice.

## Sequence scope (the comparable number)

ARCS `total_seq`, against PgRC2's default-mode archive (mean of 3 runs):

| dataset | ARCS | PgRC2 | margin |
|---------|------|-------|--------|
| E. coli (SRR2584863, 1M reads) | 3,366,583 | 3,779,401 | **-10.92%** |
| P. falciparum (SRR37283774, 500k) | 2,963,614 | 3,115,115 | **-4.86%** |
| S. cerevisiae (DRR976266, 1M) | 4,495,009 | 4,898,620 | -8.24% (see below) |

## PgRC2 loses reads on yeast, and the loss varies run to run

Its default (no `-o`) mode is **not reliably lossless**. Same input, same
binary, repeated runs — decompressed read counts against 1,000,000 in:

    999,118  |  999,118  |  998,250  |  997,382  |  999,986  |  1,000,000

Between 0 and 2,618 reads silently discarded (up to 0.26%). Only 14 of the
882 lost in one inspected run contained `N`; the rest were low-complexity /
homopolymer-rich reads. With `-o` it was EXACT in 9/9 runs across all three
datasets, so the defect is specific to default mode, and among these three
datasets only to yeast — E. coli and P. falciparum preserved all reads in 3/3
default-mode runs each.

**Consequence:** the yeast sequence-scope margin above is not like-for-like.
PgRC2 encoded 882-1,750 fewer reads than ARCS did in every run. Do not quote
-8.24% without this caveat.

## Archive size is nondeterministic for both tools — PgRC2 more so

Spread across 3 runs of the same input:

| tool / mode | spread |
|-------------|--------|
| PgRC2 E. coli, default | 12,479 B (**0.33%**) |
| PgRC2 E. coli, `-o` | 2,326 B (0.038%) |
| PgRC2 P. falciparum, `-o` | 944 B (0.022%) |
| ARCS yeast (real input) | 539 B (**0.004%**) |

ARCS's own nondeterminism comes from the growth loop claiming reads by atomic
compare-exchange; `ARCS_DETERMINISTIC=1` makes archives byte-identical at a
wall-clock cost (see `src/vodbg_pg.cpp`). Every ARCS run is fully lossless
regardless. PgRC2's varies what data survives.

## Whole archive on scope-matched inputs (read with the caveat above)

**STALE RSS WARNING (2026-08-28):** the peak-RSS column below was measured at
`4c29ee0`, which predates `49113e0` and `6cd1e82` -- two commits that each cut
assembly peak memory roughly in half with byte-identical output. Every ARCS RSS
figure here, and the "~19x worse on peak RSS" conclusion at the end, therefore
overstates current ARCS. Re-measure before quoting. Size and time columns are
unaffected. Raw logs for this table are retained in `scope/results/`.


| dataset | tool | bytes | compress | decompress | peak RSS |
|---------|------|-------|----------|------------|----------|
| yeast | ARCS | 7,066,356 | 19.34 s | 7.01 s | 4.34 GB |
| | ARCS (pre-Method-B, 137911a) | 10,943,983 | 95.54 s | 31.47 s | 12.6 GB |
| | PgRC2 `-o` | 7,063,554 | 4.03 s | 0.28 s | 232 MB |
| | SPRING | 7,240,146 | 3.04 s | 0.40 s | 707 MB |
| | Genozip | 32,002,997 | 1.25 s | 0.51 s | 885 MB |
| E. coli | ARCS | 5,997,369 | 21.38 s | 5.47 s | 4.32 GB |
| | ARCS (137911a) | 7,744,034 | 59.15 s | 22.02 s | 12.6 GB |
| | PgRC2 `-o` | 6,148,738 | 2.92 s | 0.22 s | 218 MB |
| | SPRING | 7,243,152 | 4.23 s | 0.48 s | 707 MB |
| | Genozip | 33,082,142 | 1.21 s | 0.53 s | 878 MB |
| P. falc. | ARCS | 4,040,397 | 9.27 s | 3.94 s | 1.90 GB |
| | ARCS (137911a) | 4,680,238 | 28.23 s | 11.17 s | 4.26 GB |
| | PgRC2 `-o` | 4,269,147 | 1.54 s | 0.18 s | 193 MB |
| | SPRING | 4,276,752 | 1.88 s | 0.37 s | 501 MB |
| | Genozip | 7,005,009 | 0.77 s | 0.41 s | 514 MB |

SPRING's archive is a POSIX tar that pads every member; the figures above are
the summed member payload. The raw file size reads 7,290,880 for BOTH yeast
and E. coli — a padding artifact, not a real tie. Always sum the members.

Where we stand honestly: ARCS beats SPRING and Genozip on size everywhere, and
beats PgRC2 on E. coli and P. falciparum even while carrying a permutation
PgRC2 does not store. ARCS loses on decompression by ~25x and on peak RSS by
~19x; both are the direct cost of the exact all-pairs-suffix-prefix overlap
that produces the size win, and both belong in T2 as measured.

## Reproducing

    # PgRC2 (GPL-3, not vendored here — see .gitignore)
    git clone https://github.com/kowallus/PgRC.git method_c
    cd method_c && mkdir build && cd build && cmake .. && make PgRC
    # pinned at v2.0.2 (c2d85e4); binary reports
    # "PgRC 2.0: Copyright (c) Tomasz Kowalski, Szymon Grabowski: 2024-11-20"

    bash benchmark/scope/make_scope_inputs.sh reads.fq scope/reads.fq
    bash benchmark/scope/bench.sh          # 5-tool table, serial
    bash benchmark/scope/pgrc_verify.sh    # PgRC2 losslessness across repeats

SPRING's `-g` flag means "input/output is gzipped" — do NOT pass it with plain
`.fq` input, or compression fails. With `--no-quality` its decompressed output
is 2-line records (id, sequence), not 4-line FASTQ, so verify with `NR%2==0`.
