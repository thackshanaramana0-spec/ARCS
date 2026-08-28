# Headroom ledger — where ARCS's sequence bytes actually are

**Purpose: prevent layer tunnel-vision.** Every proposed optimisation gets
checked against this table BEFORE any code is written. If a layer is worth
<2% of `total_seq`, it does not get worked on, however interesting it is.

Measured 2026-08-28 at `57fe7df` on the three scope inputs
(`scope/make_scope_inputs.sh`, raw logs in `scope/results/`). One run each,
serial. Reproduce with `ARCS_VODBG_TIMING=1 arcs compress <in.fq> <out.arcs>`.

## Where the bytes are

| dataset | pg_blob | pos_blob | aux_blob | total_seq |
|---------|---------|----------|----------|-----------|
| yeast (S. cerevisiae, 1M reads) | 3,218,276 **71.6%** | 590,308 13.1% | 687,388 15.3% | 4,495,972 |
| E. coli (1M) | 2,235,689 **66.4%** | 454,152 13.5% | 676,364 20.1% | 3,366,205 |
| P. falciparum (500k) | 2,586,871 **87.3%** | 150,780 5.1% | 226,236 7.6% | 2,963,887 |

**pg_blob is 66-87% of everything.** pos and aux together are 13-34%. Any
layer that does not touch pg_blob is playing for at most a third of the total.

## Against PgRC2

| dataset | ARCS total_seq | PgRC2 archive | margin |
|---------|----------------|---------------|--------|
| yeast | 4,495,972 | 4,898,620 | **-8.22%** |
| E. coli | 3,366,205 | 3,779,401 | **-10.93%** |
| P. falciparum | 2,963,887 | 3,115,115 | **-4.85%** |

## The floor, and how much slack is left

`lit` = pg after repeat-elim. Genome sizes: yeast 12.1 Mbp, E. coli 4.64 Mbp,
P. falciparum 23.3 Mbp.

| dataset | pg/genome | lit/genome | repeat-elim matches | lit as % of pg |
|---------|-----------|------------|---------------------|----------------|
| yeast | 1.35x | **0.98x** | 72,771 | 72.5% |
| E. coli | 4.30x | **1.21x** | 205,714 | 28.2% |
| P. falciparum | 0.68x | 0.32x | 203,415 | 47.2% |

**`lit/genome` converges to ~1x.** That is the ledger's most important line.
After repeat-elim the pseudogenome is approximately one genome copy — which is
the information content actually present — so pg_blob's floor is roughly
`genome_size x bpb/8`:

| dataset | pg_blob | est. floor | slack |
|---------|---------|-----------|-------|
| yeast | 3,218,276 | 2,962,988 | **255,288** |
| E. coli | 2,235,689 | 1,131,983 | **1,103,706** |
| P. falciparum | 2,586,871 | (n/a) | (n/a) |

The floor estimate is an idealisation: it assumes the pg should be exactly one
genome copy, and ignores that a real pg must also carry sequencing error and
true polymorphism. Treat it as an order-of-magnitude bound, not a target. It
does not apply to P. falciparum at all — 500k reads cover 0.32x of a 23.3 Mbp
genome, so there is no genome copy to converge to and the formula returns a
floor above the actual size.

**Trap avoided.** Yeast's repeat-elim finds only 72,771 matches against E.
coli's 205,714, which reads like a defect worth fixing. It is not. Yeast's pg
is already 1.35x its genome before elimination; there is almost nothing to
remove. E. coli's is 4.30x, so elimination has a great deal to find. Weak
repeat-elim on yeast is correct behaviour, and chasing it would have been days
spent for nothing.

## Per-read costs

| dataset | pos bits/read | aux bits/read | mismatches/read | bits per mismatch |
|---------|---------------|---------------|-----------------|-------------------|
| yeast | 4.72 | 5.50 | 0.43 | **12.83** |
| E. coli | 3.63 | 5.41 | 0.54 | 9.94 |
| P. falciparum | 2.41 | 3.62 | 0.42 | 8.54 |

A mismatch entry needs a position in the read (log2(150) ~ 7.2 bits if
uniform) plus a base excluding the reference symbol (log2(3) = 1.58 bits) ~
8.8 bits. Real mismatch positions cluster toward read ends, so the true
entropy is below that. Yeast spends 12.83.

Reads are emitted in global pg-position order (`3cc53a4`), so pos is a delta
stream: 1M reads over a 16.3M pg is a mean gap of 16.3, ~4.03 bits if the gaps
were geometric. Yeast spends 4.72.

## Ranked headroom

| # | layer | est. yeast value | est. E. coli value | confidence |
|---|-------|------------------|--------------------|------------|
| 1 | pg_blob: residual assembly redundancy above 1x genome | ~255 KB | ~1,100 KB | medium |
| 2 | aux_blob: ~4 bits/mismatch above entropy | ~215 KB | ~75 KB | high |
| 3 | DNA coder: 1.96 bpb -> ~1.85 | ~180 KB | ~125 KB | medium |
| 4 | pos_blob: ~0.7 bits/read above a geometric-gap model | ~87 KB | ~0 KB | high |

Layers 2 and 4 are coding problems with known floors — high confidence, small.
Layers 1 and 3 are larger and less certain. Nothing here is worth more than
~25% of total_seq on its own.

## PgRC2 component scoring — what to take and what to refuse

From the PgRC2 paper (Bioinformatics 41(3), btaf101) and its source at v2.0.2.

| their component | ARCS equivalent | verdict |
|-----------------|-----------------|---------|
| 1 QualDivision — Phred at 0.12x read len from end, <=2 -> LQ | HQ-TRIAL sweep, picks frac minimising pg_len over 3 trials | **refuse** — theirs is a speed heuristic (1.8x), ours is data-driven and picks a different frac per dataset (0.15/0.45/0.65 observed) |
| 2 PgGenDivision | HQ frac selection | parity |
| 3 Pg(HQ) — greedy SCS, hash + sparse sampling, lock-free | SA/APSP exact via libsais | **refuse** — their own paper concedes hash collisions across threads cause missed matches; ours is exact. pg 16.35M vs their 22.60M |
| 4 ReadsMatching — allows floor(read_len/3) mismatches | fallback tolerance read_len/8 (`c9837c9`) | **candidate** — they are 2.7x more relaxed |
| 5 Pg(LQ&N) — second and third pseudogenomes | singleton append | **candidate**, but the ledger caps it: yeast's whole survivor cost lives inside pos+aux = 1.28 MB |
| 6 RC redundancy removal — SimplePgMatcher, hash + sparse sampling, serial | repeat_elim, SA/LPF over pg + RC(pg), exact, parallel | **refuse** — ours already handles RC *and* has a completeness guarantee theirs lacks |
| 7a VarLenDNACoder — 242-phrase dictionary, 1-4 bases/byte, then LZMA | ARCS-DNA context model, 1.78-1.96 bpb | **refuse** — their own paper calls it a speed optimisation with <0.1% ratio loss; it exists to shrink LZMA input 3x, not to compress better |
| 7b mismatch bases: order-0 rerank by frequency, then exclusion of the pg symbol's rank | aux v1 rank-of-3 (`7747505`) | **partial** — we have exclusion; the order-0 frequency rerank needs checking |
| 7c mismatch positions/bases: PPMd, range coder, FSE and LZMA all tried, best kept | quality stream picks between static-rANS and adaptive-CM; pos/aux do not bake off | **candidate** — cheapest real win on the board |

Three of their seven components are things ARCS should explicitly not adopt,
and in each case ARCS's existing approach is the stronger one. The genuine
candidates are 4, 5, 7b and 7c — and 7c is the only one that is pure coding
with no assembly risk.

## Anti-goals

- Do not chase stage 5 (their LQ/N pseudogenomes) as a size win. The ledger
  bounds the entire yeast survivor cost at 1.28 MB, and ARCS already wins
  overall on all three datasets.
- Do not port VarLenDNACoder. It is a speed trick and would cost ratio.
- Do not "fix" yeast's repeat-elim match count. See the trap above.
- Do not optimise a layer worth <2% of total_seq.

---

# Layer 7c — measured: per-stream coder selection

`ARCS_CODER_PROBE=1` (print-only, alters no blob; `src/encoder.cpp`). Both
streams have only ever been through LZMA-9. Tried libbsc (BWT+QLFC), already
in-tree, against it.

| dataset | aux single-LZMA | aux single-BSC | aux per-column best-of | pos LZMA | pos BSC |
|---------|-----------------|----------------|------------------------|----------|---------|
| yeast | 688,380 | 673,162 (-2.21%) | **652,822 (-5.17%)** | 590,140 | 593,316 (+0.54%) |
| E. coli | 676,440 | 650,932 (-3.77%) | **631,594 (-6.63%)** | 453,936 | 452,364 (-0.35%) |
| P. falciparum | 227,016 | 214,020 (-5.72%) | **204,800 (-9.79%)** | 150,892 | 144,752 (-4.07%) |

BSC wins 9 of 11 aux columns on every dataset. The two it loses are a 44-byte
header and an empty column, both of which store raw anyway.

**Effect on total_seq:** yeast -35,558 (**-0.79%**), E. coli -46,418
(**-1.38%**), P. falciparum -28,356 (**-0.96%**).

Container framing (~90 B/column x 11 = ~1 KB) is what killed the earlier
separation-only experiment. It does not kill this one: BSC's gain is 22-45 KB,
which dwarfs it.

**Time cost: none measurable.** yeast compress 19.45 s baseline vs 19.26 s with
all 24 extra coder passes — the probe run was faster, so the difference is
noise. The `ARCS_AUX_SPLIT` note's "1.03s for ten LZMA-9 passes" does not
reproduce.

## Verdict against this ledger's own bar

The bar in Anti-goals is 2% of total_seq. **This layer is worth 0.79-1.38%, so
it fails the bar.** Recording that plainly rather than rationalising it: BSC on
these streams is a real, free, universal ~1%, but it is not where the paper's
argument lives, and it must not be allowed to consume the attention that
Layer 1 deserves.

Recommended anyway, because the bar exists to ration *effort*, and here the
effort is nearly spent already: the 0xA2 per-column container exists, the codec
is in-tree, and keep-the-smaller makes regression structurally impossible. It
needs a per-column coder tag and decoder support, and nothing else.

**Scale check before proceeding:** Layer 1 is worth ~1,100 KB on E. coli.
Layer 7c is worth ~46 KB on E. coli. Layer 1 is 24x larger. If attention has to
go one place, it goes there.
