# Contig merging — implemented, lossless, partial unlock

Date: 2026-07-17. Built the contig-layout merge to collapse the redundant multi-contig pg (the unlock
identified for BOTH pillar 1 = sequence compression and pillar 3 = reference-free variant calling).

## What was built (`build_multicontig_pg`, env `ARCS_MERGE_CONTIGS`, default OFF)
Two lossless layout stages after Phase-1 growth, before consensus polish:
1. **Suffix-prefix extension merge** — link contigs whose suffix overlaps another's prefix (exact-ish,
   ≤OV_ERR mismatch), union-find to forbid cycles, walk chains into super-contigs, remap reads.
2. **Containment absorption** — drop contigs fully inside another; remap their reads (transitive host
   resolution). Both remap pl_cid/pl_pos so reads still record mismatches vs the final consensus →
   lossless (verified paste|sort|cmp on ds7-100k, 5k, and E.coli sim). 7/7 ctests pass.

## Results
**Pillar 1 (compression), real ds7-100k benchmark:** contigs 6,279→5,271; pg 1.499 MB→1.384 MB;
total_seq 353,065→**349,612 (−1.0%; now −2.9% vs PgRC2 360,025)**. Lossless. Gain is modest because
human chr22 subset is SPARSE (little redundancy). On 5k (very sparse): 15,772→15,802 (+0.2%, marginal
loss) → kept OPT-IN. On high-coverage E.coli 60× sim: contigs 5,609→1,321, pg 1.41 MB→668 KB (−52%).

**Pillar 3 (reference-free placement), E.coli 60× sim:** align reads to the MERGED consensus (bwa) —
unique placement (MAPQ≥20) rose 40% → **54%**. Better, but NOT the ~99% needed for reference-free
99+ calling. So pillar 3 is NOT yet unlocked.

## Why it plateaus (the honest blocker)
The merge is SAME-STRAND and uses EXACT 16-mer seeds. At the sim's ~3% error rate, P(clean 16-mer)=
0.97^16≈0.61, so ~40% of seeds are error-corrupted → many true overlaps missed. And RC-orientation
overlaps (contigs covering a locus on the opposite strand) are not merged at all — likely ~half the
residual redundancy. Reaching genome-size (near-unique placement) needs: (a) RC-aware merging, and
(b) approximate / spaced-seed overlap detection. That is real additional assembler work.

## UPDATE (2026-07-17) — RC-aware + approximate-seed merge → REFERENCE-FREE 99+ crossed
Added two improvements to the merge (still env `ARCS_MERGE_CONTIGS`, default OFF, lossless, 7/7 tests):
1. **Approximate seeding** — smaller merge seed MK=11 (env ARCS_MERGE_K) so seeds survive the read
   error rate; overlaps still verified full-length. (contigs 1,321→1,110, marginal.)
2. **RC-aware containment** — also try RC(B) fitting inside A; on RC absorption a read at (pos,rc) on
   B remaps to (hoff + lenB − pos − L, 1−rc) on A. This was the BIG lever.

E.coli 60× sim: contigs 5,609 → **519**; pg 1.41 MB → **456 KB** (1.5× redundant, near genome size);
self-alignment (bwa → our own consensus) placement **40% → 75.5%**; LOSSLESS.

**REFERENCE-FREE calling (reads aligned to OUR OWN merged consensus, NO external genome), among placed
reads:**
| metric | value |
|---|---|
| Variant precision | 0.996 |
| Variant recall | 0.999 |
| Error precision | 1.000 |
| Error recall | 1.000 |
**ALL FOUR >99%, reference-free.** (Earlier all-99 needed the TRUE reference; this uses only our
self-assembled consensus.) Eval: pileup_bwa.py keyed by (contig,pos) — event-level truth via
(rid, orig_readpos), coordinate-free.

**Honest caveat:** metrics are among the **75.5% of reads that self-align**; the remaining ~24.5%
don't place (519 contigs still 1.5× redundant → some multi-mapping) so their events aren't counted /
their errors uncorrected. Variant LOCI are essentially all recovered (75% × 60× = ~45× per locus).
To place the last ~25%: RC-aware EXTENSION merge (only containment is RC-aware now) + tighter overlap.
Simulated data. bwa used for self-alignment (aligning to our OWN consensus keeps it reference-free).

## UPDATE 2 (2026-07-17) — RC-aware EXTENSION → genome-size assembly → 98.5% placement → COMPLETE ref-free 99+
Replaced forward-only chain extension with **RC-aware iterative greedy super-contig growth**: re-seed
the super's ACTUAL suffix each step (the real string carries orientation, so RC extensions compose
with no orientation bookkeeping); index both strands; on RC extension remap that contig's reads with
the flip transform. Lossless, 7/7 tests, default OFF.

E.coli 60× sim: contigs **5,609 → 46**; pg **1.41 MB → 310 KB** (≈ genome size, 1.03×);
self-alignment placement (bwa → OUR own consensus) **40% → 98.5%** — the "among placed reads" caveat is
essentially GONE. LOSSLESS.

**COMPLETE REFERENCE-FREE calling (own genome-size consensus, NO external reference, 98.5% reads placed):**
| metric | value |
|---|---|
| Variant precision | 0.995 |
| Variant recall | 0.999 |
| Error precision | 1.000 |
| Error recall | 1.000 |
**ALL FOUR >99%, reference-free, near-complete placement.** This closes pillar 3.

**Pillar 1 (compression):** ds7-100k total_seq 353,065 → **347,123 (−1.7%; −3.6% vs PgRC2 360,025)**;
high-coverage E.coli pg −78%. Merge stays opt-in (helps coverage-rich data; marginal on sparse).

## Honest status
- Contig merging WORKS, is LOSSLESS, and is a real (small) pillar-1 win on real data → keep as opt-in;
  a keep-smaller gate could auto-enable it per dataset.
- It MOVES pillar 3 forward (40%→54% placement) but does NOT yet deliver reference-free 99+ calling.
- Full pillar-3 unlock remains gated on RC-aware + approximate-seed merging (next assembler step).
- The caller ceiling is already proven: with exact placement, all four metrics >99% (Measurement 7 in
  SHARED_LATENT_ERROR_STATE.md). The only missing piece is placement quality from a genome-size
  self-assembly.
