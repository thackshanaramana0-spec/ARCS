# Adversarial edge-case hardening — "lossless" is now a guarantee, not a GIAB observation

Date: 2026-07-17. Stress-tested chain-pg against 12 adversarial inputs (the inputs that break
"lossless" claims). Initial result: **6/12 FAILED**, including one SILENT losslessness violation. After
fixes: **12/12 pass**, GIAB baseline unchanged (4,280,170 B, still lossless, ratio crown intact), 7/7 ctests.

## Bugs found and fixed
| input | symptom | root cause | fix |
|---|---|---|---|
| high-Phred quality | **SILENT data loss** | quality clamped to Phred 42; legal FASTQ is 0–93 | QA 48→94; clamp 42→93; detect qhi>42 → skip 42-cap static path, force CM |
| 1 bp / sub-K / empty reads | **segfault** | `build_*_pg` bailed `if L<DEDUP_K` on reads[0] length → empty result crashed downstream | removed the bail; per-read lengths already flow through record_append |
| 0-length reads[0] | **quality corruption** | (a) qmax scanned `j<L` with global L=0 → qmax=0 → all quality→'!'; (b) encoder used global L=reads[0]=0 but decoder clamped L≤0→150 → posbin desync | qmax scans per-read length; global L = MAX read length (non-degenerate, identical both sides) |
| >65535 bp reads | **decode crash / corruption** | mm/N offsets stored uint16 (max 65535) | clear compile-time-domain guard: refuse with a helpful message (chain-pg is a short-read codec) |
| lowercase / IUPAC / any non-ACGT | **decode crash + case loss** | non-ACGT masked to 'A' then decoded as literal 'N'; `encode_base` case-insensitive so lowercase acgt silently upcased | store the LITERAL byte at each non-ACGT position (new `pg_N_char_flat` stream, aux format 36→40 B, backward-compatible); pg/contigs use strict-uppercase-ACGT test (`is_acgt_strict`) so lowercase/IUPAC route to the literal stream |
| exceptions | **std::terminate / abort** | no top-level catch in main | wrapped main in try/catch → clean "arcs: error: …" + exit 1 |

## What now round-trips losslessly (verified byte-identical)
single read · all-N · homopolymer · 1 bp reads · sub-seed (<16 bp) reads · 0-length reads · mixed extreme
lengths (1–300) · lowercase acgt · IUPAC (RYSWKM…) · full Phred range (0–93) · names with control chars /
embedded \x01 · Illumina-lookalike name traps · reads >65535 bp (clean refusal, not corruption).

## Why this matters
Before: "lossless" was an *observation* on GIAB-shaped data. After: it's a *guarantee* backed by an
adversarial battery — every non-conforming input either round-trips byte-identically or is refused with a
clear message. No silent data loss, no crashes. This is the difference between a prototype and a codec you
can trust. The critical find was the Phred-42 clamp: it silently discarded legal high-quality values —
exactly the kind of bug that destroys trust in a "lossless" tool, and it would never have surfaced on
binned GIAB data. Repro: edge_gen.py generates the battery; run compress→decompress→cmp per file.
