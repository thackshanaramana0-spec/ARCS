# Evaluation & analysis scripts

These are the Python scripts used to measure and validate ARCS. They are **not part of the codec** —
the compressor is pure C++ and runs without them. They require external test data (GIAB / NA12878 /
E. coli, not stored in this repo) placed alongside a built `arcs` binary.

## Reference-free variant calling (the byproduct evaluation)
- `sim_ecoli.py` — simulate E. coli reads with injected errors/variants + ground truth.
- `eval_caller.py` / `eval_pr2.py` — genotype-likelihood caller precision/recall on the sim.
- `rf_eval.py` / `rf_eval_hg001.py` — reference-free calling vs GIAB truth (HG002 / NA12878).
- `rf_dual.py` / `rf_frozen2.py` — the bubble-cleanliness caller (dual-region tune + frozen held-out test).
- `kmer2snp_eval_argv.py` — score Kmer2SNP output on the same slice for head-to-head.
- `pileup_bwa.py` — build the self-alignment pileup keyed by (contig, pos).

## Compression-headroom analysis (why quality/names are at their floor)
- `flowcell.py` — held-out cross-entropy of quality with/without flowcell-2D context (negative result).
- `namefloor.py` / `nametok.py` — read-name information floor and the tokenized-name model measurement.
- `qent.py` — per-context conditional entropy of quality (run-length lever).

## Tooling
- `edge_gen.py` — generates the 12 adversarial edge-case FASTQ files (lossless stress test).
- `arcs_format.py` — parses and dumps the structure of a `.arcs` archive (header + blob table + streams).

Usage examples and exact data-slicing commands are in the corresponding `docs/*.md` result files.
