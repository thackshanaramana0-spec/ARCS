# Changelog

## v1.0.0 — 2026-06-23

Initial public release accompanying submission to *PLOS Computational Biology*.

### Compression

- Reference-free lossless FASTQ compression (single-end and paired-end)
- AKC pilot scoring for automatic WGS / amplicon regime detection
- WGS path: LSH-accelerated k-NN graph → Kruskal MST → DFS ordering → context-rANS quality codec (200 contexts)
- Amplicon path: PCR dedup → cluster-order sort → LZMA-9 quality codec
- Three lossy quality modes: `--lossy` (Novaseq 4-bin), `--lossy8` (8-bin), `--lossy2` (2-bin)

### Performance

- 11 speed optimisations in `src/mst_encoder.cpp`; DS-1 end-to-end: 38.8s → 5.1s (7.6×)
- SIMD Hamming distance (SSE2), OpenMP parallel k-NN, BFS level-parallel delta encoding
- Deterministic output: identical input always produces bit-identical archives

### Results (100k read subsets vs SPRING v1.1.1)

| Dataset | ARCS-CPP | Improvement | Speed |
|---------|----------|-------------|-------|
| DS-1 E. coli 150 bp | 6.37 MB | +3.8% smaller | 1.4× faster |
| DS-2 Human 126 bp | 5.04 MB | +4.8% smaller | 3.2× faster |
| DS-4 M. tuberculosis 50 bp | 2.95 MB | +9.7% smaller | 2.3× faster |
| DS-5 SARS-CoV-2 amplicon 220 bp | 1.22 MB | +20.8% smaller | only tool that works |

### CLI subcommands

`compress`, `decompress`, `compress_pe`, `info`, `akc`

### Platforms tested

- Windows 11 (MinGW-w64 GCC 13.2.0, CMake Release)
- Linux Ubuntu 24.04 (GCC 13.3.0)
