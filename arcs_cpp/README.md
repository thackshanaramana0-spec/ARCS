# ARCS-CPP

[![CI](https://github.com/thackshanaramana/arcs-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/thackshanaramana/arcs-cpp/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Paper:** Thackshanaramana B. *ARCS-CPP: reference-free lossless FASTQ compression via LSH-accelerated MST ordering and context-rANS quality coding.* PLOS Computational Biology, 2026 (submitted). Manuscript: [`paper/PAPER_ARCS_CPP.tex`](paper/PAPER_ARCS_CPP.tex)

**ARCS-CPP** is a reference-free lossless FASTQ compressor written in C++17. It achieves state-of-the-art compression on both whole-genome sequencing (WGS) and amplicon datasets while guaranteeing byte-exact round-trip fidelity of all reads, quality scores, and read names.

## Key results (100 k read subsets)

| Dataset | ARCS-CPP | SPRING v1.1.1 | Improvement | Speed |
|---------|----------|---------------|-------------|-------|
| DS-1 E. coli 150 bp (SRR2584863) | 6.37 MB | 6.50 MB | **+2.0%** | 1.35× faster |
| DS-2 Human 126 bp (SRR2962669) | 5.04 MB | 5.29 MB | **+4.7%** | 2.3× faster |
| DS-4 M. tuberculosis 50 bp (ERR11820348) | 2.98 MB | 3.27 MB | **+8.7%** | 1.6× faster |
| DS-5 SARS-CoV-2 amplicon 220 bp (ERR5181310) | 1.22 MB | 1.55 MB† | **+20.8%** | only tool that works |
| DS-7 Human 30× WGS (SRR1238539) | 64.98 MB | 67.56 MB | **+3.8%** | — |

† SPRING v1.1.1 crashes in lossless mode on 220 bp amplicon data; published value used for comparison.

Full verified results are in [`experiments/VERIFIED_RESULTS_V7star_FINAL.md`](experiments/VERIFIED_RESULTS_V7star_FINAL.md).

## How it works

1. **AKC pilot scoring** — 5,000 reads are sampled (configurable via `--pilot`) and compressed with zlib-9 to compute a compressibility score. Scores below 0.15 indicate amplicon data.
2. **Amplicon path** (score < 0.15) — PCR deduplication, cluster-order sort, flat quality string + permutation array encoding, LZMA-9 quality codec.
3. **WGS MST path** (score ≥ 0.15) — LSH-accelerated k-NN graph → Kruskal MST → DFS read ordering → delta-encoded sequences → context-rANS quality codec (860 contexts: 43 Phred × 2 is_dev × 10 pos bins).
4. **ARCS v1 container** — big-endian blob table with named blobs: `GENOME/MST_TREE`, `MST_DELTAS`, `QUALITY_DATA`, `QUALITY_PERM`, `NAMES`, counts.

See [`docs/algorithm.md`](docs/algorithm.md) for full technical details.

## Requirements

| Dependency | Version |
|------------|---------|
| C++ compiler | GCC ≥ 11, Clang ≥ 14, or MSVC 2022 |
| CMake | ≥ 3.16 |
| zlib | ≥ 1.2 |
| liblzma (xz-utils) | ≥ 5.2 |
| OpenMP (optional) | any |

## Building

### Linux / macOS

```bash
sudo apt-get install -y zlib1g-dev liblzma-dev   # Debian/Ubuntu
# or: sudo dnf install -y zlib-devel xz-devel     # Fedora/RHEL

git clone https://github.com/thackshanaramana/arcs-cpp.git
cd arcs-cpp
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Binary: `build/arcs`

### Windows (MinGW-w64 via MSYS2)

```bash
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake \
          mingw-w64-x86_64-zlib mingw-w64-x86_64-xz

mkdir build && cd build
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
mingw32-make -j4
```

Binary: `build/arcs.exe`

## Usage

```bash
# Compress (single-end)
arcs compress input.fastq.gz output.arcs

# Compress (paired-end)
arcs compress_pe R1.fastq.gz R2.fastq.gz output.arcs

# Decompress
arcs decompress output.arcs restored.fastq

# Inspect archive metadata
arcs info output.arcs

# Check sequencing regime before compressing
arcs akc input.fastq.gz
```

Both gzip-compressed (`.fastq.gz`) and plain (`.fastq`) input are accepted. AKC routing (WGS vs amplicon) is fully automatic.

**Key options** (placed before subcommand):

| Option | Default | Description |
|--------|---------|-------------|
| `--lossy` | off | Novaseq 4-bin quality binning (2,12,23,37) |
| `--lossy8` | off | Uniform 8-bin quality quantisation |
| `--lossy2` | off | Binary 2-bin (maximum quality compression) |
| `--knn <int>` | 40 | k-NN neighbours in MST graph |
| `--pilot <int>` | 5000 | Reads used for AKC pilot scoring |
| `--no-names` | off | Omit read names from archive |

Full CLI reference: [`docs/cli_reference.md`](docs/cli_reference.md)

## Reproducing the benchmarks

Download datasets:
```bash
bash scripts/download_datasets.sh datasets/
```

Run all benchmarks:
```bash
bash scripts/run_benchmark.sh ./build/arcs datasets/
```

Exact commands, SRA accession numbers, and SHA-256 checksums used in the paper are in [`paper/SUPPLEMENTARY_COMMANDS.md`](paper/SUPPLEMENTARY_COMMANDS.md).

## Repository layout

```
src/              C++17 source files (28 files)
third_party/      Header-only dependencies (rans_byte.h)
tests/            Unit and round-trip validation tests
datasets/         Dataset accession numbers and MD5 checksums
experiments/      Benchmark results, competitor comparison, ablation study
figures/          Paper figures (Fig 1–3)
scripts/          Dataset download and benchmark automation scripts
docs/             Technical documentation (algorithm, container format)
paper/            Manuscript (PAPER_ARCS_CPP.tex), supplementary commands
.github/          CI workflow (Ubuntu + Windows build + test)
```

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use ARCS-CPP in your research, please cite:

> Thackshanaramana B. ARCS-CPP: reference-free lossless FASTQ compression via
> LSH-accelerated MST ordering and context-rANS quality coding.
> *PLOS Computational Biology*, 2026. (submitted)

```bibtex
@article{thackshanaramana2026arcs,
  author  = {Thackshanaramana, B.},
  title   = {{ARCS-CPP}: reference-free lossless {FASTQ} compression via
             {LSH}-accelerated {MST} ordering and context-{rANS} quality coding},
  journal = {PLOS Computational Biology},
  year    = {2026},
  note    = {submitted},
  url     = {https://github.com/thackshanaramana/arcs-cpp}
}
```

## Author

**Thackshanaramana B.**  
Student Researcher, Department of Computational Intelligence  
SRM Institute of Science and Technology, Chennai, India  
ORCID: [0009-0005-9453-316X](https://orcid.org/0009-0005-9453-316X)  
Email: tb2138@srmist.edu.in
