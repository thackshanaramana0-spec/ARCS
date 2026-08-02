# scripts/

Evaluation scripts used to validate ARCS compression and variant calling results.
These are **not part of the codec** — the compressor is pure C++ and runs without them.
They require external test data (GIAB / NA12878 / E. coli, not stored in this repo)
placed alongside a built `arcs` binary.

## Variant calling evaluation

| Script | Purpose |
|---|---|
| `eval_caller.py` | Precision/recall evaluation of `arcs call` output vs truth VCF |
| `sim_indel.py` | Simulate reads with injected indels + ground truth VCF |
| `sim_indel_bench.py` | Run indel benchmark across multiple seeds |
| `sim_polyploid.py` | Simulate polyploid reads (triploid/tetraploid) with injected variants |
| `run_giab_indel.sh` | Run indel benchmark on real GIAB HG002 chr20 region |
| `run_indel_bench.sh` | Full indel benchmark sweep |
| `run_polyploid_bench.sh` | Full polyploid benchmark sweep |

## Compression evaluation

| Script | Purpose |
|---|---|
| `sim_ecoli.py` | Simulate E. coli reads for compression testing |

## Usage

All scripts expect:
1. A built `arcs` binary in PATH or the same directory
2. Test FASTQ data downloaded separately (GIAB, NA12878, E. coli)

See individual `docs/*.md` result files for exact data-slicing commands and download instructions.

## Data sources

- **GIAB HG002 / HG001 (NA12878)**: [https://www.nist.gov/programs-projects/genome-bottle](https://www.nist.gov/programs-projects/genome-bottle)
- **E. coli K-12**: NCBI accession NC_000913.3
