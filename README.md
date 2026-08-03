# ARCS: Unified Lossless FASTQ Compression and Reference-Free Variant Calling

[![CI](https://github.com/thackshanaramana0-spec/ARCS/actions/workflows/ci.yml/badge.svg)](https://github.com/thackshanaramana0-spec/ARCS/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/badge/release-v2.2.0-blue.svg)](https://github.com/thackshanaramana0-spec/ARCS/releases/tag/v2.2.0)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)](#build)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Lossless](https://img.shields.io/badge/lossless-byte--exact-brightgreen.svg)](#results)

> Compress Illumina short reads to the smallest lossless archive and call heterozygous SNVs and indels — single pass, single binary, no reference genome.

---

## The Problem

Whole-genome sequencing now generates petabytes of FASTQ data annually. Existing lossless
compressors (SPRING, Genozip, fqzcomp) treat compression as a storage-only problem and
discard the internal assembly structures they build. Separately, reference-free variant
callers (DiscoSNP++, Kmer2SNP) must reprocess the same reads from scratch, doubling
computation time and memory. No existing tool simultaneously achieves the smallest lossless
archive and performs reference-free variant calling.

## Proposed Solution

ARCS unifies lossless compression and heterozygous variant calling into a single pass over
the reads. The self-assembled pseudogenome built during compression is itself a
haplotype-aware graph. Bubbles in that graph are heterozygous SNV and indel candidates,
making variant calling a zero-cost byproduct of the best-ratio compression algorithm, with
no external aligner, no reference genome, and no additional runtime.

**Figure 1. ARCS two-phase encoding and variant calling pipeline.**

![ARCS Architecture](docs/figures/fig2.png)

*Phase 1 (De Novo Assembly): k-mer graph construction, contig building, and read placement.
Phase 2 (Parallel Encoding Streams): sequence through pseudogenome construction and FCM
arithmetic coder; quality through block-parallel context model; names through Illumina
tokenizer and LZMA-9. All three streams merge into a single lossless `.arcs` binary archive.
The placement index simultaneously feeds a pileup and bubble detection step, producing
reference-free VCF output at zero additional cost.*

## Key Features

- **Smallest lossless archive** in the field: beats Genozip by 14% and SPRING by 23% on GIAB
- **Reference-free variant calling** as a zero-cost byproduct of compression (het-SNV F1 = 0.936)
- **Byte-exact lossless** on all input types including CRLF line endings, full Phred 0-93, lowercase/IUPAC bases
- **Single binary, no runtime dependencies** beyond zlib, liblzma, and pthreads
- **Fast-decode mode** (`--fast-decode`): 8-30x faster decompress at only 2% archive overhead

---

## Contents

- [Results](#results)
- [Quick test](#quick-test)
- [Build](#build)
- [Usage](#usage)
- [Design](#design)
- [Reproducing benchmarks](#reproducing-benchmarks)
- [Citation](#citation)
- [Author](#author)
- [License](#license)

---

## Results

### Compression

**Table 1.** Lossless compression benchmark on GIAB HG002 chr20 (39 MB, 113,987 reads,
150 bp Illumina, native Linux ext4, byte-lossless roundtrip verified). Best-of-2
wall-clock times.

| Tool | Archive size | Compress | Decompress | Peak RAM | Byte-lossless |
|---|---|---|---|---|---|
| **ARCS** | **4.31 MB** | **4.4 s** | **1.9 s** | 366 MB | **yes** |
| Genozip 15 | 5.03 MB | 4.8 s | 0.45 s | 382 MB | yes |
| SPRING | 5.60 MB | 5.0 s | 4.3 s | 319 MB | yes |
| PgRC2 | 0.245 MB† | 4.5 s | 0.07 s | 88 MB | no (lossy quality) |
| fqzcomp v4.6 | 7.26 MB | 0.6 s | 0.9 s | 60 MB | no (quality remapping) |

†PgRC2 default applies lossy quality binning and is not a lossless comparison.

ARCS produces the smallest fully-lossless archive: **14% smaller than Genozip, 23% smaller
than SPRING**. Wins on all 8 public datasets spanning 4 organisms (E. coli, H. sapiens,
M. tuberculosis, SARS-CoV-2) and read lengths from 51 bp to 221 bp. Byte-exact lossless including CRLF line endings, which SPRING silently alters.

With `--fast-decode`: 0.28 s decompress on GIAB (vs 1.9 s default), archive 4.41 MB
(still 12% smaller than Genozip, 21% smaller than SPRING).

See [docs/RESULTS.md](docs/RESULTS.md) for full 8-dataset tables and per-stream breakdowns.

### Reference-free variant calling

**Table 2.** Het-SNV calling benchmark on 4 GIAB individuals (HG002-HG005, 3 ancestries),
5 chr20 windows per individual, evaluated with rtg vcfeval against NIST high-confidence
call sets.

| Tool | Het-SNV F1 | Precision | Recall | Notes |
|---|---|---|---|---|
| **ARCS** | **0.936** | **0.991** | **0.886** | zero-cost byproduct of compression |
| DiscoSNP++ | 0.918 | | | dedicated reference-free caller |
| Kmer2SNP | 0.532 | | | dedicated reference-free caller |

`arcs call reads.fq out.vcf` — one command, no aligner, no reference genome required.

ARCS also supports:
- **Het-indels**: GIAB HG002 chr20 F1 = 0.505; wins vs DiscoSNP++ in homopolymer (0.538 vs 0.417) and tandem-repeat (0.383 vs 0.324) strata
- **Polyploid calling**: triploid synthetic data F1 = 0.994-1.000 (`ARCS_PLOIDY=3`)

See [docs/REFFREE_COMPARISON.md](docs/REFFREE_COMPARISON.md) for complete calling benchmarks.

---

## Quick test

Sample FASTQs (50 reads each, 10 organisms) are included in [`test_data/`](test_data/).
No downloads required:

```bash
./build/arcs compress test_data/ecoli_sample.gz out.arcs
./build/arcs decompress out.arcs restored.fastq
diff <(zcat test_data/ecoli_sample.gz) restored.fastq && echo "LOSSLESS OK"
```

All 10 species:

```bash
for f in test_data/*.gz; do
    ./build/arcs compress "$f" /tmp/t.arcs
    ./build/arcs decompress /tmp/t.arcs /tmp/t.fq
    diff <(zcat "$f") /tmp/t.fq && echo "OK: $f" || echo "FAIL: $f"
done
```

Full lossless regression suite (7 CTests): `ctest --test-dir build`

---

## Build

```bash
git clone https://github.com/thackshanaramana0-spec/ARCS.git
cd ARCS
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc)
ctest --test-dir build
```

**Dependencies:** C++17 compiler, CMake 3.16 or later, zlib, liblzma. libzstd is optional but required for `--fast-decode` mode.

```bash
# Ubuntu / Debian
sudo apt-get install cmake build-essential zlib1g-dev liblzma-dev libzstd-dev

# macOS (Homebrew)
brew install cmake xz
```

**Install system-wide** (puts `arcs` on your PATH):

```bash
sudo cmake --install build
```

**Portable build** that runs on any x86-64 machine (SSE4.1, 2010 onwards):

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DARCS_NATIVE=OFF -DARCS_PORTABLE=ON
```

By default the build targets your own CPU (`-march=native`) for maximum performance.

---

## Usage

```bash
# Compress (best ratio, default)
arcs compress reads.fastq out.arcs

# Compress with fast-decode (8-30x faster decompress, ~2% larger archive)
arcs compress --fast-decode reads.fastq out.arcs

# Decompress (byte-exact; format auto-detected)
arcs decompress out.arcs reads_restored.fastq

# Reference-free variant calling (standalone)
arcs call reads.fastq variants.vcf

# Fused compress and call (single pass, byte-identical archive)
arcs compress --call reads.fastq out.arcs

# Archive inspection
arcs info out.arcs
```

### Environment knobs

All defaults are tuned for best lossless ratio. These are for advanced use only.

| Variable | Effect |
|---|---|
| `ARCS_PG_BLOCKS=1` | Single-stream pg codec for absolute best ratio, slower decode |
| `ARCS_QUAL_BLOCKS=N` | Quality coder thread and block count override |
| `ARCS_NAMES_BLOCKS=N` | Name coder block count override |
| `ARCS_MERGE_OFF=1` | Disable contig merge step (enabled by default) |
| `ARCS_ORDER_FREE=1` | Multiset mode: smaller archive but loses read order |
| `ARCS_PAR_SHARDS=N` | Parallel shard count for assembly (default: min(4, cores)) |
| `ARCS_ENC_TIMING=1` | Per-phase encode timing |
| `ARCS_DEC_TIMING=1` | Per-phase decode timing |

---

## Design

**Figure 2. ARCS detailed algorithm: from FASTQ input to lossless archive and VCF output.**

![ARCS Algorithm](docs/figures/fig1.png)

*(A) FASTQ reads are loaded and pilot-scored. (B) A de novo string graph is assembled
without a reference. (C) Reads are placed against assembled contigs. (D) A compact
placement index stores read-to-contig mappings. (E) The pseudogenome is constructed
from linearised contigs; variant evidence is simultaneously extracted via pileup base
counts. (F) Sequences are entropy-coded with an FCM arithmetic coder; heterozygous
bubbles in the assembly graph are detected as variant candidates. (G) The final output
is a lossless `.arcs` binary archive and a reference-free VCF file, produced in a
single pass over the input.*

### Algorithmic contributions

- **Multi-contig self-assembled pseudogenome.** A greedy k-NN chain walk builds a pseudogenome
  without a reference. Sequence compresses approximately 10x smaller than k-mer models alone.
- **Adaptive quality coder.** Conditioned on local sequence context and run-length; reaches
  the conditional-entropy floor of binned Illumina quality (1.574 bpq, matching fqzcomp).
- **Tokenized read-name model.** Splits Illumina names into template and binary-packed X:Y
  flowcell coordinates; names reach their information-theoretic floor.
- **Block-parallel codecs.** Sequence, quality, and name coders auto-scale with core count.
- **Parallel shard assembly.** Reads split round-robin into N shards, assembled in parallel,
  merged with globally remapped indices; 4x compress speedup on multi-core hardware.
- **Fast-decode mode.** 2-bit packed DNA (zstd-6) and a precomputed 2D quality table deliver
  Genozip-class decode speed while retaining the smallest lossless archive in the field.
- **Adversarial robustness.** Empty reads, sub-seed reads, lowercase and IUPAC bases,
  full Phred 0-93, CRLF line endings, and control-character names all round-trip correctly.

---

## Repository layout

```
src/           C++17 source: encoder, decoder, chain assembler, quality coder, variant caller
tests/         7 CTest roundtrip and unit tests
test_data/     50-read sample FASTQs for 10 organisms (36 KB total)
third_party/   rans_byte.h (rANS, public domain) + libbsc (BWT+QLFC, MIT)
docs/          Benchmark results, algorithm notes, pipeline figures
scripts/       Evaluation scripts: variant calling, indel bench, polyploid sim
benchmark/     Parametric Linux benchmark runner
.github/       CI: build and ctest on Ubuntu 22.04 and macOS 13
```

---

## Reproducing benchmarks

Full benchmark datasets are not stored in this repository. Download separately:

- **GIAB HG002 / HG001 (NA12878):** [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
- **E. coli K-12:** NCBI accession NC_000913.3

See [docs/RESULTS.md](docs/RESULTS.md) for full multi-dataset tables, accession numbers,
and dataset descriptions.

Parametric benchmark script for Linux:

```bash
bash benchmark/run_linux.sh ~/fastq_datasets ./build/arcs ./results
```

---

## Citation

If you use ARCS in your research, please cite:

> Thackshanaramana B (2026). *ARCS: reference-free lossless FASTQ compression with
> integrated heterozygous variant calling.* Manuscript in preparation.

---

## Author

**Thackshanaramana B**  
SRM Institute of Science and Technology  
tb2138@srmist.edu.in · thackshanaramanab@gmail.com  
[ORCID: 0009-0005-9453-316X](https://orcid.org/0009-0005-9453-316X)

---

## License

MIT. See [LICENSE](LICENSE).

`third_party/libbsc/` is distributed under the Apache 2.0 license. See `third_party/libbsc/LICENSE`.  
`third_party/rans_byte.h` is public domain (Fabian Giesen).
