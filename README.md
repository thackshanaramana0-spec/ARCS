# ARCS — Lossless FASTQ Compression + Reference-Free Variant Calling

ARCS is a reference-free, lossless compressor for Illumina FASTQ sequencing data.
It builds a **self-assembled pseudogenome** from the reads via a greedy k-NN chain walk,
compresses sequence against it with an embedded context-mixing DNA coder, and compresses
quality values with an adaptive range coder conditioned on local sequence context.
As a **zero-cost byproduct** of the compression assembly it also performs reference-free
heterozygous SNV and indel calling — no external aligner or k-mer tool required.

Everything is written in C++17 (no runtime dependencies beyond zlib, liblzma, pthreads)
and compiles to a single self-contained binary — the same class as SPRING, PgRC2, Genozip,
and fqzcomp.

---

## Results

### Compression (GIAB HG002, chr20 region, 39 MB, 113,987 reads, native ext4, byte-lossless-verified)

| Tool | Archive | Compress | Decompress | RAM | Lossless |
|---|---|---|---|---|---|
| **ARCS** | **4.31 MB** | **4.4 s** | **1.9 s** | 366 MB | **yes** |
| Genozip 15 | 5.03 MB | 4.8 s | 0.45 s | 382 MB | yes |
| SPRING | 5.60 MB | 5.0 s | 4.3 s | 319 MB | yes |
| PgRC2 | 0.245 MB† | 4.5 s | 0.07 s | 88 MB | no (lossy quality) |
| fqzcomp v4.6 | 7.26 MB | 0.6 s | 0.9 s | 60 MB | no (quality quirk) |

†PgRC2 default is lossy-quality — not a comparable lossless number.

**ARCS produces the smallest fully-lossless archive of any tested tool: −14% vs Genozip, −23% vs SPRING.**
Wins 8/8 public datasets. Byte-exact lossless including CRLF line endings (SPRING silently corrupts CRLF files).
See [docs/RESULTS.md](docs/RESULTS.md) for full numbers.

### Reference-free variant calling (4 GIAB individuals, 3 ancestries, rtg vcfeval gold standard)

| Tool | Het-SNV F1 | Notes |
|---|---|---|
| **ARCS** | **0.936** | zero-cost byproduct of compression |
| DiscoSNP++ | 0.918 | dedicated reference-free caller |
| Kmer2SNP | 0.532 | dedicated reference-free caller |

`arcs call reads.fq out.vcf` — one command, no aligner, no reference, no extra runtime.
Also supports het-indels (F1 ≈ 0.51 on real GIAB, ties DiscoSNP++) and polyploid calling
(triploid synthetic F1 = 0.994). See [docs/REFFREE_COMPARISON.md](docs/REFFREE_COMPARISON.md).

---

## Build

```bash
# Linux / macOS
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc)
ctest --test-dir build          # 7/7 lossless roundtrip tests
```

**Dependencies:** C++17 compiler, CMake ≥ 3.16, zlib, liblzma.

```bash
# Ubuntu / Debian
sudo apt-get install cmake build-essential zlib1g-dev liblzma-dev

# macOS (Homebrew)
brew install cmake xz
```

---

## Usage

```bash
# Compress (chain-pg pseudogenome mode — default)
./build/arcs compress reads.fastq out.arcs

# Decompress (fully lossless, byte-exact)
./build/arcs decompress out.arcs reads_restored.fastq

# Reference-free variant calling
./build/arcs call reads.fastq variants.vcf
```

### Environment knobs (all optional — defaults are tuned)

| Variable | Effect |
|---|---|
| `ARCS_PG_BLOCKS=1` | Single-stream sequence codec for absolute best ratio (slightly slower decode) |
| `ARCS_QUAL_BLOCKS=N` | Override quality coder thread/block count |
| `ARCS_NAMES_BLOCKS=N` | Override name coder block count |
| `ARCS_PG_FAST_DECODE=1` | LZMA-ASCII sequence codec: 2.8× faster decompress, +0.23% archive (GIAB) |
| `ARCS_MERGE_OFF=1` | Disable contig merge step (default on; merge collapses fragmented contigs for better ratio) |
| `ARCS_ORDER_FREE=1` | Opt into multiset mode (smaller archive, loses read order) |
| `ARCS_ENC_TIMING=1` | Per-phase encode timing |
| `ARCS_DEC_TIMING=1` | Per-phase decode timing |

---

## Repository layout

```
src/          C++17 source — encoder, decoder, chain assembler, quality coder, variant caller
tests/        7 CTest roundtrip and unit tests
test_data/    50-read sample FASTQs for 10 organisms (36 KB total) — run immediately after build
third_party/  rans_byte.h (rANS entropy coder, public domain) + libbsc (BWT+QLFC, MIT)
docs/         Benchmark results and algorithm notes
scripts/      Evaluation scripts (variant calling, indel bench, polyploid sim)
benchmark/    Linux benchmark runner
.github/      CI: build + ctest on Ubuntu 22.04 and macOS 13
```

---

## Quick test

Sample FASTQs (50 reads, 10 organisms) are included in [`test_data/`](test_data/):

```bash
./build/arcs compress test_data/ecoli_sample.gz out.arcs
./build/arcs decompress out.arcs restored.fastq
diff <(zcat test_data/ecoli_sample.gz) restored.fastq && echo "LOSSLESS OK"
```

---

## Reproducing benchmarks

Full benchmark data is not stored in this repository. Download separately:

- **GIAB HG002 / HG001 (NA12878)** — [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
- **E. coli K-12** — NCBI accession NC_000913.3

See [docs/RESULTS.md](docs/RESULTS.md) for full multi-dataset tables and dataset descriptions.

---

## Design

- **Multi-contig self-assembled pseudogenome** — greedy k-NN chain walk builds a pseudogenome
  from reads without a reference; sequence compresses ~10× smaller than k-mer models.
- **Adaptive quality coder** — conditioned on local sequence context + run-length; sits at
  the conditional-entropy floor of binned Illumina quality.
- **Tokenized read-name model** — splits Illumina names into template + binary-packed X:Y coordinates.
- **Block-parallel codecs** — sequence, quality, and name coders auto-scale with core count.
- **Byte-lossless on adversarial input** — empty reads, sub-seed reads, lowercase/IUPAC bases,
  full Phred 0–93, CRLF line endings, control-char names all round-trip correctly.

---

## Citation

If you use ARCS in your research, please cite:

> Thackshanaramana B et al. (2026). *ARCS: reference-free lossless FASTQ compression
> with integrated heterozygous variant calling.* (manuscript in preparation)

---

## License

MIT — see [LICENSE](LICENSE).
