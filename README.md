# ARCS — a lossless FASTQ compressor

ARCS is a reference-free, lossless compressor for FASTQ sequencing data. It builds a
**self-assembled pseudogenome** from the reads (multi-contig greedy assembler), compresses the
sequence against it with an embedded context-mixing DNA coder, and compresses quality values with an
adaptive range coder conditioned on the local sequence. As a zero-cost byproduct it also performs
reference-free error correction and SNV calling.

Everything is written in C++ (no runtime dependencies beyond zlib / liblzma / pthread) and compiles to
a single self-contained binary — the same class as SPRING, PgRC2, Genozip, and fqzcomp.

## Results (real GIAB HG002, chr20 region, 39 MB, native ext4, all byte-lossless-verified)

| tool | archive | compress | decompress | comp RAM | lossless |
|---|---|---|---|---|---|
| **ARCS** | **4.28 MB** | **4.4 s** | 1.9 s | 366 MB | **yes** |
| Genozip 15 | 5.03 MB | 4.8 s | 0.45 s | 382 MB | yes |
| SPRING | 5.60 MB | 5.0 s | 4.3 s | 319 MB | yes |
| PgRC2 | 0.245 MB* | 4.5 s | 0.07 s | 88 MB | no (lossy quality) |
| fqzcomp v4.6 | 7.26 MB | 0.6 s | 0.9 s | 60 MB | no (quality quirk) |

*PgRC2's default is lossy-quality — not a comparable lossless number.

**ARCS produces the smallest fully-lossless archive of any tool tested (−15% vs Genozip, −24% vs SPRING),
while matching them on compress speed and RAM.** Thread-for-thread (all at physical-core count) it wins
ratio and compress speed together. See [docs/REAL_BENCHMARK_GIAB.md](docs/REAL_BENCHMARK_GIAB.md).

Reference-free SNV calling (NA12878/HG001, held-out region): ARCS F1 **0.84** vs Kmer2SNP 0.82,
DiscoSNP++ 0.71 — competitive-to-ahead as a free byproduct. See
[docs/REFFREE_NA12878_COMPARISON.md](docs/REFFREE_NA12878_COMPARISON.md).

## Build

```bash
cd arcs_cpp
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
ctest            # 7/7 lossless roundtrip tests
```

Requires a C++17 compiler, CMake ≥ 3.16, zlib, and liblzma.

## Usage

```bash
# compress (default multi-contig pseudogenome mode)
./arcs --chain-pg compress reads.fastq out.arcs

# decompress (fully lossless)
./arcs decompress out.arcs reads.out.fastq
```

Useful environment knobs (all optional; defaults are tuned):
- `ARCS_PG_BLOCKS=1` — single-stream sequence codec for absolute-best ratio (slightly slower decode)
- `ARCS_QUAL_BLOCKS=N` / `ARCS_NAMES_BLOCKS=N` — override thread/block auto-scaling
- `ARCS_USE_GECO3=1` — use the external GeCo3 sequence coder (~0.5% smaller, much slower)
- `ARCS_ENC_TIMING=1` / `ARCS_DEC_TIMING=1` — per-phase timing
- `ARCS_MERGE_CONTIGS=1` — collapse redundant contigs (used by the reference-free calling path)

## Repository layout

```
arcs_cpp/     the codec (C++17, one binary) — src/, tests/, CMakeLists.txt
docs/         results and design notes (the paper trail); docs/notes/ = prior-art surveys
scripts/      Python evaluation / analysis scripts (need external test data — see below)
```

## Reproducing the benchmarks

The benchmarks use public data not stored in this repo:
- **GIAB HG002 / HG001 (NA12878)** — NIST Genome in a Bottle (BAMs, truth VCFs, confident BEDs).
- **E. coli K-12** reference (NCBI NC_000913.3).

The `scripts/` folder contains the evaluation harness (reference-free calling, quality-headroom
analysis, the adversarial edge-case generator, and an archive-format inspector). Each expects the test
files alongside the `arcs` binary; see the individual result docs for the exact slicing commands.

## Design highlights

- **Multi-contig self-assembled pseudogenome** → sequence compresses ~10× smaller than k-mer models.
- **Adaptive quality coder** conditioned on local sequence + run-length; sits at the conditional-entropy
  floor of binned Illumina quality (see [docs/QUALITY_HEADROOM_WALL.md](docs/QUALITY_HEADROOM_WALL.md)).
- **Tokenized read-name model** — splits Illumina names into template + binary-packed coordinates.
- **Block-parallel** sequence, quality, and names codecs that auto-scale with cores.
- **Byte-lossless on adversarial input** (empty reads, sub-seed reads, lowercase/IUPAC bases, full Phred
  0–93, control-char names); see [docs/EDGE_HARDENING.md](docs/EDGE_HARDENING.md).

## License

MIT — see [LICENSE](LICENSE).
