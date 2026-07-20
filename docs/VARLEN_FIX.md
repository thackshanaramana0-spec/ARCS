# Variable-length read support (bug fix)

Date: 2026-07-17. Fixed the pre-existing bug where chain-pg assumed a FIXED read length
(L = reads[0].seq.size()), so real trimmed FASTQ (variable-length reads) crashed the decoder
("chain-pg: pg position out of bounds", decoder.cpp:853) — losslessness was broken on real data.

## What was wrong
chain-pg used one global L everywhere: placement, consensus polish, merge RC remaps, the aux stream,
the decoder's pg extraction, and the quality coders. Real 2×250 Illumina reads after adapter/quality
trimming have variable length (35..250) → position math + per-read loops overran → decode threw.
Fixed-length FASTQ (DS7, E.coli) was always fine (7/7 tests) — the bug just was never exercised.

## The fix (per-read length end to end)
- `ChainEncodeResult.pg_readlen` — new per-read length vector.
- `record_mapped` / `record_append` — use the read's ACTUAL length (target.size()), store pg_readlen.
- `build_multicontig_pg` placement — shadow `L = orig.size()` per read; seed offsets bounds-guarded
  (`off+K>L`); short reads (<K) fall through to a new contig (no data loss).
- Merge RC remaps (extension + containment) — use each read's length in the flip transform.
- Consensus polish + pileup dump — per-read length.
- AUX serialization — new `col_readlen` stream (uint16 LE/read); header grew 32→36 bytes (9 fields).
- Decoder — `readlen_of(i)` accessor; per-read Li for pg extraction / mismatch / N / RC / dev_sets /
  surprise; passes per-read lengths (`rlens`) to the quality coder.
- Quality — `qual_cm_encode` loops each read's actual quality length; `qual_cm_decode` takes `rlens`
  and produces per-read-length output. Position-binning context still uses global L (identical on
  both sides → lossless). Variable-length data FORCES the CM coder (mode 0x07); the fixed-L static
  rANS path is skipped for var-length (kept for fixed-length via keep-smaller gate).

## Validation (all lossless)
- Synthetic variable-length (80–150 bp, 3000 reads): seq+qual roundtrip LOSSLESS.
- **Real GIAB HG002 2×250 (chr20, lengths 35–250, real trimmed): seq+qual roundtrip LOSSLESS.**
- Fixed-length regression: DS7-100k total_seq 353,065→353,049 (≈unchanged; readlen stream is ~free,
  it LZMA-compresses to nothing when all lengths are equal); 7/7 ctests pass.

## Notes
- Backward compatible: `readlen_of` falls back to global L if the readlen stream is absent (old
  archives). New archives always include it.
- The compressor now handles real trimmed FASTQ — a required production capability.

## Read-name-after-space fix (2026-07-17)
The FASTQ reader (fastq_io.cpp) used to TRIM each name at the first space ("keep only identifier"),
discarding the comment (barcode / mate info / `1:N:0:...`) irrecoverably. Removed the trim → the FULL
header line is kept. It flows through unchanged: the name encoder stores names newline-delimited
(spaces are safe; only newline is forbidden and can't appear in a header) and the decoder splits on
newline. Fast-path PREFIX.N encoding simply doesn't match names-with-comments → LZMA fallback stores
them fully. VALIDATED: synthetic names with comment ("@INSTR:… 1:N:0:ACGTACGT") roundtrip
BYTE-IDENTICAL; real GIAB HG002 full-record (name+seq+qual) LOSSLESS; 7/7 ctests pass. Read names are
now fully lossless.
