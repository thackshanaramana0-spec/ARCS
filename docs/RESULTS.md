# trial20 — verified results (re-run 2026-07-16)

WSL binary `/c/Temp/arcs_test/arcs`. Datasets on disk: `test5k.fastq` (5k), `ds7_100k.fastq` (100k).
No DS1-100k exists on disk — only DS7-100k + the 5k set.

## Sequence track: total_seq = pg_blob + pos_blob + aux_blob

| Assembler stage | 5k (PgRC2 = 18,379) | DS7-100k (PgRC2 = 360,025) |
|---|---|---|
| chain (legacy) | 27,749 | 591,882 |
| dedup | 20,802 | 468,255 |
| **multi-contig (default)** | **15,772  (−14.2%)** | **353,065  (−1.9%)** |

### DS7-100k multi-contig internals
- pg_blob 222,633 B (GeCo3 won over ARCS-DNA 251,975 by 29,343 B)
- pos_blob 51,116 B (zigzag-delta varint)
- aux_blob 79,316 B (RC + mismatches + N; qmm derived, not stored)
- 6,279 contigs, 93,721 reads mapped (93.7%), pg_len 1,498,887 B (9.87× < raw reads)
- full archive: 5,695,933 B

## Losslessness
Multi-contig is **record-lossless** — every name+seq+quality preserved exactly
(verified `paste - - - - | sort | cmp` = MULTISET_LOSSLESS_OK, 100000/100000 records).
It **reorders** reads into pg-position order, so a raw `cmp` on the file "differs" purely
on ordering (no data loss). Same contract as PgRC2 (reorder-allowed), so total_seq is
apples-to-apples. Strict input-order would cost extra permutation bytes.

## ARCS-DNA compressor (standalone)
0.469 bpb on a realistic 9.87M-base pg; 12s / 679MB vs GeCo3 41s / 2.4GB (0.396 bpb).
~18% behind GeCo3 on ratio, 3.4× faster, 3.5× lighter. `compress_pg` keeps the smaller
of ARCS-DNA (codec 0x04) / GeCo3 (0x03) / LZMA (0x01) per file.

## Context
Sequence is ~11% of the FASTQ archive; quality is ~80% (~4.8 MB of the 5.7 MB DS7-100k archive).
Beating PgRC2 on sequence saves ~3% of the total → quality (idea B) is the real lever.
