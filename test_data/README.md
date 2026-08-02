# test_data

Small sample FASTQ files (50 reads each, gzip-compressed) for quick correctness testing
across 10 organisms and read-length profiles. Total size: ~36 KB.

| File | Organism | Read length |
|---|---|---|
| `ecoli_sample.gz` | *E. coli* K-12 | 100 bp |
| `giab_hg002_sample.gz` | Human (GIAB HG002) | 150 bp |
| `human100_sample.gz` | Human WGS | 100 bp |
| `human127_sample.gz` | Human WGS | 127 bp |
| `mouse_sample.gz` | *Mus musculus* | 100 bp |
| `celegans_sample.gz` | *C. elegans* | 100 bp |
| `drosophila_sample.gz` | *D. melanogaster* | 100 bp |
| `zebrafish_sample.gz` | *D. rerio* | 100 bp |
| `mtb_sample.gz` | *M. tuberculosis* | 100 bp |
| `sars2_sample.gz` | SARS-CoV-2 | 150 bp |

## Quick test

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc)

./build/arcs compress test_data/ecoli_sample.gz out.arcs
./build/arcs decompress out.arcs restored.fastq
diff <(zcat test_data/ecoli_sample.gz) restored.fastq && echo "LOSSLESS OK"
```

For all 10 samples:

```bash
for f in test_data/*.gz; do
    ./build/arcs compress "$f" /tmp/t.arcs
    ./build/arcs decompress /tmp/t.arcs /tmp/t.fq
    diff <(zcat "$f") /tmp/t.fq && echo "OK: $f" || echo "FAIL: $f"
done
```

Full lossless regression suite: `ctest --test-dir build`.
