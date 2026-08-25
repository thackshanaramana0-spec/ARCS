# CLAUDE.md — ARCS Benchmark Project

This file governs all AI-assisted work in this repository. Read it fully before touching any script, dataset, or result.

---

## Project Identity

**ARCS** = Assemble → Retain → Compress → Serve  
**Paper title:** "ARCS: a retained assembly for unified lossless FASTQ compression and reference-free genomic analysis"  
**Target journal:** Nature Methods  
**Server:** Ubuntu 24.04, SDC3 Chennai, 12 vCPU, ~90 GB RAM, 250 GB SSD  
**Repo path on server:** `/root/arcs-clean`

---

## Absolute Rules — Never Break These

1. **Do not change the 10 benchmark accessions.** They are locked in `benchmark/DATASET_LOCKED.md`. Do not substitute, re-verify, or re-derive them. If a download fails, report which accession failed and stop.
2. **ARCS compression uses zero flags.** Always: `arcs compress INPUT.fq OUTPUT.arcs` — nothing else. No `--fast`, no `ARCS_PAR_SHARDS`, no `ARCS_CHUNK_THREADS`. Auto-chunk for >2 GB inputs is expected default behaviour, not a problem.
3. **Never run two timed benchmark jobs concurrently.** One tool at a time. Concurrent jobs contaminate timing and RAM measurements.
4. **Never silently alter methodology.** If something fails, stop, print the error, and wait for instructions.
5. **Claim 1 must be fully lossless before Claim 2 starts.** If any `lossless=LOSSY` appears in Claim 1 output, halt immediately.
6. **Projected numbers are sanity checks only.** Fresh server measurements are authoritative. Do not reject a server result because it differs from a projection.
7. **ARCS must pass 7/7 ctests before any benchmark phase runs.**

---

## The 10 Locked Datasets

Source of truth: `benchmark/DATASET_LOCKED.md`. Reproduced here for fast reference.

| # | Accession | Organism | Kingdom | _1 size (uncompressed) | Auto-chunk |
|---|-----------|----------|---------|------------------------|-----------|
| 1 | SRR2584863 | E. coli B REL606 | Bacteria | ~576 MB | No |
| 2 | ERR552797 | M. tuberculosis H37Rv | Bacteria | ~217 MB | No |
| 3 | SRR554369 | P. aeruginosa PAO1 | Bacteria | ~334 MB | No |
| 4 | ERR5181310 | SARS-CoV-2 | Virus | ~30 MB | No |
| 5 | GIAB HG002 chr20 | H. sapiens | Animalia | ~37 MB | No |
| 6 | SRR16357346 | C. elegans N2 (CeNDR) | Animalia | ~1.42 GB | No |
| 7 | DRR976266 | S. cerevisiae | Fungi | ~1.67 GB | No |
| 8 | SRR1945765 | Arabidopsis thaliana | Plantae | ~1.95 GB | No |
| 9 | SRR36741279 | Leishmania major | Protista | ~1.70 GB | No |
| 10 | SRR37283774 | P. falciparum | Protista | ~669 MB | No |

**All 10 files are under 2 GB. Auto-chunk does NOT fire on any dataset.**

**Banned accessions (never use):** SRR390728, SRR988075, SRR065390, SRR327342, SRR1663585, SRR870667, SRR1296601, ERR015526, SRR1294122, ERR174310

**Disk budget:** ~98 GB peak on 250 GB SSD — safe.

**Cross-validation:** Server SPRING bpb for SRR554369 must be ~0.2416 ±2%, for SRR870667 must be ~1.2621 ±2% (PgRC2 2025 gold values). If they deviate more, something is wrong with the benchmark setup.

---

## 3 Claims and 6 Tables

### Claim 1 — COMPACT (T1 + T2)
- **T1:** Archive size (bytes) — ARCS vs SPRING vs Genozip, all 10 datasets
- **T2:** Compress time (s) + decompress time (s) + peak RAM (VmHWM KB), all tools, all datasets
- ARCS expected to win ratio on all 10; speed slower than SPRING/Genozip (assembly cost is expected)

### Claim 2 — FAITHFUL (T3 + T4 + T5)
- **T3:** Het-SNV F1 — ARCS vs DiscoSNP++ vs Kmer2SNP on HG002/HG003/HG004/HG005 chr20
- **T4:** Coverage sweep — ARCS F1 at 10×/15×/20×/30× (HG002)
- **T5:** Het-indel F1 — ARCS vs DiscoSNP++ on HG002 chr20
- Expected: ARCS avg F1 ≈ 0.936 > DiscoSNP++ 0.918 > Kmer2SNP 0.532

### Claim 3 — ADDRESSABLE (T6)
- **T6:** `arcs export` vs SPAdes, `arcs coverage` vs BWA+mosdepth, `arcs query` (unique)
- Expected: export ~50-200× speedup vs SPAdes, coverage ~2-5× speedup vs BWA+mosdepth

---

## Phase Structure

Run phases strictly in order. Do not proceed to the next phase if the current one has failures.

### Phase 0 — Build and test (5-10 min)

```bash
cd /root/arcs-clean
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc)
ctest --test-dir build
```

**Pass condition:** `7/7 tests passed`  
**Stop if:** any test fails — do not proceed

```bash
# Tool check
for t in prefetch fasterq-dump aws spring genozip; do
    command -v "$t" && echo "$t OK" || echo "$t MISSING"
done
```

**Stop if:** any tool is MISSING — install before continuing

### Phase 1 — Download (30-90 min)

```bash
bash benchmark/download.sh /data/fastq 2>&1 | tee /data/download.log
tail -30 /data/download.log
```

**Pass condition:** Last line contains `DOWNLOAD COMPLETE`  
**Stop if:** any line contains `FAIL` — paste the failing phase output for diagnosis  
**Note:** Large downloads — C. elegans (~3.2 GB SRA), T. cacao (~7.7 GB SRA) — take longest

### Phase 2 — Claim 1: Compression benchmark (2-4 hours)

```bash
bash benchmark/benchmark.sh claim1 /data/fastq /root/arcs-clean/build/arcs ./results 2>&1 | tee ./results/claim1.log
```

After completion, run both checks:
```bash
grep -h "TOOL=ARCS" ./results/claim1/*.log
grep "LOSSY" ./results/claim1/*.log || echo "ALL LOSSLESS"
```

**Pass condition:** Every log line ends with `lossless=LOSSLESS`  
**Stop if:** any `LOSSY` — this is a bug, do not proceed to Claim 2  
**Cross-validate:** SPRING bpb for SRR554369 and SRR870667 must match published values ±2%

### Phase 3 — Claim 2: Variant calling (2-3 hours)

```bash
export CHR20_FA=~/refs/chr20.fa
export CHR20_SDF=~/refs/chr20.sdf
export GIAB_TRUTH=~/giab_truth
export DISCO_DIR=~/DiscoSnp
export CONDA_ENV=kmer2snp_r

bash benchmark/benchmark.sh claim2 /data/fastq /root/arcs-clean/build/arcs ./results 2>&1 | tee ./results/claim2.log
cat ./results/claim2/t3_snv_f1.csv
```

**Pass condition:** T3 CSV has F1 values for all 4 GIAB individuals  
**Expected range:** ARCS F1 0.90–0.95, DiscoSNP++ F1 0.88–0.94, Kmer2SNP F1 0.48–0.56  
**Stop if:** ARCS F1 < 0.85 — check VCF liftover and truth region

### Phase 4 — Claim 3: Archive analysis (1-2 hours)

```bash
bash benchmark/benchmark.sh claim3 /data/fastq /root/arcs-clean/build/arcs ./results 2>&1 | tee ./results/claim3.log
cat ./results/claim3/t6_results.csv
```

**Pass condition:** T6 CSV shows export speedup ≥ 40× for E. coli vs SPAdes  
**Stop if:** `arcs export` output is empty FASTA or `arcs coverage` TSV has zero rows

---

## Critical Command Syntax

```bash
# ARCS — zero flags always
arcs compress reads.fq out.arcs                          # compress
arcs decompress out.arcs decoded.fq                      # decompress
arcs compress --call reads.fq out.arcs out.vcf           # compress + call variants
arcs export out.arcs contigs.fa                          # export pseudogenome
arcs coverage out.arcs coverage.tsv                      # per-contig coverage
arcs query out.arcs 0-100000 region.fq                  # NOTE: hyphen, not colon

# SPRING — always needs -g for FASTQ output on decompress
spring -c -i in.fq -o out.spring -t $(nproc) -g         # compress
spring -d -i out.spring -o decoded.fq -t $(nproc) -g    # decompress — -g required

# Genozip
genozip --force -o out.genozip in.fq
genounzip --force -o out.fq in.genozip

# DiscoSNP++ — -G flag is mandatory
run_discoSnp++.sh ... -G chr20.fa                        # -G required, no exceptions

# Lossless check — order-free, CRLF-safe
paste - - - - < orig.fq | tr -d '\r' | sort > /tmp/a
paste - - - - < dec.fq  | tr -d '\r' | sort > /tmp/b
cmp -s /tmp/a /tmp/b && echo LOSSLESS || echo LOSSY
```

---

## What to Paste Back for Debugging

| Situation | Paste this |
|-----------|-----------|
| Phase 0 test failure | Full `ctest` output |
| Download failure | Lines around `FAIL` in download.log |
| LOSSY result in Claim 1 | The full `__ARCS.log` for the failing dataset |
| ARCS F1 < 0.85 in Claim 2 | Last 50 lines of claim2.log + t3_snv_f1.csv |
| SPAdes failure in Claim 3 | Last 30 lines of spades.log |
| Any unexpected output | The exact phase line that failed + 10 lines above it |

---

## Correct Binary Path

```bash
# CORRECT
/root/arcs-clean/build/arcs

# WRONG — old binary, do not use
/usr/local/bin/arcs
```

Always build from source on the server. Never use a system-installed arcs binary.
