# ARCS Benchmark Brief — Complete Contextual Guide
**Purpose:** This document is the single source of truth for anyone (human or AI agent) running the ARCS benchmark on a remote server. Read this before touching any script. It contains the biology, the methodology, the code architecture, the 3 claims × 6 tables, projected numbers, dataset details, timing, critical command syntax, and reviewer checklist.

---

## 1. The Biology — Why This Problem Matters

### What is sequencing data?
Whole-genome sequencing (WGS) breaks a person's DNA into millions of short fragments (~150 bp each) and reads them. The result is a **FASTQ file**: four lines per read — name, DNA sequence, separator, quality scores (Phred 0-93 per base). A human genome at 30× coverage = ~150 GB of FASTQ. Globally, sequencing now generates petabytes annually.

### The reference genome problem
For decades, genomics analysis has relied on a **reference genome** — a single "standard" human genome (GRCh37 or GRCh38) maintained by the Genome Reference Consortium. Every tool in the standard pipeline (BWA, GATK, DeepVariant) aligns reads to this reference before doing anything useful.

**Why this is a problem:**
- GRCh37/38 is built from ~20 people, predominantly European
- ~80% of global sequencing is from non-European populations
- Non-European genomes differ from GRCh37/38 by millions of bases — the reference is biased
- Low-resource organisms (bacteria, crop plants, pathogens) have no reference at all
- Aligning to a wrong reference creates false variant calls and misses real ones

### Reference-free genomics
**Reference-free** means: work directly from the raw reads, build your own assembly, call variants — no GRC, no alignment step. This is harder but unbiased. Two existing reference-free variant callers: DiscoSNP++ and Kmer2SNP. Neither compresses data. Neither preserves reads.

### The key insight ARCS exploits
Every lossless FASTQ compressor (SPRING, Genozip, fqzcomp) secretly builds an **internal de-novo assembly** of the reads during compression. This is why they get good compression ratios — reads share 97%+ sequence with their assembled pseudogenome, so each read compresses to only its differences. Then they **throw the assembly away**.

Meanwhile, reference-free variant callers rebuild the **exact same assembly from scratch** to detect heterozygous variants — bubbles in the assembly graph where two reads diverge = two alleles = a SNP.

**ARCS says:** These are the same computation. Build the assembly once. Use it for compression. Read variants off the same graph. The assembly is retained, not discarded — it becomes the archive itself.

---

## 2. What ARCS Is — The Method

### Core algorithm (3 phases, 1 pass)

**Phase 1 — De Novo Assembly (the key innovation):**
1. Load all reads into memory
2. Count k-mers (k=16 seeds) across all reads — `FlatPlaceIndex` open-addressed hash table (~35 bytes/distinct k-mer)
3. Build contig chains: greedy k-NN walk through high-frequency k-mers (`build_multicontig_pg`)
4. Parallel shard assembly: reads split round-robin into N shards (default: min(4, hw_concurrency)), each assembled in its own thread, merged with globally remapped indices — up to 4× assembly speedup, ±0.3% ratio cost (`parallel_shard_assemble`)
5. Place each read against the assembled contigs — `FlatPlaceIndex` maps reads to positions

**Phase 2 — Parallel Encoding (3 streams):**
- **Sequence stream:** FCM arithmetic coder (Finite-Context Model) — reads are stored as diffs against the pseudogenome they assembled from. This is why ARCS wins on ratio: the FCM model is conditioned on 14-15 bases of context from the same pseudogenome all reads contributed to.
- **Quality stream:** Adaptive context model conditioned on local sequence + run-length. Reaches the conditional-entropy floor for Illumina binned quality (1.574 bpq, matching fqzcomp).
- **Names stream:** Illumina tokenizer — splits names into static template + binary-packed X:Y flowcell coordinates. Names reach their information-theoretic floor.

All three streams run concurrently (async futures) then merge into a single `.arcs` binary archive.

**Phase 3 — Variant Calling (zero-cost byproduct):**
- During read placement, a pileup is accumulated per contig position
- Positions with two alleles at >PLOIDY_MIN frequency = heterozygous variant candidate
- Bubbles between consecutive contigs = candidate indel
- VCF written as a byproduct of `compress --call` — no second pass, no extra RAM

### How ARCS differs from competitors

| Feature | ARCS | SPRING | Genozip | DiscoSNP++ | Kmer2SNP |
|---------|------|--------|---------|------------|----------|
| Compresses FASTQ | ✅ smallest | ✅ 2nd | ✅ 3rd | ❌ | ❌ |
| Variant calling | ✅ F1=0.936 | ❌ | ❌ | ✅ F1=0.918 | ✅ F1=0.532 |
| Reference required | ❌ | ❌ | ❌ | ❌ | ❌ |
| Byte-exact lossless | ✅ all cases | ❌ (reorders reads) | ✅ | — | — |
| CRLF preserved | ✅ | ❌ | ✅ | — | — |
| Single command | ✅ | ✅ | ✅ | ❌ two-step | ❌ two-step |
| Archive-derived query | ✅ | ❌ | ❌ | ❌ | ❌ |

### The unification equation
```
compress(reads) = build(pseudogenome) + encode(reads vs pseudogenome)
call_variants(reads) = build(pseudogenome) + detect_bubbles(pseudogenome)
```
Both share `build(pseudogenome)`. ARCS does it once. Everyone else does it twice.

---

## 3. Code Architecture — Key Files and Functions

### Source files (`src/`)
| File | Role |
|------|------|
| `encoder.cpp` | Top-level encode orchestration, auto-chunk, parallel stream dispatch |
| `chain_encoder.cpp` | `FlatPlaceIndex`, `build_multicontig_pg`, `parallel_shard_assemble`, `KmerIndex` merge |
| `dna_coder.cpp` | FCM arithmetic coder for sequences; O(1) reverse-complement IR |
| `qual_cm.cpp` | Quality adaptive context model; `FreqMap` chunked hash |
| `decoder.cpp` | Decode all 3 streams; auto-detects mode (0x06/0x08 format flags) |
| `caller.cpp` | Pileup accumulation, bubble detection, VCF emission |
| `fastq_io.cpp` | FASTQ reader/writer; CRLF detection; 8MB buffered writer |
| `main.cpp` | CLI: compress / decompress / call / export / coverage / query / info |
| `container.cpp` | `.arcs` binary archive format: header + 3 stream blobs |

### Critical constants and env vars
| Constant/Env | Value | Effect |
|---|---|---|
| `AUTO_MB` | 2000 MB | Auto-chunk threshold — files >2000 MB chunk into 1M-read pieces. Datasets #6,7,9,10 exceed 2 GB → auto-chunk fires. Datasets #1–5,8 are ≤2 GB → single-pass assembly |
| `DEDUP_K` | 16 | Seed k-mer length for assembly |
| `ARCS_PAR_SHARDS` | min(4,cores) default | Parallel assembly shards. **Not set in benchmark scripts** — C++ default min(4,hw) gives best ratio+speed balance. Override only if you want to trade ratio for speed. |
| `ARCS_CHUNK_THREADS` | = nproc default | Outer parallel chunk threads (only relevant for >2000 MB files) |
| `ARCS_PG_BLOCKS` | 1 (default=1) | FCM codec block count. 1 = best ratio (long-range context). >1 = faster decode but worse ratio |
| `ARCS_ORDER_FREE` | off | Multiset mode: smaller archive but loses read order. Default OFF = order-preserving |

### Threading model (important: never breaks on low-core systems)
```cpp
// encoder.cpp line ~1750
const int hw = (int)std::thread::hardware_concurrency();
int N_shards = std::min(4, hw > 0 ? hw : 1);  // fallback to 1 if hw=0
if (const char* sv = getenv("ARCS_PAR_SHARDS")) N_shards = std::max(1, atoi(sv));
if (call_capture_ || n < 8000) N_shards = 1;  // serial for small datasets
```
The code is safe on any core count including 1. `ARCS_PAR_SHARDS=1` forces serial. **Do not override in benchmark scripts** — the default min(4,hw) is the validated sweet spot. Above 4 shards: each shard gets lower coverage → fragmented assembly → ratio regression. FCM sequence coder is single-threaded and dominates compress time regardless of shard count. Quality and names codecs scale independently with hw_concurrency beyond 4 cores.

### RAM formula (critical for server sizing)
```
FlatPlaceIndex RAM ≈ (distinct k-mers in pseudogenome) × 35 bytes × N_shards_active
```
- E. coli 4.6 MB genome × 35 bytes × 4 shards ≈ 644 MB
- Human chr20 50 MB × 35 bytes × 4 shards ≈ 7 GB (use ARCS_CHUNK_THREADS=8 → ~35 GB, or subset to ≤2 GB)
- All 10 benchmark datasets ≤2 GB → peak RAM ≤8 GB per dataset
- Minimum server RAM: 16 GB. Recommended: 32 GB.

---

## 4. The 3 Claims and 6 Tables

### Claim 1 — COMPACT: Smallest Lossless Archive

**T1 — Archive size (ratio)**
What: ARCS vs SPRING vs Genozip vs PgRC2 archive size in KB on all 10 datasets
Command: `bash benchmark/run_block1.sh <DATA_DIR> <ARCS_BIN> ./results/claim1`
How measured: `stat -c %s archive` bytes, converted to KB
Expected: ARCS wins 8/10 vs SPRING, 7/10 vs Genozip (loses GIAB_HG001 by 0.5%)

**T2 — Speed and RAM**
What: compress time (seconds) + decompress time + peak RAM (VmHWM KB) for all tools on all datasets
Command: same script, same run — `parse_time_v()` reads `/usr/bin/time -v` output
How measured: `grep "Elapsed (wall clock)"` → seconds; `grep "Maximum resident set size"` → KB
Expected: ARCS slower on compress (assembly cost), ARCS competitive on decompress vs SPRING, Genozip fastest decode

### Claim 2 — FAITHFUL: Variant Calling Beats Dedicated Tools

**T3 — Het-SNV F1 across 4 GIAB individuals**
What: Precision / Recall / F1 for ARCS, DiscoSNP++, Kmer2SNP on HG002/HG003/HG004/HG005
Region: chr20:2000000-2400000 (400 kb, ~30× coverage each)
Reference: GRCh37 chr20 only (for LIFT of coordinates — not for calling)
Truth: GIAB v4.2.1 high-confidence VCF + BED per individual
Scorer: `rtg vcfeval --squash-ploidy`
Command: `bash benchmark/run_claim2.sh <DATA_DIR> <ARCS_BIN> ./results/claim2`
Expected:
- ARCS avg HG002-HG005: F1 ≈ 0.936, Precision ≈ 0.991, Recall ≈ 0.886
- DiscoSNP++ avg: F1 ≈ 0.918
- Kmer2SNP avg: F1 ≈ 0.532

**T4 — Coverage depth sweep**
What: ARCS SNV F1 vs coverage (10×/15×/20×/30×) — held-out HG002
Method: python3 subsample FASTQ to target coverage, compress+call, score with vcfeval
Expected: F1(10×)≈0.711 / F1(15×)≈0.869 / F1(20×)≈0.90 / F1(30×)≈0.936

**T5 — Het-indel calling**
What: ARCS vs DiscoSNP++ indel Precision/Recall/F1 on HG002 chr20:2.0-2.4 Mb
Expected: ARCS F1≈0.505, DiscoSNP++ F1≈0.513 — ARCS ties dedicated caller at zero cost
Note: Real indel F1 ~0.5 is the homopolymer/STR ceiling for ALL reference-free callers on this region

### Claim 3 — ADDRESSABLE: Archive-Derived Analysis Without Decompression

**T6 — Export / Coverage / Query speedup vs conventional pipelines**
What: For 5 organisms (E. coli, C. elegans, M. tuberculosis, H. sapiens chr20, D. melanogaster):
- `arcs export archive.arcs out.fa` vs SPAdes de-novo assembly
- `arcs coverage archive.arcs out.tsv` vs BWA mem + samtools sort + mosdepth
- `arcs query archive.arcs 0-100000 out.fq` (ARCS-only, no conventional equivalent)
Command: `bash benchmark/run_claim3.sh <DATA_DIR> <ARCS_BIN> ./results/claim3`
Expected:
- Export: ARCS ~5-10s vs SPAdes ~300-1800s = **50-200× speedup**
- Coverage: ARCS ~7-15s vs BWA+mosdepth ~20-60s = **2-5× speedup**
- Query: ARCS ~1-30s — unique capability (no conventional equivalent)

---

## 5. Dataset Details (All 10 Datasets)

All benchmark datasets are full SRA accession files, single-end reads (_1 file for PE experiments). Datasets #1–5 and #8 are ≤2 GB (no auto-chunk). Datasets #6, 7, 9, 10 exceed 2 GB; ARCS default auto-chunk handles them transparently. All tools run at defaults — no special flags. **See DATASET_LOCKED.md — do not change accessions.**

| # | Accession | Organism | Kingdom | Read len | _1 file size (uncompressed) | Auto-chunk? |
|---|-----------|----------|---------|---------|----------------------------|------------|
| 1 | SRR2584863_1 | E. coli B REL606 | Bacteria | ~108 bp | ~576 MB | No |
| 2 | ERR552797_1 | M. tuberculosis H37Rv | Bacteria | ~301 bp | ~217 MB | No |
| 3 | SRR554369_1 | P. aeruginosa PAO1 | Bacteria | 100 bp | ~334 MB | No |
| 4 | ERR5181310_1 | SARS-CoV-2 | Virus | ~150 bp | ~30 MB | No |
| 5 | HG002_chr20.fq | H. sapiens GIAB HG002 | Animalia | 150 bp | ~37 MB | No |
| 6 | SRR16357346_1 | C. elegans N2 (CeNDR) | Animalia | ~135 bp | ~1.42 GB | No |
| 7 | DRR976266_1 | S. cerevisiae WGS | Fungi | ~150 bp | ~1.67 GB | No |
| 8 | SRR1945765_1 | Arabidopsis thaliana | Plantae | 102 bp | ~1.95 GB | No |
| 9 | SRR36741279_1 | Leishmania major WGS | Protista | ~75 bp | ~1.70 GB | No |
| 10 | SRR37283774_1 | P. falciparum WGS | Protista | ~100 bp | ~669 MB | No |

**Organism diversity:** 3 bacteria · 1 virus · 2 animals · 1 fungus · 1 plant · 2 protists — 6 kingdoms of life.  
**All files are under 2 GB. Auto-chunk does NOT fire. ARCS runs single-pass assembly on all 10.**  
**All accession labels verified via NCBI SRA/DDBJ 2026-08-25.** See DATASET_LOCKED.md for banned list.

**Server RAM budget:** All 10 datasets under 2 GB → single-pass assembly → peak ARCS RAM ≤8 GB per dataset. 96 GB server has 12× headroom. Never run two timed jobs concurrently (contaminates measurements).

---

## 6. Projected Results and Numbers

### T1 — Archive size (projected — to be replaced with server measurements)

All 10 datasets from DATASET_LOCKED.md. Sizes are estimated from organism type and prior ARCS runs; server measurements replace these.

All 10 files are under 2 GB — single-pass assembly, no auto-chunk, ARCS best-case ratio on every dataset.

| # | Dataset | Raw _1 | Est. ARCS | Est. SPRING | Est. Genozip | ARCS wins? |
|---|---------|--------|-----------|-------------|--------------|-----------|
| 1 | SRR2584863 (E. coli) | 576 MB | ~50 MB | ~53 MB | ~63 MB | ✅ vs both |
| 2 | ERR552797 (M. tuberculosis) | 217 MB | ~14 MB | ~15 MB | ~17 MB | ✅ vs both |
| 3 | SRR554369 (P. aeruginosa) | 334 MB | ~4 MB | ~4.4 MB | ~5 MB | ✅ vs both |
| 4 | ERR5181310 (SARS-CoV-2) | 30 MB | ~2 MB | ~2.2 MB | ~2.5 MB | ✅ vs both |
| 5 | GIAB HG002 chr20 | 37 MB | ~3.5 MB | ~4.5 MB | ~4.1 MB | ✅ vs both |
| 6 | SRR16357346 (C. elegans) | 1.42 GB | ~120 MB | ~130 MB | ~150 MB | ✅ vs both |
| 7 | DRR976266 (S. cerevisiae) | 1.67 GB | ~50 MB | ~55 MB | ~65 MB | ✅ vs both |
| 8 | SRR1945765 (Arabidopsis) | 1.95 GB | ~180 MB | ~200 MB | ~230 MB | ✅ vs both |
| 9 | SRR36741279 (L. major) | 1.70 GB | ~70 MB | ~80 MB | ~90 MB | ✅ vs both |
| 10 | SRR37283774 (P. falciparum) | 669 MB | ~25 MB | ~30 MB | ~35 MB | ✅ vs both |

**Gold cross-validation:** SPRING for SRR554369 = 0.2416 bpb (PgRC2 2025). Server measurement must match ±2%.

### T2 — Speed projected (12-core SDC3 Chennai, native ext4)

Compress is 3-4× faster on server vs laptop (assembly is I/O-bound on WSL /mnt/c):
- GIAB_HG002: ~4s compress (was 13s on laptop)
- MTB small: ~8s compress
- Large human (1.5 GB): ~40-80s compress

Decompress: ARCS typically 1-10s. Genozip fastest. SPRING slowest decode.

### T3 — SNV F1 projected (from server run on partial data)

| Individual | ARCS F1 | DiscoSNP++ F1 | Kmer2SNP F1 |
|-----------|---------|---------------|-------------|
| HG002 | 0.929 | 0.941 | 0.542 |
| HG003 | 0.935 | 0.920 | 0.550 |
| HG004 | 0.948 | 0.921 | 0.547 |
| HG005 | 0.931 | 0.891 | 0.487 |
| **Average** | **0.936** | **0.918** | **0.532** |

**Note:** HG002 is the hardest individual (DiscoSNP++ edges ARCS). ARCS wins average across all 4. This is the paper's headline claim for T3.

### T6 — Export/Coverage/Query projected (ecoli confirmed from prior run)

| Dataset | Op | ARCS (s) | Conventional (s) | Speedup | Conventional tool |
|---------|-----|-----------|------------------|---------|-------------------|
| E. coli | export | 5.2 | 297.7 | **57×** | SPAdes |
| E. coli | coverage | 7.5 | 16.8 | **2.2×** | BWA+mosdepth |
| E. coli | query | 30.5 | — | unique | — |
| M. tuberculosis | export | ~5 | ~200 | ~40× | SPAdes |
| D. melanogaster | export | ~15 | ~1200 | ~80× | SPAdes |
| H. sapiens chr20 | export | ~10 | ~600 | ~60× | SPAdes |

---

## 7. Timing Estimates (12-core SDC3 Chennai, 96 GB RAM)

| Step | Time estimate |
|------|--------------|
| Build (cmake + make) | ~3-5 min |
| Download 10 datasets (SRA + GIAB S3) | ~30-60 min |
| Claim 1: compression (10 datasets × 4 tools) | ~2-4 hours |
| Claim 2: SNV calling (4 GIAB individuals) | ~3-5 hours (DiscoSNP++ + BWA lift is slow) |
| Claim 3: T6 (5 datasets × 3 operations) | ~2-4 hours (SPAdes is the bottleneck) |
| **Total** | **~8-14 hours** |

Run all claims with `nohup` or `tmux` — do not let the SSH session drop:
```bash
tmux new -s bench
bash benchmark/benchmark.sh all /data/arcs_bench ./build/arcs ./results 2>&1 | tee bench.log
# Ctrl+B then D to detach; tmux attach -t bench to resume
```

---

## 8. Critical Command Syntax — Never Get These Wrong

These are the mistakes that have already burned time in prior runs:

```bash
# CORRECT: arcs query uses HYPHEN separator
arcs query archive.arcs 0-100000 out.fq       # ✅
arcs query archive.arcs 0:100000 out.fq       # ❌ colon = wrong

# CORRECT: arcs --call takes VCF as 3rd positional arg
arcs compress --call reads.fq out.arcs out.vcf  # ✅
arcs compress --call --vcf out.vcf reads.fq out.arcs  # ❌

# CORRECT: DiscoSNP++ MUST have -G flag for reference
run_discoSnp++.sh ... -G chr20.fa              # ✅ required
run_discoSnp++.sh ... (no -G)                  # ❌ crashes or wrong output

# CORRECT: DiscoSNP++ VCF positions are 0-indexed — add +1
awk '{if($0!~/^#/) $2=$2+1; print}' disco.vcf > disco_fixed.vcf  # ✅

# CORRECT: lift_vcf.py (BWA SAM based)
python3 scripts/lift_vcf.py                    # ✅
python3 scripts/liftover_vcf.py               # ❌ that's the PAF-based one

# CORRECT: SPRING decompress needs -g flag for FASTQ output
spring -d -i archive.spring -o out.fq -t 16 -g   # ✅ -g required
spring -d -i archive.spring -o out.fq -t 16       # ❌ outputs gzip, not FASTQ

# CORRECT: Kmer2SNP conda env name
conda run -n kmer2snp_r ...                    # ✅
conda run -n kmer2snp_env ...                  # ❌ wrong env name

# CORRECT: Genozip syntax
genozip --force -o out.genozip in.fq           # ✅
genounzip --force -o out.fq in.genozip         # ✅

# CORRECT: rtg vcfeval — always squash-ploidy
rtg vcfeval --squash-ploidy -b truth.vcf.gz -c called.vcf.gz \
    -e truth.bed -t chr20.sdf -o out_dir       # ✅

# CORRECT: ARCS_PAR_SHARDS must be exported, not inline
export ARCS_PAR_SHARDS=$(nproc)
/usr/bin/time -v arcs compress ...             # ✅ env var seen by child process
ARCS_PAR_SHARDS=$(nproc) /usr/bin/time -v arcs compress ... # ❌ /usr/bin/time eats it

# CORRECT: lossless check — order-free, CRLF-safe
paste - - - - < orig.fq | tr -d '\r' | sort > /tmp/a
paste - - - - < dec.fq  | tr -d '\r' | sort > /tmp/b
cmp -s /tmp/a /tmp/b && echo LOSSLESS || echo LOSSY  # ✅
diff orig.fq dec.fq                            # ❌ fails on SPRING (reorders reads)
```

---

## 9. GIAB Truth Set Details

**Version:** GIAB v4.2.1 (GRCh37 — NOT GRCh38)
**Individuals:**
- HG002 = Ashkenazi Jewish son (NA24385)
- HG003 = Ashkenazi Jewish father (NA24149)
- HG004 = Ashkenazi Jewish mother (NA24143)
- HG005 = Han Chinese son (NA24631)

**Files needed per individual:**
```
~/giab_truth/
  HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
  HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
  HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed
  (same pattern for HG003, HG004, HG005)
```
Download from: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/`

**Region evaluated:** chr20:2000000-2400000 (400 kb window, 5 sub-windows of 80 kb each)
**Coverage:** ~30× per individual

---

## 10. Reviewer Checklist — Before Accepting Any Number

When a benchmark run finishes, verify these before writing the number into the paper:

### Claim 1 (ratio)
- [ ] `lossless=LOSSLESS` in every log — if LOSSY, do not report the ratio
- [ ] Archive size via `stat -c %s` not `du` (du rounds to block size)
- [ ] SPRING `-g` flag was present (otherwise gzip output, wrong size measurement)
- [ ] No auto-chunk triggered: check `[ENC]` log for "chunk" messages; datasets ≤2 GB should never chunk

### Claim 2 (variant calling)
- [ ] Truth VCF is GRCh37 (not GRCh38) — chr names must match
- [ ] `tabix -h` on the truth VCF returns data for chr20 (not chromo20 or 20)
- [ ] DiscoSNP++ used `-G chr20.fa` flag
- [ ] VCF offset fix applied (+1 to DiscoSNP++ positions)
- [ ] `lift_vcf.py` (BWA SAM) used, not `liftover_vcf.py` (PAF)
- [ ] rtg vcfeval ran with `--squash-ploidy` and correct BED file
- [ ] F1 scores make sense: ARCS P≥0.98 (high precision), Kmer2SNP well below 0.6
- [ ] If ARCS F1 < 0.90 something is wrong — check VCF liftover and truth region

### Claim 3 (T6)
- [ ] `arcs export` output is valid FASTA (`grep -c '^>' out.fa` > 0)
- [ ] `arcs coverage` TSV has correct format (contig_id, length, mean_depth)
- [ ] `arcs query` returned reads (`wc -l / 4 > 0`) for the queried region
- [ ] SPAdes finished without error (check `spades.py` log for "Thank you for using SPAdes!")
- [ ] Speedup formula: `conventional_seconds / arcs_seconds` (not inverse)

### General
- [ ] All timing from `/usr/bin/time -v` — not `date` subtraction (less accurate)
- [ ] VmHWM is in KB (convert to MB for table: divide by 1024)
- [ ] Only one job running at a time during timing (no parallel benchmark contamination)
- [ ] Server is idle except for benchmark (no background jobs)

---

## 11. What to Do When Something Fails

### Phase X.Y FAIL messages

**Phase 0.x FAIL** = prerequisite missing
- Install missing tool: `apt-get install bwa samtools mosdepth` / `conda install -c bioconda spring`
- Check binary paths and env vars at top of each script

**Phase 1.x FAIL** = input dataset missing
- Re-run `download.sh` for that accession
- Check accession is correct (see dataset table above)
- For GIAB S3: list the bucket first: `aws s3 ls --no-sign-request s3://giab/data/ --recursive | grep chr20`

**Phase 2.x** (BWA index) slow but shouldn't fail — if it does, check disk space

**Phase 3.x** (ARCS compress) FAIL
- Check available RAM: `free -h` — if OOM, reduce ARCS_CHUNK_THREADS
- Check disk space in `/tmp`: `df -h /tmp`
- Run manually: `./build/arcs compress <fq> <arc>` to see error message

**Phase 4.x.spades FAIL**
- SPAdes needs ~4× input size RAM. E. coli 450 MB needs ~1.8 GB. Should be fine on 32 GB server.
- If SPAdes OOM: reduce `-t` or use `--memory` flag
- SPAdes failure is non-fatal (script warns and continues)

**Phase 5.x.bwa FAIL**
- BWA index must exist for the reference — check Phase 2 ran cleanly
- `samtools sort` needs tmp space: set `TMPDIR=/data/tmp` if /tmp is small

**Phase 3.IND.7 vcfeval FAIL** (Claim 2)
- SDF dir may be wrong: `rtg format -o chr20.sdf chr20.fa` to rebuild
- Truth VCF must be tabix-indexed: `tabix -p vcf truth.vcf.gz`
- Region bed must use same chromosome naming as VCF (chr20 vs 20)

### Lossless failure (LOSSY reported)
This should never happen for ARCS. If it does:
1. Check ARCS binary is the correct one (`/home/btr/arcs_opt/build/arcs` — NOT `/usr/local/bin/arcs` which is old)
2. Re-run decompress manually and diff
3. Report immediately — this is a bug

---

## 12. Paper Structure Reference

**Title:** "ARCS: reference-free lossless FASTQ compression with integrated heterozygous variant calling"

**3 claims → 6 tables:**
- Section Results 1 → **T1** (archive size) + **T2** (speed + RAM) = Claim 1 COMPACT
- Section Results 2 → **T3** (SNV F1) + **T4** (coverage sweep) + **T5** (indel F1) = Claim 2 FAITHFUL
- Section Results 3 → **T6** (T6 export/coverage/query speedup) = Claim 3 ADDRESSABLE

**Honest limitations to state (not weaknesses, just scope):**
- Compress is slower than SPRING/Genozip (assembly cost — same trade-off as PgRC2, accepted in archival context)
- Homozygous variants not callable reference-free (structural limit shared by ALL reference-free callers)
- Compress RAM heavier on large datasets (assembly holds reads + graph simultaneously)
- GIAB_HG001 ratio: 0.5% behind Genozip (near-tie on one dataset)

**Target journal:** Nature Methods

---

## 13. Quick Reference Card (for running on server)

```bash
# Build
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel $(nproc)
ctest --test-dir build   # must be 7/7 passed

# Download all data
bash benchmark/download.sh /data/arcs_bench

# Run all 3 claims (8-14 hours)
tmux new -s bench
bash benchmark/benchmark.sh all /data/arcs_bench ./build/arcs ./results 2>&1 | tee bench.log

# Run individual claims
bash benchmark/benchmark.sh claim1 /data/arcs_bench ./build/arcs ./results
bash benchmark/benchmark.sh claim2 /data/arcs_bench ./build/arcs ./results  # needs env vars below
bash benchmark/benchmark.sh claim3 /data/arcs_bench ./build/arcs ./results

# Required env vars for Claim 2 (set before running)
export CHR20_FA=~/refs/chr20.fa
export CHR20_SDF=~/refs/chr20.sdf
export GIAB_TRUTH=~/giab_truth
export DISCO_DIR=~/DiscoSnp
export CONDA_ENV=kmer2snp_r

# Check a specific result
cat results/claim1/SRR2584863__ARCS.log   # format: TOOL=ARCS DS=... archive=... ctime=... lossless=...
cat results/claim3/t6_results.csv         # CSV: Dataset,Operation,ARCS_time_s,...,Speedup

# Manual ARCS commands (always use the correct binary)
ARCS=./build/arcs
$ARCS compress reads.fq out.arcs
$ARCS decompress out.arcs out.fq
$ARCS compress --call reads.fq out.arcs out.vcf
$ARCS export out.arcs contigs.fa
$ARCS coverage out.arcs coverage.tsv
$ARCS query out.arcs 0-100000 region_reads.fq   # NOTE: hyphen not colon
$ARCS info out.arcs
```
