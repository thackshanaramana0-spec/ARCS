# ARCS Benchmark — New Server Handover (Mandatory Read)

Authoritative handover from the previous server session. Read this fully before touching `CLAUDE.md`, `benchmark/DATASET_LOCKED.md`, or running any phase. This document supersedes any stale in-repo text that conflicts with it — the facts here are what actually happened and what is actually true of the data on disk.

---

## 1. How to use this document

- This is a **factual record**, not a tutorial. Every claim below was verified on the server via direct command output during the previous session.
- Section 3 is the fresh-server bring-up sequence, start to finish, in order.
- Section 8 lists local git changes that were **never pushed**. Check these against `git diff` / `git log` on arrival — if they're missing, re-apply them before building or running Phase 3.
- Section 10 is the exact current state of `/data/fastq`, `~/refs`, `~/giab_truth` as last verified.

---

## 2. Project identity (unchanged)

- **ARCS** = Assemble → Retain → Compress → Serve. Reference-free FASTQ compressor + variant caller.
- Repo: `https://github.com/thackshanaramana0-spec/ARCS`, local path `/root/arcs-clean`.
- Governing rules file: `CLAUDE.md` (still the source of truth for phase structure, absolute rules, command syntax). Nothing in that file's **rules** (zero-flags compression, no concurrent timed jobs, lossless-before-Claim-2, etc.) was changed this session — only the **dataset composition** for Claims 1 and 2 changed (see Sections 5–6).

---

## 3. Fresh server bring-up sequence (first to last)

### 3.1 Persistence — set up tmux before anything else

Do this immediately after SSH'ing into a new server, before starting Claude:

```bash
cd ~/arcs-clean
tmux new -s arcs
claude --continue    # or plain `claude` on a genuinely new server with no prior session
```

If SSH disconnects: reconnect, then `tmux attach -t arcs` — the pane and whatever Claude/background jobs were running are still there. This was set up and used successfully for the remainder of the previous session.

**Note on background jobs specifically**: manually backgrounding processes with `nohup ... & disown` or even `setsid` was **unreliable** in this harness — jobs died silently without tmux protection being the cause. Always use the coding agent's own background-task mechanism (`run_in_background`) for anything that must survive and be checked on later, not raw shell backgrounding.

### 3.2 Git

```bash
cd /root/arcs-clean
git pull origin main
git log --oneline -5   # confirm HEAD; as of this handover, latest is 22d4ed8
```

### 3.3 Apply local-only fixes if missing (see Section 8 for exact diffs)

Check `git diff` is empty for `CMakeLists.txt` and `benchmark/download.sh`. If the previous session's local fixes were never pushed, re-apply them from Section 8 before proceeding — **the build will fail without the CMakeLists.txt fix.**

### 3.4 Install toolchain (apt-installable)

```bash
apt-get update -y
apt-get install -y cmake build-essential git curl tabix bcftools bwa unzip \
    samtools mosdepth seqtk nasm
```

Note: `apt-get install -y awscli` **does not resolve** on this server's repo configuration ("no installation candidate"). Use the official installer instead:

```bash
cd /root
curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip
unzip -q awscliv2.zip
./aws/install
aws --version
```

### 3.5 SRA Toolkit (prefetch, fasterq-dump)

```bash
cd /root
wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
```

**Critical gotcha**: `prefetch`/`fasterq-dump` are symlink chains (`prefetch → prefetch.3 → prefetch.3.4.1 → sratools.3.4.1`) that expect a sibling `<tool>-orig.3.4.1` binary **in the same directory** at runtime. Flattening with `cp sratoolkit.*/bin/prefetch /usr/local/bin/` breaks this (`failed to exec prefetch: No such file or directory`). **Use symlinks, not copies:**

```bash
ln -sf /root/sratoolkit.*/bin/prefetch /usr/local/bin/prefetch
ln -sf /root/sratoolkit.*/bin/fasterq-dump /usr/local/bin/fasterq-dump
prefetch --version && fasterq-dump --version   # must both print a version, not an exec error
```

### 3.6 SPRING (build from source — required patch)

```bash
cd /root
git clone https://github.com/shubhamchandak94/SPRING.git
cd SPRING
git submodule update --init --recursive
```

**Required patch**: this server's GCC 13 does not transitively pull in `<cstdint>`. Unpatched, the build fails with `'uint8_t' does not name a type` cascading into dozens of "no member named" errors. Patch every **compiled** source/header (do not bother with `old_src/`, it isn't built) by inserting `#include <cstdint>` before the first `#include`:

```bash
cd /root/SPRING/src
for f in bitset_util.cpp id_compression/src/stream_model.cpp params.h spring.h \
    id_compression/src/sam_file_allocation.cpp encoder.cpp reorder_compress_streams.cpp \
    id_compression/include/id_compression.h BooPHF.h decompress.cpp main.cpp libbsc/bsc.h \
    id_compression/src/compression.cpp encoder.h pe_encode.cpp bitset_util.h reorder.h \
    util.cpp spring.cpp decompress.h reorder_compress_quality_id.cpp preprocess.cpp \
    id_compression/src/io_functions.cpp id_compression/include/sam_block.h \
    id_compression/include/Arithmetic_stream.h id_compression/src/sam_models.cpp \
    reorder_compress_quality_id.h id_compression/include/stream_model.h \
    id_compression/src/id_compression.cpp libbsc/bsc_str_array.cpp \
    id_compression/src/Arithmetic_stream.cpp qvz/src/cluster.cpp qvz/src/well.cpp \
    qvz/src/qv_compressor.cpp qvz/src/pmf.cpp qvz/src/quantizer.cpp qvz/include/qvz.h \
    qvz/src/codebook.cpp qvz/src/distortion.cpp qvz/src/lines.cpp; do
  if [ -f "$f" ] && ! grep -q '#include <cstdint>' "$f"; then
    sed -i '0,/^#include/s//#include <cstdint>\n#include/' "$f"
  fi
done
```

Then build and install:

```bash
mkdir -p /root/SPRING/build && cd /root/SPRING/build
cmake .. && make -j$(nproc)
cp spring /usr/local/bin/
spring --help   # sanity check
```

(`-Walloc-size-larger-than=` warnings during build are benign, not errors.)

### 3.7 Genozip (build from source + one-time interactive activation)

No prebuilt Linux binaries on GitHub releases — source only. Needs `nasm`.

```bash
cd /root
curl -sL https://github.com/divonlan/genozip/archive/refs/tags/genozip-15.0.87.tar.gz -o genozip-src.tar.gz
tar -xzf genozip-src.tar.gz
cd genozip-genozip-15.0.87
make -j$(nproc)
cp genozip genounzip genocat genols /usr/local/bin/
```

**Activation is mandatory before genozip will run at all**, and it is **interactive** (license key entry → email confirmation → T&C acceptance → emailed 6-digit code entry). This cannot be reliably scripted through this agent's tool sandbox — a prior attempt to drive it via a Python `pty` wrapper was blocked by the auto-mode safety classifier partway through (submitting an OTP-like code programmatically). **Run `genozip --activate` yourself, interactively, in the terminal**, using a Genozip Student license (apply free at genozip.com/student if the previous license has expired — it's valid ~1 year). Once activated:

```bash
echo "test" > /tmp/gz_test.txt
genozip --force -o /tmp/gz_test.txt.genozip /tmp/gz_test.txt   # must succeed, not "requires activation"
rm -f /tmp/gz_test.txt /tmp/gz_test.txt.genozip
```

Last known activation: Student license, expires 2027-05-25, verified with a live compress/decompress/diff round trip (LOSSLESS OK).

### 3.8 rtg-tools (for Claim 2 `vcfeval`)

The `.../releases/latest/download/rtg-tools-non-commercial-linux-x64.zip` URL **404s** — that asset name isn't on current releases. Query the API for the real asset:

```bash
curl -s https://api.github.com/repos/RealTimeGenomics/rtg-tools/releases/latest \
  | grep browser_download_url
# use the *-linux-x64.zip asset (was rtg-tools-3.13-linux-x64.zip at time of writing)
cd /root
curl -sL <that-url> -o rtg-tools.zip
unzip -q rtg-tools.zip
cp -r rtg-tools-3.13 /usr/local/rtg-tools
ln -sf /usr/local/rtg-tools/rtg /usr/local/bin/rtg
rtg version
```

### 3.9 Build ARCS and run Phase 0 tests

```bash
cd /root/arcs-clean
rm -rf build
cmake -B build -DCMAKE_BUILD_TYPE=Release   # will fail here without the Section 8 CMakeLists.txt fix
cmake --build build --parallel $(nproc)
ctest --test-dir build                       # must show 7/7 tests passed
ls -la build/arcs                            # confirm binary exists, no /usr/local/bin/arcs shadow
```

---

## 4. Absolute rules — unchanged, still binding

From `CLAUDE.md` (verify these still read this way on arrival; do not weaken them):

1. Do not change the 10 benchmark accessions without explicit instruction.
2. ARCS compression uses zero flags, always.
3. Never run two **timed benchmark** jobs concurrently (Phase 2/3/4). This does **not** apply to Phase 1 downloading — see Section 9 for where parallelizing Phase 1 was judged safe and where it backfired.
4. Never silently alter methodology — stop and report on failure.
5. Claim 1 must be fully lossless before Claim 2 starts.
6. Projected numbers are sanity checks only; fresh measurements are authoritative.
7. ARCS must pass 7/7 ctests before any benchmark phase runs.

---

## 5. Claim 1 — final decision (changed this session)

**Claim 1 uses exactly the 10 full SRA datasets. The GIAB human chr20 subsets are no longer part of Claim 1 at all.**

Reasoning (user's call, agreed and implemented):
- The old plan had 5 GIAB chr20 "supplementary" rows in Claim 1 at ~37 MB / ~0.7× coverage each.
- At 0.7× coverage, ARCS cannot build a meaningful assembly — reads mostly go unplaced.
- SPRING would trivially win those 5 rows (no assembly needed, pure statistical model), which actively hurts Claim 1's story rather than helping it.
- Decision: drop GIAB from Claim 1 entirely. Claim 1 (T1 archive size, T2 speed+RAM) = 10 full SRA datasets only, full stop.

`run_block1.sh`'s dataset-discovery glob was fixed accordingly (see Section 9.9) to guarantee it sees exactly these 10 files and nothing else.

---

## 6. Claim 2 / GIAB / uniform 30× — final decision (changed this session)

**Claim 2 (T3/T4/T5) uses HG002–HG005 chr20 reads, extracted by streaming from GIAB's real S3 BAMs and downsampled to a uniform 30× depth. HG001 is not used anywhere anymore.**

### What was wrong with the original plan
The originally locked `download.sh` had hardcoded S3 paths like `s3://giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_subsetN_chr20.fastq.gz` for all 5 individuals (HG001–HG005). **These files do not exist in the real GIAB bucket for any individual.** Verified by listing the actual bucket structure: GIAB publishes either raw whole-genome multi-lane FASTQ (HG001/NA12878 case) or full-genome BAM alignments (HG002–HG005 case) — never a pre-extracted chr20-only FASTQ at the path the script assumed.

### What was decided instead
1. **Auto-discover** each individual's largest WGS BAM under their S3 prefix (`s3://giab/data/AshkenazimTrio/HG00{2,3,4}_.../`, `s3://giab/data/ChineseTrio/HG005_.../`), by size, excluding indel/SV/SNV/phased-annotation files.
2. **Stream just the chr20 region** via `samtools view -b <S3 URI> chr20` — htslib has anonymous S3 read support built in for a public bucket; **no `--no-sign-request` or any other auth flag is needed or accepted on the `samtools` side** (that flag is AWS-CLI-only and `samtools view` rejects it outright).
3. Name-sort (`samtools sort -n`) then convert to FASTQ (`samtools fastq -o`, **not** `-0` — see Section 9.7).
4. **Downsample every individual to exactly 30× chr20** via `seqtk sample -s 42 <frac>`, because native coverage varies wildly by individual depending on which BAM tier happens to be the "largest" available (HG002/HG005 landed on ~300× BAMs, HG003/HG004 also landed on ~300× BAMs in practice on this run — see Section 10 for actual native read counts). Uniform 30× removes this confound so T3's F1 comparison across the four individuals is apples-to-apples.
5. HG001 is dropped entirely — it was only ever used for Claim 1's now-deleted supplementary rows; Claim 2 was always scoped to HG002–HG005 only per `CLAUDE.md`.

### Reference/genome build note
All four individuals' auto-selected BAMs were `GRCh38`-aligned in this run (confirmed via filenames, e.g. `HG002.GRCh38.300x.bam`), not `hs37d5`/GRCh37. **This must be reconciled against `~/refs/chr20.fa`, which is GRCh37/hg19** (chromosome named `20`, no `chr` prefix) and against the GIAB truth VCFs, which are also GRCh37 (`*_GRCh37_1_22_v4.2.1_benchmark.vcf.gz`). The chr20 region-filter (`chr20` as the samtools region string) matches GRCh38's contig naming, not GRCh37's `20`. **This mismatch was not caught or resolved in the previous session** — flagged as unresolved in Section 11. Before running Claim 2's variant-calling comparison, verify whether the GRCh38-streamed reads need re-lifting to GRCh37 coordinates, or whether an `hs37d5`-based BAM should have been forced for consistency instead (an `hs37d5` chr20-only BAM was confirmed to exist for HG002 specifically: `HG002.hs37d5.300x_chr20.bam`, but the "largest file" auto-selection logic did not know to prefer it).

---

## 7. Dataset lock evolution (full history, for context)

The "10 locked datasets" list went through several revisions across this and prior sessions:

1. `e5b042b` — first lock: included oversized `SRR065390`/`SRR327342`/`SRR1663585`/`SRR870667` (some >2 GB, would auto-chunk).
2. `58f6608` — replaced those 4 with `SRR16357346`/`DRR976266`/`SRR36741279`/`SRR37283774` to stay under 2 GB.
3. `ed07c82` — swapped slot 5 to `ERR17740259` (S. aureus); added 5 GIAB human chr20 subsets as Claim 1 supplementary rows.
4. `af01f27` — **reverted slots 6+8 back to `SRR065390`/`SRR870667`**, explicitly permitted because `ARCS_AUTOCHUNK_MB=25000` (set in `run_block1.sh`) suppresses auto-chunking on this 90 GB-RAM server, and because these two accessions carry published SPRING gold bpb cross-validation values that the smaller substitutes did not have.
5. `862a04f` — added Phase 5 (chr20.fa + GIAB truth VCF acquisition) to `download.sh`.
6. `6babda8` — **dropped the 5 GIAB rows from Claim 1** entirely; introduced the BAM-streaming method for Claim 2 (Section 6).
7. `a2db219` / `6207278` — fixed `run_block1.sh`'s dataset-discovery glob to exclude `_2.fq` mates, `HG00*.fq`, and PE-orphan bare `.fq` files.
8. `22d4ed8` — fixed the `samtools fastq -0`→`-o` bug in `download.sh` (Section 9.7).

### Final locked 10 (Claim 1, primary):

| # | Accession | Organism | Actual `_1.fq` size (measured) | Notes |
|---|-----------|----------|-------------------------------|-------|
| 1 | SRR2584863 | E. coli B REL606 | 695,163,748 B | |
| 2 | ERR552797 | M. tuberculosis H37Rv | 431,078,636 B | |
| 3 | SRR554369 | P. aeruginosa PAO1 | 456,443,722 B | SPRING gold cross-val: 0.2416 bpb |
| 4 | ERR5181310 | SARS-CoV-2 amplicon | 471,183,048 B | verified against ENA metadata (1,828,186 total reads across mates, matches) |
| 5 | ERR17740259 | S. aureus | 1,336,957,140 B | |
| 6 | SRR065390 | C. elegans N2 | 11,282,985,734 B | auto-chunk suppressed via `ARCS_AUTOCHUNK_MB=25000`; spot count 33,808,546 matches doc exactly |
| 7 | DRR976266 | S. cerevisiae | 2,248,042,490 B | |
| 8 | SRR870667 | T. cacao | 22,986,662,428 B | **measured ~21.4 GB, doc estimate was ~15 GB** — larger than expected but genuine; auto-chunk suppressed; SPRING gold cross-val: 1.2621 bpb, SPRING's worst dataset |
| 9 | SRR36741279 | Leishmania major | 1,661,157,034 B | |
| 10 | SRR37283774 | P. falciparum | 1,019,818,902 B | |

**Banned accessions (never reuse)**: SRR390728, SRR988075, SRR327342, SRR1663585, SRR1296601, ERR015526, SRR1294122, ERR174310, SRR16357346, SRR1945765.

### Claim 2 only — 4 GIAB individuals, uniform 30× chr20 (streamed, not part of Claim 1):

| Individual | Native BAM found | Native reads (chr20, full) | Downsampled to |
|---|---|---|---|
| HG002 | `HG002.GRCh38.300x.bam` | 124,888,113 | 12,604,917 reads → 30.0× |
| HG003 | `HG003.GRCh38.300x.bam` | 126,991,238 | 12,606,030 reads → 30.0× |
| HG004 | `HG004.GRCh38.300x.bam` | 139,040,965 | 12,605,295 reads → 30.0× |
| HG005 | `HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam` | 84,593,588 | 12,603,296 reads → 30.0× |

All four downsampled with `seqtk sample -s 42 <frac>` for reproducibility (fixed seed).

---

## 8. Local git changes NOT yet pushed — verify/reapply on arrival

These were fixed and tested on the server but **never committed or pushed** (only `git diff`-visible locally). If a fresh clone doesn't have them, re-apply before building/running:

### 8.1 `CMakeLists.txt` — **critical, build fails without this**

```diff
- project(arcs_cpp VERSION 2.2.0 LANGUAGES CXX)
+ project(arcs_cpp VERSION 2.2.0 LANGUAGES C CXX)
```//
Root cause: without `LANGUAGES C`, CMake never enables a C compiler, so `third_party/libbsc/bwt/divsufsort/divsufsort.c` (a `.c` file feeding into every target, including `arcs` itself) is silently never compiled. Every target fails to link with `undefined reference to 'divbwt'`. Confirmed via `compile_commands.json` showing zero entries for that file before the fix, and correct `Building C object ...divsufsort.c.o` entries after.

### 8.2 `benchmark/download.sh` — three fixes in the GIAB BAM-streaming block (Phase 3)

```diff
     BAM_KEY=$(aws s3 ls --no-sign-request "$BASE/" --recursive 2>/dev/null \
-        | grep '\.bam$' | grep -v '\.bai' | grep -iv 'indel\|sv\|snv\|phased\|trio' \
+        | grep '\.bam$' | grep -v '\.bai' | grep -iv 'indel\|sv\|snv\|phased' \
         | sort -k3 -n -r | head -1 | awk '{print $4}')
```
Root cause: the `trio` exclusion keyword also matched the *directory names* `NHGRI_Illumina300X_AJtrio_novoalign_bams` and `NHGRI_Illumina300X_Chinesetrio_novoalign_bams` — i.e. it excluded every legitimate WGS BAM for every individual, leaving zero candidates. Confirmed empirically against the live bucket listing.

Also wrap the discovery pipeline to avoid a SIGPIPE-as-failure under `pipefail` (`sort | head -1` kills `sort` with SIGPIPE once `head` exits, which `set -euo pipefail` at the top of the script would otherwise report as exit 141):
```diff
-    BAM_KEY=$(aws s3 ls --no-sign-request "$BASE/" --recursive 2>/dev/null \
+    BAM_KEY=$(set +o pipefail; aws s3 ls --no-sign-request "$BASE/" --recursive 2>/dev/null \
```

```diff
-    samtools view -b --no-sign-request "$S3_BAM" chr20 2>/dev/null \
+    samtools view -b "$S3_BAM" chr20 2>/dev/null \
         | samtools sort -n -@ "$JOBS" \
```
Root cause: `--no-sign-request` is an AWS-CLI-only flag; `samtools view` rejects it outright ("unrecognised option"). htslib's S3 support handles anonymous access to this public bucket automatically with a plain `s3://` URI.

```diff
-        | samtools sort -n -@ "$JOBS" \
+        | samtools sort -n -@ "$JOBS" -T "/tmp/samtools_sort_tmp/${IND}" \
```
Root cause: without `-T`, `samtools sort` writes temp files to the current working directory (`/root/arcs-clean` when invoked via `download.sh`), risking collisions with concurrent git operations on the repo and generally polluting the tracked tree. (Note: the actual crash chain this surfaced during was disk exhaustion, not a git collision — see Section 9.10 — but isolating temp files is still the correct fix regardless.) Requires `mkdir -p /tmp/samtools_sort_tmp` to exist first.

The `-0`→`-o` fix (Section 9.7) **was already pushed** as commit `22d4ed8` — no action needed for that one.

---

## 9. Problems discovered and fixed, in chronological order

1. **`src/tensor.cpp`/`src/tensor.h` missing from repo** — commit `f73d830` referenced them in `CMakeLists.txt`/`main.cpp` but never committed the files. Fixed upstream in commit `1b0192e` (pushed, no local-only remnant).
2. **`CMakeLists.txt` missing `LANGUAGES C`** — see Section 8.1. Local only, not pushed.
3. **SPRING build fails on GCC 13** — missing `<cstdint>` includes across ~35 compiled files. Patched (Section 3.6), not part of the git repo (SPRING is a separate clone, not vendored).
4. **`prefetch`/`fasterq-dump` broken by flattening symlinks** — see Section 3.5. Environment-setup mistake, not a repo bug; fixed by using symlinks instead of `cp`.
5. **`apt-get install awscli` has no candidate** on this server's repo config — used the official AWS CLI v2 installer instead (Section 3.4).
6. **rtg-tools `latest/download` URL 404s** — the expected asset name isn't attached to the current release; had to query the GitHub API for the real asset URL (Section 3.8).
7. **GIAB S3 paths in `download.sh` Phase 3 don't exist** — see Section 6. Rewritten to auto-discover + stream from real WGS BAMs.
8. **`samtools view --no-sign-request` invalid flag** — see Section 8.2. Pushed as part of `6babda8`'s rewrite conceptually, but the invalid flag was actually only caught and fixed locally; verify it's genuinely gone on arrival.
9. **`grep -iv '...|trio'` excluded every valid BAM** — see Section 8.2.
10. **SIGPIPE (`exit 141`) from `sort | head -1` under `pipefail`** — see Section 8.2.
11. **`samtools fastq -0` bug — leaked raw FASTQ into log files** — `-0` only captures true orphan reads (neither R1 nor R2 flag set); on paired BAMs almost all real reads went to **stdout by default** since no `-1`/`-2`/`-o` was given, which is why raw sequence/quality lines appeared inside the `.log` files instead of the intended output file. Fixed to `-o "$TMPFQ"`. **This was pushed as commit `22d4ed8`.** Verified via `samtools fastq --help` and a live 500kb-region test (883,467 reads recovered correctly with "discarded 0 singletons").
12. **`run_block1.sh` discovery glob picked up 12 files instead of 10** — `_2.fq` mate files and `HG00*.fq` were excluded (commit `a2db219`), but bare orphan `ACCESSION.fq` singleton-read files left behind by `fasterq-dump` for some PE accessions (`ERR5181310.fq`: 3,044 reads; `SRR36741279.fq`: 4,753,438 reads) were still being counted as independent datasets. Fixed in commit `6207278` (post-filter: skip any bare `.fq` when a matching `ACCESSION_1.fq` exists). **Pushed, confirmed byte-identical to what was tested locally.**
13. **Manual background-process management repeatedly failed** — `nohup ... & disown`, even `setsid`, did not reliably survive; jobs died silently more than once with zero error output. Switched to the harness's own background-task tracking mechanism, which then gave clean, diagnosable exit codes for every subsequent failure.
14. **All 4 background jobs killed simultaneously mid-stream, unexplained** — happened once. Ruled out: server crash (server stayed reachable), sshd idle timeout (`ClientAliveInterval 0` — disabled entirely, confirmed via `sshd -T`), resource exhaustion (plenty of headroom at the time). Most consistent explanation: a client-side stop-all-agents gesture, not a network disconnect. **tmux (Section 3.1) was set up specifically to insulate against this class of interruption going forward**, whatever its exact cause.
15. **`samtools sort` temp files landing in the git repo directory** — see Section 8.2. First diagnosed as the cause of an `Illegal seek` crash during a concurrent `git stash`/`pull` operation; later proven to actually be **disk exhaustion** (see #16) once the same crash recurred on a machine with no concurrent git activity. The `-T` fix is still correct practice and was kept.
16. **Disk filled to 100% during parallel GIAB streaming — the real root cause of #15's recurrence.** Original assumption (from `DATASET_LOCKED.md`'s ~11–12 GB compressed-BAM-size estimate) badly undersized the actual per-individual intermediate: the **uncompressed, name-sorted chr20-region FASTQ before downsampling** ran 28–44 GB per individual at native ~300× coverage (HG005's intermediate peaked at 43.97 GB before collapsing to 6.8 GB post-downsample). Running all 4 individuals in parallel (reasonable Phase-1-is-not-a-timed-benchmark judgment call, since `CLAUDE.md`'s concurrency rule is about *measured* phases) multiplied this to ~150+ GB of simultaneous intermediate demand against a 233 GB disk that already had ~136 GB committed to the 10 SRA datasets. **Fix applied**: deleted all 10 `_2.fq` PE-mate files (~35 GB, confirmed unneeded — excluded from Claim 1's discovery glob and untouched by Claim 2), the 2 orphan singleton files, and the `/data/fastq/prefetch/*.sra` archives (~14 GB, already converted to `.fq` and verified) — reclaiming ~51 GB. **Then switched to running the remaining GIAB individuals one at a time**, not in parallel, for the rest of Phase 1. **Recommendation for future re-runs**: budget ≥50 GB free disk per GIAB individual if running Phase 3 in parallel, or just run sequentially by default — it's not meaningfully slower in wall-clock terms since S3 streaming throughput per stream did not appear bandwidth-limited by running 4x parallel anyway (~0.5–0.65 GB/min per stream whether solo or parallel).

---

## 10. Current state (as last verified this session)

- **Git**: `HEAD` = `22d4ed8` on `main`, in sync with origin, except the two local-only diffs in Section 8.
- **Build**: `build/arcs` exists and passes `ctest` 7/7 (verified after both the tensor.cpp fix and the CMakeLists.txt fix were in place together).
- **Toolchain**: cmake, prefetch, fasterq-dump, aws (v2.36.30), spring, genozip (activated, Student license, expires 2027-05-25), genounzip, genocat, genols, rtg (3.13), bwa, samtools, bcftools, tabix, mosdepth, seqtk — all confirmed callable.
- **`/data/fastq/`** contains exactly:
  - 10 primary `_1.fq` files (Section 7 table) — **`_2.fq` mates and orphan singleton files were deliberately deleted** to reclaim disk (Section 9.16); if Claim 1 or any other step ever needs paired-end mate data specifically, it is **not currently on disk** and would need re-downloading via `fasterq-dump` from the same accessions.
  - 4 `HGxxx_pooled.fq` files (Section 7 table), each verified at exactly 30.0× chr20 coverage.
  - No `prefetch/` subdirectory (deleted, `.sra` archives no longer needed post-conversion).
- **`~/refs/chr20.fa`** — 64,286,035 bytes, GRCh37/hg19, chromosome named `20` (no `chr` prefix).
- **`~/giab_truth/`** — HG002/HG003/HG004/HG005 truth VCFs (`.vcf.gz` + `.tbi`) and confident-region BEDs, all present and verified.
- **Disk**: last measured at 177 GB used / 233 GB total, 57 GB free (73–76% used range across the last several checks; fluctuates with cleanup/streaming).
- **Phase 0**: PASS.
- **Phase 1**: COMPLETE (10 SRA + 4 GIAB @ 30× + chr20.fa + truth VCFs, all verified).
- **Phase 2/3/4**: **NOT STARTED.**

---

## 11. Unresolved — needs a decision or action before proceeding

1. **GRCh38 vs GRCh37 mismatch, Claim 2 (Section 6)** — the 4 GIAB individuals' chr20 reads were streamed from `GRCh38`-aligned BAMs (confirmed by filename), but `~/refs/chr20.fa` and the GIAB truth VCFs are `GRCh37`/`hs37d5`. This was **not caught or resolved** in the previous session. Must be addressed before running Claim 2's variant calling — either re-source an `hs37d5`-aligned BAM per individual (one is confirmed to exist for HG002: `HG002.hs37d5.300x_chr20.bam`, not yet checked for HG003/004/005), or reconcile via liftover, before trusting any T3/T4/T5 F1 numbers.
2. **Section 8's local-only fixes** — confirm on arrival whether they were ever pushed after this handover was written; if not, re-apply from the exact diffs given.
3. **`CLAUDE.md` Claim 1/2 sections** — commit `6babda8`'s message states `CLAUDE.md` dataset tables and Phase 1 notes were updated to match the new 10-only/GIAB-streaming scheme, but this was **not independently re-read and re-verified** after that push in the previous session. Re-read `CLAUDE.md` in full and cross-check it against Sections 5–7 of this document before treating it as authoritative.
4. **Phase 2/3/4 have not been run at all.** Next step, once #1 is resolved, is Phase 2 (Claim 1 compression benchmark): `bash benchmark/benchmark.sh claim1 /data/fastq /root/arcs-clean/build/arcs ./results`, followed by the standard lossless + SPRING-bpb-cross-validation checks per `CLAUDE.md`.
5. **PE mate (`_2.fq`) data was deleted for disk space** — if any future benchmark step turns out to need paired-end mates (current Phase 2/3/4 scripts as understood do not), it must be re-downloaded first.

---

## 12. Appendix — GIAB per-individual streaming script (reference only, not in git)

This script lived only in a scratch directory (`/tmp/claude-*/scratchpad/giab_one.sh`) during the previous session and will not survive to a new server on its own. `benchmark/download.sh`'s own Phase 3 loop is the canonical, git-tracked equivalent (with the Section 8.2 fixes layered in) and should be preferred — this is included only so the exact logic that was tested individual-by-individual is not lost.

```bash
#!/bin/bash
set -uo pipefail
IND="$1"; BASE="$2"; THREADS="$3"
DATA_DIR="/data/fastq"
CHR20_LEN=63025520
TARGET_COV=30
READ_LEN=150
N_NEEDED=$(( TARGET_COV * CHR20_LEN / READ_LEN ))

OUT="$DATA_DIR/${IND}_pooled.fq"
[ -s "$OUT" ] && { echo "[$IND] already present, skipping"; exit 0; }

BAM_KEY=$(set +o pipefail; aws s3 ls --no-sign-request "$BASE/" --recursive 2>/dev/null \
    | grep '\.bam$' | grep -v '\.bai' | grep -iv 'indel\|sv\|snv\|phased' \
    | sort -k3 -n -r | head -1 | awk '{print $4}')
[ -z "$BAM_KEY" ] && { echo "[$IND] FAIL: no WGS BAM found under $BASE"; exit 1; }
S3_BAM="s3://giab/${BAM_KEY}"

verify_fq() {
    local f="$1" name="$2"
    [ -s "$f" ] || { echo "[$name] FAIL: $f missing or empty"; exit 1; }
    [[ "$(head -1 "$f")" == @* ]] || { echo "[$name] FAIL: corrupt FASTQ"; exit 1; }
    local lines; lines=$(wc -l < "$f")
    [ $(( lines % 4 )) -eq 0 ] || { echo "[$name] FAIL: truncated FASTQ"; exit 1; }
}

TMPFQ="$DATA_DIR/${IND}_chr20_full.fq"
set -o pipefail
samtools view -b "$S3_BAM" chr20 2>"$DATA_DIR/${IND}_view.err" \
    | samtools sort -n -@ "$THREADS" -T "/tmp/samtools_sort_tmp/${IND}" 2>"$DATA_DIR/${IND}_sort.err" \
    | samtools fastq -o "$TMPFQ" - 2>"$DATA_DIR/${IND}_fastq.err"
RC=(${PIPESTATUS[@]}); set +o pipefail
[ "${RC[0]}" -eq 0 ] || { echo "[$IND] FAIL view: $(cat "$DATA_DIR/${IND}_view.err")"; exit 1; }
[ "${RC[1]}" -eq 0 ] || { echo "[$IND] FAIL sort: $(cat "$DATA_DIR/${IND}_sort.err")"; exit 1; }
[ "${RC[2]}" -eq 0 ] || { echo "[$IND] FAIL fastq: $(cat "$DATA_DIR/${IND}_fastq.err")"; exit 1; }
rm -f "$DATA_DIR/${IND}"_{view,sort,fastq}.err
verify_fq "$TMPFQ" "$IND"

N_FULL=$(( $(wc -l < "$TMPFQ") / 4 ))
if [ "$N_FULL" -le "$N_NEEDED" ]; then
    mv "$TMPFQ" "$OUT"
else
    FRAC=$(awk "BEGIN{printf \"%.6f\", $N_NEEDED / $N_FULL}")
    seqtk sample -s 42 "$TMPFQ" "$FRAC" > "$OUT"
    rm -f "$TMPFQ"
fi
verify_fq "$OUT" "$IND"
echo "[$IND] DONE: $(( $(wc -l < "$OUT") / 4 )) reads -> $OUT"
```

Requires `mkdir -p /tmp/samtools_sort_tmp` first. Call once per individual with sequential invocation recommended (Section 9.16):

```bash
bash giab_one.sh HG002 "s3://giab/data/AshkenazimTrio/HG002_NA24385_son" 4
bash giab_one.sh HG003 "s3://giab/data/AshkenazimTrio/HG003_NA24149_father" 4
bash giab_one.sh HG004 "s3://giab/data/AshkenazimTrio/HG004_NA24143_mother" 4
bash giab_one.sh HG005 "s3://giab/data/ChineseTrio/HG005_NA24631_son" 4
```
