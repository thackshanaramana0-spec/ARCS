#!/bin/bash
# Download all 10 datasets for ARCS benchmark
#
# Usage:
#   bash benchmark/download.sh <DATA_DIR> [JOBS]
#
# DATA_DIR : destination directory for FASTQ files
# JOBS     : parallel sra-tools jobs (default: 4)
#
# Tools required on PATH:
#   prefetch, fasterq-dump  (SRA Toolkit >= 3.0)
#   aws                      (AWS CLI, for GIAB S3)
#   pigz or gzip             (for .fastq.gz → .fq conversion)
#
# After download, each file is verified with: head -4 | check first char = '@'
# All outputs land in DATA_DIR as plain .fq files (no gzip — ARCS reads uncompressed)

set -euo pipefail

DATA_DIR="${1:?Usage: $0 <DATA_DIR> [JOBS]}"
NPROC=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
JOBS="${2:-$NPROC}"
mkdir -p "$DATA_DIR"

# ── Phase helpers ─────────────────────────────────────────────────────────────
_PH=""
phase() { _PH="$1"; echo ""; echo "[Phase $_PH] $2"; }
pdone() { echo "[Phase $_PH] DONE — $*"; }
pfail() { echo "[Phase $_PH] FAIL — $*"; exit 1; }
pskip() { echo "[Phase $_PH] SKIP — $*"; }
pinfo() { echo "[Phase $_PH]   → $*"; }

verify_fq() {
    local f="$1"
    local name="$2"
    if [ ! -s "$f" ]; then
        pfail "$name: file missing or empty — $f"
    fi
    local first
    first=$(head -1 "$f" 2>/dev/null || true)
    if [[ "$first" != @* ]]; then
        pfail "$name: first line is not '@' — corrupt FASTQ? first: $first"
    fi
    local n
    n=$(( $(wc -l < "$f") / 4 ))
    local sz
    sz=$(stat -c %s "$f")
    pinfo "$name — $n reads, ${sz} bytes — OK"
}

# ── Dataset table ─────────────────────────────────────────────────────────────
# 10 full SRA datasets (complete unsubsampled runs)
# Name               Accession       Organism                    Approx _1 size  Notes
# SRR2584863_1       SRR2584863      E. coli B REL606 WGS        ~576 MB         no auto-chunk
# ERR552797_1        ERR552797       M. tuberculosis H37Rv WGS   ~217 MB         no auto-chunk
# SRR554369_1        SRR554369       P. aeruginosa PAO1 WGS      ~334 MB         no auto-chunk  [SPRING gold 0.2416 bpb]
# ERR5181310_1       ERR5181310      SARS-CoV-2 amplicon WGS     ~30 MB          no auto-chunk
# ERR17740259_1      ERR17740259     S. aureus WGS ~148bp        ~970 MB         no auto-chunk  (Firmicutes, 33% GC)
# SRR065390_1        SRR065390       C. elegans N2 WGS 100bp PE  ~11 GB          ARCS_AUTOCHUNK_MB=25000 required  [SPRING+PgRC2 benchmark]
# DRR976266_1        DRR976266       S. cerevisiae WGS 150bp     ~1.67 GB        no auto-chunk
# SRR870667_1        SRR870667       T. cacao WGS 69bp SE        ~15 GB          ARCS_AUTOCHUNK_MB=25000 required  [SPRING gold 1.2621 bpb]
# SRR36741279_1      SRR36741279     Leishmania major WGS ~75bp  ~1.70 GB        no auto-chunk
# SRR37283774_1      SRR37283774     P. falciparum WGS ~100bp    ~669 MB         no auto-chunk
#
# 5 human chr20 subset datasets — from GIAB S3 (no account needed)
# HG001_pooled.fq   HG001 NA12878   European CEU chr20           ~37 MB          no auto-chunk
# HG002_pooled.fq   HG002 NA24385   Ashkenazi son chr20          ~37 MB          no auto-chunk
# HG003_pooled.fq   HG003 NA24149   Ashkenazi father chr20       ~37 MB          no auto-chunk
# HG004_pooled.fq   HG004 NA24143   Ashkenazi mother chr20       ~37 MB          no auto-chunk
# HG005_pooled.fq   HG005 NA24631   Han Chinese son chr20        ~37 MB          no auto-chunk
# Human full WGS (~150 GB) exceeds single-pass assembly; chr20 subsets shown as supplementary.

SRA_IDS=(
    SRR2584863
    ERR552797
    SRR554369
    ERR5181310
    ERR17740259
    SRR065390
    DRR976266
    SRR870667
    SRR36741279
    SRR37283774
)

# ── Phase 1: Tool checks ─────────────────────────────────────────────────────
phase "1.1" "Check required tools"
for t in prefetch fasterq-dump aws; do
    if command -v "$t" &>/dev/null; then
        pinfo "$t found: $(command -v "$t")"
    else
        pfail "$t not found — install SRA Toolkit + AWS CLI"
    fi
done
pdone "all tools present"

# ── Phase 2: SRA downloads ────────────────────────────────────────────────────
phase "2" "Download SRA datasets (${#SRA_IDS[@]} accessions)"
for ACC in "${SRA_IDS[@]}"; do
    OUT1="$DATA_DIR/${ACC}_1.fq"
    OUT2="$DATA_DIR/${ACC}_2.fq"
    phase "2.${ACC}" "Download $ACC"
    if [ -s "$OUT1" ]; then
        pskip "$ACC already downloaded: $OUT1"
        continue
    fi
    pinfo "prefetch $ACC → $DATA_DIR/prefetch/$ACC"
    mkdir -p "$DATA_DIR/prefetch"
    prefetch --output-directory "$DATA_DIR/prefetch" "$ACC" 2>&1 | tail -5 || pfail "$ACC prefetch failed"

    pinfo "fasterq-dump $ACC → $DATA_DIR"
    fasterq-dump \
        --outdir "$DATA_DIR" \
        --threads "$JOBS" \
        --progress \
        --skip-technical \
        "$DATA_DIR/prefetch/$ACC/$ACC.sra" 2>&1 | tail -5 || pfail "$ACC fasterq-dump failed"

    # fasterq-dump outputs .fastq; rename to .fq
    for suffix in "" "_1" "_2"; do
        src="$DATA_DIR/${ACC}${suffix}.fastq"
        dst="$DATA_DIR/${ACC}${suffix}.fq"
        if [ -s "$src" ] && [ ! -s "$dst" ]; then
            mv "$src" "$dst"
            pinfo "renamed: ${ACC}${suffix}.fastq → ${ACC}${suffix}.fq"
        fi
    done

    # For SE datasets fasterq-dump emits ${ACC}.fastq (no _1 suffix); normalise to _1
    SE_SRC="$DATA_DIR/${ACC}.fq"
    if [ ! -s "$OUT1" ] && [ -s "$SE_SRC" ]; then
        mv "$SE_SRC" "$OUT1"
        pinfo "SE read: renamed ${ACC}.fq → ${ACC}_1.fq"
    fi

    # Verify _1 exists
    if [ -s "$OUT1" ]; then
        verify_fq "$OUT1" "$ACC"
        pdone "$ACC download complete"
    else
        pfail "$ACC: expected $OUT1 but not found after fasterq-dump"
    fi
done

# ── Phase 3: GIAB HG001-HG005 chr20 from S3 ─────────────────────────────────
phase "3" "Download GIAB HG001-HG005 chr20 reads from S3 (5 human subset datasets)"
pinfo "Using s3://giab — public bucket, no credentials needed"
pinfo "These 5 chr20 subsets are used for BOTH Claim 1 compression AND Claim 2 variant calling"

# HG001-HG005 chr20 subsets from GIAB S3
# If any path 404s, list bucket to find correct path:
#   aws s3 ls --no-sign-request s3://giab/data/AshkenazimTrio/HG002_NA24385_son/ --recursive | grep chr20

declare -A GIAB_S3
GIAB_S3[HG001]="s3://giab/data/NA12878_HG001/NIST_HiSeq_HG001_Homogeneity-10953946/HG001_HiSeq300x_subsetN_chr20.fastq.gz"
GIAB_S3[HG002]="s3://giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_subsetN_chr20.fastq.gz"
GIAB_S3[HG003]="s3://giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-10953946/HG003_HiSeq300x_subsetN_chr20.fastq.gz"
GIAB_S3[HG004]="s3://giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-10953946/HG004_HiSeq300x_subsetN_chr20.fastq.gz"
GIAB_S3[HG005]="s3://giab/data/ChineseTrio/HG005_NA24631_son/NIST_HiSeq_HG005_Homogeneity-10953946/HG005_HiSeq300x_subsetN_chr20.fastq.gz"

for IND in HG001 HG002 HG003 HG004 HG005; do
    phase "3.${IND}" "Download $IND chr20 FASTQ"
    OUT="$DATA_DIR/${IND}_pooled.fq"
    if [ -s "$OUT" ]; then
        pskip "$IND already downloaded: $OUT"
        continue
    fi
    S3PATH="${GIAB_S3[$IND]}"
    pinfo "S3: $S3PATH"
    # Try direct gz download; fall back to listing if path wrong
    TMPGZ="$DATA_DIR/${IND}_chr20.fastq.gz"
    if aws s3 cp --no-sign-request "$S3PATH" "$TMPGZ" 2>/dev/null; then
        pinfo "decompressing $TMPGZ"
        if command -v pigz &>/dev/null; then
            pigz -d -c "$TMPGZ" > "$OUT"
        else
            gzip -d -c "$TMPGZ" > "$OUT"
        fi
        rm -f "$TMPGZ"
        verify_fq "$OUT" "$IND"
        pdone "$IND chr20 ready: $OUT"
    else
        pinfo "Direct path failed — listing S3 to find correct path"
        pinfo "Run: aws s3 ls --no-sign-request s3://giab/data/ --recursive | grep -i 'chr20.*fastq'"
        pfail "$IND: S3 download failed — check path in download.sh GIAB_S3 table"
    fi
done

# ── Phase 4: Verify all expected files ───────────────────────────────────────
phase "4" "Final verification of all dataset files"
EXPECTED=(
    "$DATA_DIR/SRR2584863_1.fq"
    "$DATA_DIR/ERR552797_1.fq"
    "$DATA_DIR/SRR554369_1.fq"
    "$DATA_DIR/ERR5181310_1.fq"
    "$DATA_DIR/ERR17740259_1.fq"
    "$DATA_DIR/SRR065390_1.fq"
    "$DATA_DIR/DRR976266_1.fq"
    "$DATA_DIR/SRR870667_1.fq"
    "$DATA_DIR/SRR36741279_1.fq"
    "$DATA_DIR/SRR37283774_1.fq"
    "$DATA_DIR/HG001_pooled.fq"
    "$DATA_DIR/HG002_pooled.fq"
    "$DATA_DIR/HG003_pooled.fq"
    "$DATA_DIR/HG004_pooled.fq"
    "$DATA_DIR/HG005_pooled.fq"
)
MISSING=0
for F in "${EXPECTED[@]}"; do
    if [ -s "$F" ]; then
        verify_fq "$F" "$(basename "$F")"
    else
        pinfo "MISSING: $F"
        MISSING=$(( MISSING + 1 ))
    fi
done

if [ "$MISSING" -gt 0 ]; then
    pfail "$MISSING file(s) missing — fix above errors before running benchmark"
fi

pdone "All $(( ${#EXPECTED[@]} - MISSING )) files verified"

# ── Phase 5: Claim 2 prerequisites — chr20.fa + GIAB truth VCFs ──────────────
phase "5" "Download Claim 2 prerequisites (chr20 reference + GIAB truth VCFs)"
REFS_DIR="${HOME}/refs"
TRUTH_DIR="${HOME}/giab_truth"
mkdir -p "$REFS_DIR" "$TRUTH_DIR"

# chr20 reference — GRCh37 (hg19), chromosome named "20" (no chr prefix)
phase "5.1" "chr20.fa reference (GRCh37 / hg19)"
CHR20_FA="$REFS_DIR/chr20.fa"
if [ -s "$CHR20_FA" ]; then
    pskip "chr20.fa already present ($(stat -c %s "$CHR20_FA") bytes)"
else
    pinfo "Downloading chr20 from UCSC hg19 (~65 MB compressed)"
    TMPGZ="$REFS_DIR/chr20.fa.gz"
    curl -L -o "$TMPGZ" \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz" \
        || pfail "chr20.fa.gz download failed"
    pinfo "Decompressing and stripping 'chr' prefix from sequence name"
    zcat "$TMPGZ" | sed 's/^>chr/>/' > "$CHR20_FA"
    rm -f "$TMPGZ"
    [ -s "$CHR20_FA" ] && pdone "chr20.fa ready: $CHR20_FA" || pfail "chr20.fa empty after download"
fi

# GIAB truth VCFs v4.2.1 GRCh37 for HG002-HG005
# HG001 is unused (Claim 2 benchmarks HG002-HG005 only)
phase "5.2" "GIAB truth VCFs v4.2.1 GRCh37 (HG002-HG005, ~2.8 GB total)"
_GIAB_FTP="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release"

declare -A VCF_URL BED_URL
VCF_URL[HG002]="$_GIAB_FTP/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED_URL[HG002]="$_GIAB_FTP/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
VCF_URL[HG003]="$_GIAB_FTP/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED_URL[HG003]="$_GIAB_FTP/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
VCF_URL[HG004]="$_GIAB_FTP/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED_URL[HG004]="$_GIAB_FTP/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
VCF_URL[HG005]="$_GIAB_FTP/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh37/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BED_URL[HG005]="$_GIAB_FTP/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh37/HG005_GRCh37_1_22_v4.2.1_benchmark.bed"

for IND in HG002 HG003 HG004 HG005; do
    VCF="$TRUTH_DIR/${IND}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
    TBI="$VCF.tbi"
    BED_FILE="$TRUTH_DIR/$(basename "${BED_URL[$IND]}")"
    phase "5.2.${IND}" "$IND truth VCF + BED"
    if [ -s "$VCF" ] && [ -s "$TBI" ]; then
        pskip "$IND VCF already present"
    else
        pinfo "Downloading $IND truth VCF (~700 MB)"
        curl -L -o "$VCF" "${VCF_URL[$IND]}" || pfail "$IND VCF download failed"
        tabix -p vcf "$VCF" || pfail "$IND VCF tabix failed"
        pdone "$IND VCF + index ready"
    fi
    if [ -s "$BED_FILE" ]; then
        pskip "$IND BED already present"
    else
        pinfo "Downloading $IND confident BED"
        curl -L -o "$BED_FILE" "${BED_URL[$IND]}" || pfail "$IND BED download failed"
        pdone "$IND BED ready"
    fi
done

pdone "All Claim 2 prerequisites ready"
echo ""
echo "========================================"
echo "  DOWNLOAD COMPLETE"
echo "  DATA_DIR  : $DATA_DIR"
echo "  REFS_DIR  : $REFS_DIR  (chr20.fa)"
echo "  TRUTH_DIR : $TRUTH_DIR (GIAB HG002-HG005 VCFs)"
echo "  Run benchmark with:"
echo "    bash benchmark/benchmark.sh claim1 $DATA_DIR <ARCS_BIN> ./results"
echo "========================================"
