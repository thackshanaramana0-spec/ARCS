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
# Name               Accession       Organism                    Approx _1 size  Notes
# SRR2584863_1       SRR2584863      E. coli B REL606 WGS        ~576 MB         no auto-chunk
# ERR552797_1        ERR552797       M. tuberculosis H37Rv WGS   ~217 MB         no auto-chunk
# SRR554369_1        SRR554369       P. aeruginosa PAO1 WGS      ~334 MB         no auto-chunk
# ERR5181310_1       ERR5181310      SARS-CoV-2 amplicon WGS     ~30 MB          no auto-chunk
# SRR16357346_1      SRR16357346     C. elegans N2 WGS ~135bp    ~1.42 GB        no auto-chunk
# DRR976266_1        DRR976266       S. cerevisiae WGS 150bp     ~1.67 GB        no auto-chunk
# SRR1945765_1       SRR1945765      Arabidopsis thaliana WGS    ~1.95 GB        no auto-chunk
# SRR36741279_1      SRR36741279     Leishmania major WGS ~75bp  ~1.70 GB        no auto-chunk
# SRR37283774_1      SRR37283774     P. falciparum WGS ~100bp    ~669 MB         no auto-chunk
#
# GIAB HG002-HG005 chr20 reads — from S3 (no account needed)
# Compression benchmark uses HG002_chr20.fq; variant calling uses HG002-HG005 individually.
#
# All tools run at ARCS defaults (no special flags). Auto-chunk is ARCS default
# behaviour for files >2 GB; no ARCS_CHUNK_THREADS override needed.

SRA_IDS=(
    SRR2584863
    ERR552797
    SRR554369
    ERR5181310
    SRR16357346
    DRR976266
    SRR1945765
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

# ── Phase 3: GIAB HG002-HG005 chr20 from S3 ─────────────────────────────────
phase "3" "Download GIAB HG002-HG005 chr20 reads from S3"
pinfo "Using s3://giab — public bucket, no credentials needed"

# HG002-HG005 novaseq 2×150 chr20 subsets
# These are pre-subset chr20 BAMs from GIAB; we use the FASTQs where available
# Strategy: download from s3://giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/
#           filter chr20 reads from sorted BAM using samtools view

declare -A GIAB_S3 GIAB_OUT
GIAB_S3[HG002]="s3://giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_subsetN_chr20.fastq.gz"
GIAB_S3[HG003]="s3://giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-10953946/HG003_HiSeq300x_subsetN_chr20.fastq.gz"
GIAB_S3[HG004]="s3://giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-10953946/HG004_HiSeq300x_subsetN_chr20.fastq.gz"
GIAB_S3[HG005]="s3://giab/data/ChineseTrio/HG005_NA24631_son/NIST_HiSeq_HG005_Homogeneity-10953946/HG005_HiSeq300x_subsetN_chr20.fastq.gz"

# Note: exact S3 paths may differ — if the above 404, list with:
#   aws s3 ls --no-sign-request s3://giab/data/AshkenazimTrio/HG002_NA24385_son/ --recursive | grep chr20
# Then update GIAB_S3 paths accordingly.

for IND in HG002 HG003 HG004 HG005; do
    phase "3.${IND}" "Download $IND chr20 FASTQ"
    OUT="$DATA_DIR/${IND}_chr20.fq"
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
    "$DATA_DIR/SRR16357346_1.fq"
    "$DATA_DIR/DRR976266_1.fq"
    "$DATA_DIR/SRR1945765_1.fq"
    "$DATA_DIR/SRR36741279_1.fq"
    "$DATA_DIR/SRR37283774_1.fq"
    "$DATA_DIR/HG002_chr20.fq"
    "$DATA_DIR/HG003_chr20.fq"
    "$DATA_DIR/HG004_chr20.fq"
    "$DATA_DIR/HG005_chr20.fq"
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
echo ""
echo "========================================"
echo "  DOWNLOAD COMPLETE"
echo "  DATA_DIR: $DATA_DIR"
echo "  Run benchmark with:"
echo "    bash benchmark/benchmark.sh claim1 $DATA_DIR <ARCS_BIN> ./results"
echo "========================================"
