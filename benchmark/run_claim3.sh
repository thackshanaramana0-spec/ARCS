#!/bin/bash
# Claim 3 — ADDRESSABLE: T6 archive-derived analysis vs conventional pipeline
#
# Usage:
#   bash benchmark/run_claim3.sh <DATA_DIR> <ARCS_BIN> [OUT_DIR]
#
# DATA_DIR : directory containing 5M-read subsets:
#            ecoli.fq  celegans_5M.fq  mtb.fq  HG002_pooled.fq  drosophila_5M.fq
# ARCS_BIN : path to arcs binary
# OUT_DIR  : results directory (default: ./results/claim3)
#
# Required env vars (with defaults):
#   REF_DIR    directory containing organism references for BWA  [~/refs]
#   SPADES     path to spades.py                                 [spades.py]
#
# References needed in REF_DIR:
#   ecoli.fa        E. coli K-12 MG1655
#   celegans.fa     C. elegans WBcel235
#   mtb.fa          M. tuberculosis H37Rv
#   hg38_chr20.fa   H. sapiens chr20
#   drosophila.fa   D. melanogaster dm6
#
# 3 operations per dataset:
#   export   : arcs export  vs  SPAdes de-novo assembly
#   coverage : arcs coverage vs  BWA + mosdepth
#   query    : arcs query (ARCS only, 0-100000 region)

set -euo pipefail

# ── Helpers ───────────────────────────────────────────────────────────────────
_PH=""
log()   { echo "[$(date +%H:%M:%S)] $*"; }
phase() { _PH="$1"; echo ""; echo "[Phase $_PH] $2"; }
pdone() { echo "[Phase $_PH] DONE — $*"; }
pfail() { echo "[Phase $_PH] FAIL — $*"; exit 1; }
pskip() { echo "[Phase $_PH] SKIP — $*"; }
pinfo() { echo "[Phase $_PH]   → $*"; }
step()  { echo ""; echo "══════════════════════════════════════════════"; echo "[STEP] $*"; echo "══════════════════════════════════════════════"; }
ok()    { echo "  [OK]   $*"; }
warn()  { echo "  [WARN] $*"; }
fail()  { echo "  [FAIL] $*"; exit 1; }
chk()   { echo "  [CHK]  $*"; }

parse_time_v() {
    local logf="$1"
    local wall vmhwm
    wall=$(grep "Elapsed (wall clock)" "$logf" | awk '{
        n=split($NF,a,":");
        if(n==3) printf "%.2f", a[1]*3600+a[2]*60+a[3];
        else if(n==2) printf "%.2f", a[1]*60+a[2];
        else printf "%.2f", a[1];
    }')
    vmhwm=$(grep "Maximum resident set size" "$logf" | awk '{print $NF}')
    echo "${wall:-0} ${vmhwm:-0}"
}

# ── Args ──────────────────────────────────────────────────────────────────────
DATA_DIR="${1:?Usage: $0 <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
ARCS_BIN="${2:?Usage: $0 <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
OUT_DIR="${3:-./results/claim3}"
REF_DIR="${REF_DIR:-$HOME/refs}"
SPADES="${SPADES:-spades.py}"

mkdir -p "$OUT_DIR"
WD="$OUT_DIR/workdir"
mkdir -p "$WD"

NPROC=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)

log "Claim 3 — ADDRESSABLE (T6)"
log "DATA_DIR : $DATA_DIR"
log "ARCS_BIN : $ARCS_BIN"
log "OUT_DIR  : $OUT_DIR"
log "REF_DIR  : $REF_DIR"
log "SPADES   : $SPADES"

# ── Dataset config ────────────────────────────────────────────────────────────
# name → (filename, reference, query_region)
declare -A DS_FQ DS_REF DS_QUERY
DS_FQ[ecoli]="ecoli.fq"
DS_REF[ecoli]="$REF_DIR/ecoli.fa"
DS_QUERY[ecoli]="0-100000"

DS_FQ[celegans]="celegans_5M.fq"
DS_REF[celegans]="$REF_DIR/celegans.fa"
DS_QUERY[celegans]="0-100000"

DS_FQ[mtb]="mtb.fq"
DS_REF[mtb]="$REF_DIR/mtb.fa"
DS_QUERY[mtb]="0-100000"

DS_FQ[HG002]="HG002_pooled.fq"
DS_REF[HG002]="$REF_DIR/hg38_chr20.fa"
DS_QUERY[HG002]="2000000-2100000"

DS_FQ[drosophila]="drosophila_5M.fq"
DS_REF[drosophila]="$REF_DIR/drosophila.fa"
DS_QUERY[drosophila]="0-100000"

DS_ORDER="ecoli celegans mtb HG002 drosophila"

RESULTS_CSV="$OUT_DIR/t6_results.csv"
echo "Dataset,Operation,ARCS_time_s,ARCS_ram_kb,Conv_tool,Conv_time_s,Conv_ram_kb,Speedup" > "$RESULTS_CSV"

# ── STEP 0: Prerequisites ─────────────────────────────────────────────────────
step "0/8 Prerequisites check"

phase "0.1" "Check arcs binary"
[ -x "$ARCS_BIN" ] && pdone "arcs found: $ARCS_BIN" || pfail "arcs binary not found: $ARCS_BIN"

phase "0.2" "Check SPAdes"
if command -v "$SPADES" &>/dev/null; then
    pdone "spades found: $(command -v $SPADES)"
elif command -v spades.py &>/dev/null; then
    SPADES=spades.py; pdone "spades.py found on PATH"
else
    pfail "SPAdes not found — install: conda install -c bioconda spades"
fi

phase "0.3" "Check bwa + samtools + mosdepth"
for t in bwa samtools mosdepth; do
    command -v "$t" &>/dev/null && pinfo "$t found: $(command -v $t)" || pfail "$t not found on PATH"
done
pdone "all tools present"

# ── STEP 1: Check input files ─────────────────────────────────────────────────
step "1/8 Check input FASTQ files"
_DSI=0
for DS in $DS_ORDER; do
    _DSI=$((_DSI+1))
    phase "1.${_DSI}" "Check $DS FASTQ"
    FQ="$DATA_DIR/${DS_FQ[$DS]}"
    if [ -s "$FQ" ]; then
        N=$(( $(wc -l < "$FQ") / 4 ))
        SZ=$(stat -c %s "$FQ")
        pinfo "$DS — $N reads, $SZ bytes"
        if [ "$DS" = "drosophila" ]; then
            LAST=$(tail -1 "$FQ")
            if [ -z "$LAST" ]; then
                warn "drosophila last line is empty — possible truncation! Verify file before continuing"
            else
                pinfo "drosophila integrity check: last line non-empty — OK"
            fi
        fi
        pdone "$DS file OK"
    else
        pfail "$DS FASTQ not found: $FQ"
    fi
done

# ── STEP 2: Check reference files and build BWA indices ──────────────────────
step "2/8 Check reference files + BWA indices"
_DSI=0
for DS in $DS_ORDER; do
    _DSI=$((_DSI+1))
    phase "2.${_DSI}" "Check/build BWA index for $DS"
    REF="${DS_REF[$DS]}"
    if [ -s "$REF" ]; then
        pinfo "$DS ref found ($(stat -c %s "$REF") bytes)"
        if [ ! -s "${REF}.bwt" ]; then
            pinfo "$DS BWA index missing — building now (1-5 min)"
            bwa index "$REF" 2>/dev/null
            pdone "$DS BWA index built"
        else
            pdone "$DS BWA index exists"
        fi
    else
        pskip "$DS reference not found at $REF — coverage step will be skipped"
    fi
done

# ── STEP 3: Compress all datasets into ARCS archives ─────────────────────────
step "3/8 Build ARCS archives (compress all datasets)"
export ARCS_PAR_SHARDS=$NPROC
_DSI=0
for DS in $DS_ORDER; do
    _DSI=$((_DSI+1))
    phase "3.${_DSI}" "ARCS compress $DS"
    FQ="$DATA_DIR/${DS_FQ[$DS]}"
    ARC="$WD/${DS}.arcs"
    if [ ! -s "$ARC" ]; then
        pinfo "compressing $DS → $ARC"
        TF="$WD/${DS}_compress_time.txt"
        /usr/bin/time -v "$ARCS_BIN" compress "$FQ" "$ARC" 2>"$TF"
        read -r CT CRAM <<< "$(parse_time_v "$TF")"
        pdone "$DS compressed — ${CT}s, ${CRAM}KB RAM, archive=$(stat -c %s "$ARC") bytes"
    else
        pskip "$DS archive already exists: $ARC"
    fi
done

# ── STEP 4: T6 — Export (ARCS vs SPAdes) ─────────────────────────────────────
step "4/8 T6 Export — arcs export vs SPAdes de-novo"
_DSI=0
for DS in $DS_ORDER; do
    _DSI=$((_DSI+1))
    ARC="$WD/${DS}.arcs"

    # ARCS export
    phase "4.${_DSI}.arcs" "$DS — arcs export (pseudogenome → FASTA)"
    EXPORT_OUT="$WD/${DS}_export.fa"
    if [ ! -s "$EXPORT_OUT" ]; then
        pinfo "running: arcs export $ARC $EXPORT_OUT"
        TF="$WD/${DS}_export_time.txt"
        /usr/bin/time -v "$ARCS_BIN" export "$ARC" "$EXPORT_OUT" 2>"$TF"
        read -r ARCS_T ARCS_RAM <<< "$(parse_time_v "$TF")"
        N_CONTIGS=$(grep -c '^>' "$EXPORT_OUT" 2>/dev/null || echo 0)
        pdone "$DS arcs export — ${ARCS_T}s, ${ARCS_RAM}KB RAM, $N_CONTIGS contigs"
    else
        ARCS_T="cached"; ARCS_RAM="cached"
        pskip "$DS arcs export already done"
    fi

    # SPAdes de-novo
    phase "4.${_DSI}.spades" "$DS — SPAdes de-novo assembly"
    FQ="$DATA_DIR/${DS_FQ[$DS]}"
    SPADES_OUT="$WD/${DS}_spades"
    if [ ! -d "$SPADES_OUT" ]; then
        pinfo "running: spades.py --se $FQ -o $SPADES_OUT -t $NPROC --careful (slow — may take 5-30 min)"
        TF="$WD/${DS}_spades_time.txt"
        /usr/bin/time -v "$SPADES" --se "$FQ" -o "$SPADES_OUT" -t "$NPROC" --careful \
            >/dev/null 2>"$TF" || { warn "[$DS] SPAdes failed — check $WD/${DS}_spades_time.txt"; }
        read -r SPADES_T SPADES_RAM <<< "$(parse_time_v "$TF")"
        pdone "$DS SPAdes done — ${SPADES_T}s, ${SPADES_RAM}KB RAM"
    else
        SPADES_T="cached"; SPADES_RAM="cached"
        pskip "$DS SPAdes output already exists"
    fi

    # Speedup + CSV
    phase "4.${_DSI}.speedup" "$DS — export speedup calculation"
    if [[ "$ARCS_T" != "cached" && "$SPADES_T" != "cached" ]]; then
        SPEEDUP=$(awk "BEGIN{printf \"%.1f\", $SPADES_T / ($ARCS_T + 0.001)}")
        pdone "$DS export: ARCS ${ARCS_T}s vs SPAdes ${SPADES_T}s = ${SPEEDUP}× speedup"
        echo "$DS,export,$ARCS_T,$ARCS_RAM,SPAdes,$SPADES_T,$SPADES_RAM,$SPEEDUP" >> "$RESULTS_CSV"
    else
        pskip "$DS export: at least one side cached — speedup not computed"
        echo "$DS,export,$ARCS_T,$ARCS_RAM,SPAdes,$SPADES_T,$SPADES_RAM,N/A" >> "$RESULTS_CSV"
    fi
done

# ── STEP 5: T6 — Coverage (ARCS vs BWA+mosdepth) ────────────────────────────
step "5/8 T6 Coverage — arcs coverage vs BWA+mosdepth"
_DSI=0
for DS in $DS_ORDER; do
    _DSI=$((_DSI+1))
    ARC="$WD/${DS}.arcs"
    REF="${DS_REF[$DS]}"
    FQ="$DATA_DIR/${DS_FQ[$DS]}"

    # ARCS coverage
    phase "5.${_DSI}.arcs" "$DS — arcs coverage"
    COV_OUT="$WD/${DS}_coverage.tsv"
    if [ ! -s "$COV_OUT" ]; then
        pinfo "running: arcs coverage $ARC $COV_OUT"
        TF="$WD/${DS}_cov_time.txt"
        /usr/bin/time -v "$ARCS_BIN" coverage "$ARC" "$COV_OUT" 2>"$TF"
        read -r ARCS_T ARCS_RAM <<< "$(parse_time_v "$TF")"
        pdone "$DS arcs coverage — ${ARCS_T}s, ${ARCS_RAM}KB RAM"
    else
        ARCS_T="cached"; ARCS_RAM="cached"
        pskip "$DS arcs coverage already done"
    fi

    # BWA + mosdepth conventional
    phase "5.${_DSI}.bwa" "$DS — BWA mem + samtools sort + mosdepth"
    if [ -s "$REF" ]; then
        BAM="$WD/${DS}_bwa.bam"
        if [ ! -s "$BAM" ]; then
            pinfo "running: bwa mem -t $NPROC $REF $FQ | samtools sort → $BAM"
            TF_BWA="$WD/${DS}_bwa_time.txt"
            /usr/bin/time -v bash -c \
                "bwa mem -t $NPROC '$REF' '$FQ' 2>/dev/null | samtools sort -@ $NPROC -o '$BAM' && samtools index '$BAM'" \
                2>"$TF_BWA"
            read -r BWA_T BWA_RAM <<< "$(parse_time_v "$TF_BWA")"
            pinfo "$DS BWA align+sort done — ${BWA_T}s, ${BWA_RAM}KB RAM"
        else
            BWA_T=0; BWA_RAM=0
            pskip "$DS BAM already exists"
        fi
        MOSD_OUT="$WD/${DS}_mosd"
        if [ ! -d "$MOSD_OUT" ]; then
            pinfo "running: mosdepth --threads $NPROC $MOSD_OUT $BAM"
            TF_MSD="$WD/${DS}_mosd_time.txt"
            /usr/bin/time -v mosdepth --threads "$NPROC" "$MOSD_OUT" "$BAM" 2>"$TF_MSD"
            read -r MSD_T MSD_RAM <<< "$(parse_time_v "$TF_MSD")"
            CONV_T=$(awk "BEGIN{printf \"%.2f\", ${BWA_T:-0} + ${MSD_T:-0}}")
            CONV_RAM=$(awk "BEGIN{printf \"%d\", (${BWA_RAM:-0} > ${MSD_RAM:-0}) ? ${BWA_RAM:-0} : ${MSD_RAM:-0}}")
            pdone "$DS mosdepth done — total conv ${CONV_T}s, peak ${CONV_RAM}KB RAM"
        else
            CONV_T="cached"; CONV_RAM="cached"
            pskip "$DS mosdepth output already exists"
        fi

        phase "5.${_DSI}.speedup" "$DS — coverage speedup calculation"
        if [[ "$ARCS_T" != "cached" && "$CONV_T" != "cached" ]]; then
            SPEEDUP=$(awk "BEGIN{printf \"%.1f\", $CONV_T / ($ARCS_T + 0.001)}")
            pdone "$DS coverage: ARCS ${ARCS_T}s vs BWA+mosdepth ${CONV_T}s = ${SPEEDUP}× speedup"
            echo "$DS,coverage,$ARCS_T,$ARCS_RAM,BWA+mosdepth,$CONV_T,$CONV_RAM,$SPEEDUP" >> "$RESULTS_CSV"
        else
            pskip "$DS coverage: at least one side cached — speedup not computed"
            echo "$DS,coverage,$ARCS_T,$ARCS_RAM,BWA+mosdepth,$CONV_T,$CONV_RAM,N/A" >> "$RESULTS_CSV"
        fi
    else
        pskip "$DS reference not found at $REF — BWA+mosdepth skipped"
        echo "$DS,coverage,$ARCS_T,$ARCS_RAM,BWA+mosdepth,N/A,N/A,N/A" >> "$RESULTS_CSV"
    fi
done

# ── STEP 6: T6 — Query (ARCS only) ───────────────────────────────────────────
step "6/8 T6 Query — arcs query (ARCS only)"
_DSI=0
for DS in $DS_ORDER; do
    _DSI=$((_DSI+1))
    phase "6.${_DSI}" "$DS — arcs query region ${DS_QUERY[$DS]}"
    ARC="$WD/${DS}.arcs"
    REGION="${DS_QUERY[$DS]}"
    QUERY_OUT="$WD/${DS}_query.fq"
    pinfo "running: arcs query $ARC $REGION $QUERY_OUT"
    TF="$WD/${DS}_query_time.txt"
    /usr/bin/time -v "$ARCS_BIN" query "$ARC" "$REGION" "$QUERY_OUT" 2>"$TF"
    read -r QT QRAM <<< "$(parse_time_v "$TF")"
    N_READS=$(( $(wc -l < "$QUERY_OUT") / 4 ))
    pdone "$DS query — ${QT}s, ${QRAM}KB RAM, $N_READS reads returned for region $REGION"
    echo "$DS,query,$QT,$QRAM,none,N/A,N/A,N/A" >> "$RESULTS_CSV"
done

# ── STEP 7: Summary table ─────────────────────────────────────────────────────
phase "7" "Results summary"
pinfo "T6 results written to: $RESULTS_CSV"
echo ""
printf "%-12s %-10s %-10s %-10s %-14s %-10s %-10s %-8s\n" \
    "DATASET" "OPERATION" "ARCS_s" "ARCS_MB" "CONV_TOOL" "CONV_s" "CONV_MB" "SPEEDUP"
echo "─────────────────────────────────────────────────────────────────────────────────────"
tail -n +2 "$RESULTS_CSV" | while IFS=, read DS OP AT AR CT CVT CVR SPD; do
    ARMB=$(awk "BEGIN{printf \"%.0f\", ${AR:-0}/1024}" 2>/dev/null || echo "?")
    CVMB=$(awk "BEGIN{printf \"%.0f\", ${CVR:-0}/1024}" 2>/dev/null || echo "?")
    printf "%-12s %-10s %-10s %-10s %-14s %-10s %-10s %-8s\n" \
        "$DS" "$OP" "$AT" "${ARMB}MB" "$CT" "$CVT" "${CVMB}MB" "$SPD"
done
pdone "summary table printed"

phase "8" "DONE"
pdone "All Claim 3 results in: $OUT_DIR"
pinfo "CSV: $RESULTS_CSV"
pinfo "Finished at: $(date)"
