#!/bin/bash
# Block 1 — Claim 1 (COMPACT): T1 ratio + T2 speed + RAM
# Runs ARCS, SPRING, Genozip, PgRC2 on all 10 datasets.
#
# Usage:
#   bash benchmark/run_block1.sh <DATA_DIR> <ARCS_BIN> [OUT_DIR]
#
# DATA_DIR : directory with FASTQ files (exact names below)
# ARCS_BIN : path to arcs binary
# OUT_DIR  : results directory (default: ./results/block1)
#
# Tools expected on PATH: spring, genozip, genounzip, pgrc2 (optional)
# Timing: single run per tool, /usr/bin/time -v for wall+VmHWM
# NEVER run two timed jobs concurrently — contaminates measurements

set -euo pipefail

DATA_DIR="${1:?Usage: $0 <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
ARCS_BIN="${2:?Usage: $0 <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
OUT_DIR="${3:-./results/block1}"
mkdir -p "$OUT_DIR"

WD=/tmp/arcs_bench_wd
mkdir -p "$WD"

NPROC=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)

# ── Progress helpers ──────────────────────────────────────────────────────────
_PHASE=""
phase()  { _PHASE="$1"; echo ""; echo "[Phase $_PHASE] $2"; }
pdone()  { echo "[Phase $_PHASE] DONE — $*"; }
pfail()  { echo "[Phase $_PHASE] FAIL — $*"; exit 1; }
pskip()  { echo "[Phase $_PHASE] SKIP — $*"; }
pinfo()  { echo "[Phase $_PHASE]   → $*"; }

# ── Phase 1: Dataset discovery ───────────────────────────────────────────────
phase "1" "Dataset discovery in $DATA_DIR"
mapfile -t FQ_FILES < <(find "$DATA_DIR" -maxdepth 1 \( -name "*.fq" -o -name "*.fastq" \) | sort)
if [ ${#FQ_FILES[@]} -eq 0 ]; then
    pfail "no .fq or .fastq files found in $DATA_DIR"
fi
for f in "${FQ_FILES[@]}"; do pinfo "found: $f"; done
pdone "${#FQ_FILES[@]} dataset(s) found"

# ── Parse /usr/bin/time -v output → wall seconds and VmHWM KB ───────────────
parse_time_v() {
    local logf="$1"
    local wall vmhwm
    # format: h:mm:ss.cc (3 fields) or m:ss.cc (2 fields)
    wall=$(grep "Elapsed (wall clock)" "$logf" | awk '{
        n=split($NF,a,":");
        if(n==3) printf "%.2f", a[1]*3600 + a[2]*60 + a[3];
        else if(n==2) printf "%.2f", a[1]*60 + a[2];
        else printf "%.2f", a[1];
    }')
    # VmHWM in KB (Maximum resident set size (kbytes))
    vmhwm=$(grep "Maximum resident set size" "$logf" | awk '{print $NF}')
    echo "${wall:-0} ${vmhwm:-0}"
}

# ── Lossless check (order-free, CRLF-safe) ───────────────────────────────────
losscmp() {
    paste - - - - < "$1" | tr -d '\r' | sort > /tmp/_orig_$$
    paste - - - - < "$2" | tr -d '\r' | sort > /tmp/_dec_$$
    local res=LOSSY
    if cmp -s /tmp/_orig_$$ /tmp/_dec_$$; then res=LOSSLESS; fi
    rm -f /tmp/_orig_$$ /tmp/_dec_$$
    echo "$res"
}

# ── Phase 2: Verify dataset files ────────────────────────────────────────────
phase "2" "Verify dataset files (organism check)"
for SRC in "${FQ_FILES[@]}"; do
    DS=$(basename "$SRC" | sed 's/\.\(fq\|fastq\)$//')
    FIRST=$(head -2 "$SRC" | tail -1 | cut -c1-30)
    RAW=$(stat -c %s "$SRC")
    N=$(( $(wc -l < "$SRC") / 4 ))
    pinfo "$DS — ${RAW}B — ${N} reads — seq_prefix: $FIRST"
done
pdone "all files readable"

# ── Phase 3+: Main benchmark loop ────────────────────────────────────────────
DS_IDX=0
for SRC in "${FQ_FILES[@]}"; do
    DS=$(basename "$SRC" | sed 's/\.\(fq\|fastq\)$//')
    DS_IDX=$(( DS_IDX + 1 ))

    phase "${DS_IDX}.0" "Dataset $DS — copy to /tmp"
    IN="$WD/$DS.fq"
    cp "$SRC" "$IN"
    RAW=$(stat -c %s "$IN")
    pdone "raw size = ${RAW}B"

    # ── ARCS ─────────────────────────────────────────────────────────────────
    phase "${DS_IDX}.1" "[$DS] ARCS compress"
    A="$WD/$DS.arcs"
    DEC="$WD/$DS.arcs.dec"
    TF=/tmp/time_arcs_c_$$
    # ARCS_PAR_SHARDS not set — C++ default min(4, hw_concurrency) gives best ratio/speed balance
    /usr/bin/time -v "$ARCS_BIN" compress "$IN" "$A" 2>"$TF"
    read -r CWALL CRAMKB <<< "$(parse_time_v "$TF")"
    ARCH=$(stat -c %s "$A")
    pinfo "compress done — archive=${ARCH}B  time=${CWALL}s  ram=${CRAMKB}KB"

    phase "${DS_IDX}.2" "[$DS] ARCS decompress"
    TF2=/tmp/time_arcs_d_$$
    /usr/bin/time -v "$ARCS_BIN" decompress "$A" "$DEC" 2>"$TF2"
    read -r DWALL DRAMKB <<< "$(parse_time_v "$TF2")"
    pinfo "decompress done — time=${DWALL}s  ram=${DRAMKB}KB"

    phase "${DS_IDX}.3" "[$DS] ARCS lossless check"
    LL=$(losscmp "$IN" "$DEC")
    rm -f "$TF" "$TF2" "$A" "$DEC"
    printf "TOOL=ARCS DS=%s raw=%s archive=%s ctime=%.2f dtime=%.2f cram_kb=%s dram_kb=%s lossless=%s\n" \
        "$DS" "$RAW" "$ARCH" "$CWALL" "$DWALL" "$CRAMKB" "$DRAMKB" "$LL" \
        | tee "$OUT_DIR/${DS}__ARCS.log"
    pdone "ARCS — archive=${ARCH}B  ctime=${CWALL}s  dtime=${DWALL}s  $LL"

    # ── SPRING ───────────────────────────────────────────────────────────────
    phase "${DS_IDX}.4" "[$DS] SPRING compress+decompress"
    if command -v spring &>/dev/null; then
        A="$WD/$DS.spring"
        DEC="$WD/$DS.spring.dec"
        TF=/tmp/time_spring_c_$$
        /usr/bin/time -v spring -c -i "$IN" -o "$A" -t "$NPROC" -g 2>"$TF"
        read -r CWALL CRAMKB <<< "$(parse_time_v "$TF")"
        ARCH=$(stat -c %s "$A")
        pinfo "SPRING compress done — archive=${ARCH}B  time=${CWALL}s  ram=${CRAMKB}KB"
        TF2=/tmp/time_spring_d_$$
        /usr/bin/time -v spring -d -i "$A" -o "$DEC" -t "$NPROC" -g 2>"$TF2"
        read -r DWALL DRAMKB <<< "$(parse_time_v "$TF2")"
        LL=$(losscmp "$IN" "$DEC")
        rm -f "$TF" "$TF2" "$A" "$DEC"
        printf "TOOL=SPRING DS=%s raw=%s archive=%s ctime=%.2f dtime=%.2f cram_kb=%s dram_kb=%s lossless=%s\n" \
            "$DS" "$RAW" "$ARCH" "$CWALL" "$DWALL" "$CRAMKB" "$DRAMKB" "$LL" \
            | tee "$OUT_DIR/${DS}__SPRING.log"
        pdone "SPRING — archive=${ARCH}B  ctime=${CWALL}s  dtime=${DWALL}s  $LL"
    else
        pskip "SPRING not on PATH"
    fi

    # ── Genozip ──────────────────────────────────────────────────────────────
    phase "${DS_IDX}.5" "[$DS] Genozip compress+decompress"
    if command -v genozip &>/dev/null; then
        A="$WD/$DS.genozip"
        DEC="$WD/$DS.gz.dec"
        TF=/tmp/time_gz_c_$$
        /usr/bin/time -v genozip --force -o "$A" "$IN" 2>"$TF"
        read -r CWALL CRAMKB <<< "$(parse_time_v "$TF")"
        ARCH=$(stat -c %s "$A")
        pinfo "Genozip compress done — archive=${ARCH}B  time=${CWALL}s  ram=${CRAMKB}KB"
        TF2=/tmp/time_gz_d_$$
        /usr/bin/time -v genounzip --force -o "$DEC" "$A" 2>"$TF2"
        read -r DWALL DRAMKB <<< "$(parse_time_v "$TF2")"
        LL=$(losscmp "$IN" "$DEC")
        rm -f "$TF" "$TF2" "$A" "$DEC"
        printf "TOOL=Genozip DS=%s raw=%s archive=%s ctime=%.2f dtime=%.2f cram_kb=%s dram_kb=%s lossless=%s\n" \
            "$DS" "$RAW" "$ARCH" "$CWALL" "$DWALL" "$CRAMKB" "$DRAMKB" "$LL" \
            | tee "$OUT_DIR/${DS}__Genozip.log"
        pdone "Genozip — archive=${ARCH}B  ctime=${CWALL}s  dtime=${DWALL}s  $LL"
    else
        pskip "Genozip not on PATH"
    fi

    # ── PgRC2 ────────────────────────────────────────────────────────────────
    phase "${DS_IDX}.6" "[$DS] PgRC2 compress+decompress"
    if command -v pgrc2 &>/dev/null; then
        A="$WD/$DS.pgrc2"
        DEC="$WD/$DS.pgrc2.dec"
        TF=/tmp/time_pgrc2_c_$$
        /usr/bin/time -v pgrc2 -c "$IN" "$A" 2>"$TF"
        read -r CWALL CRAMKB <<< "$(parse_time_v "$TF")"
        ARCH=$(stat -c %s "$A")
        pinfo "PgRC2 compress done — archive=${ARCH}B  time=${CWALL}s  ram=${CRAMKB}KB"
        TF2=/tmp/time_pgrc2_d_$$
        /usr/bin/time -v pgrc2 -d "$A" "$DEC" 2>"$TF2"
        read -r DWALL DRAMKB <<< "$(parse_time_v "$TF2")"
        LL=$(losscmp "$IN" "$DEC")
        rm -f "$TF" "$TF2" "$A" "$DEC"
        printf "TOOL=PgRC2 DS=%s raw=%s archive=%s ctime=%.2f dtime=%.2f cram_kb=%s dram_kb=%s lossless=%s\n" \
            "$DS" "$RAW" "$ARCH" "$CWALL" "$DWALL" "$CRAMKB" "$DRAMKB" "$LL" \
            | tee "$OUT_DIR/${DS}__PgRC2.log"
        pdone "PgRC2 — archive=${ARCH}B  ctime=${CWALL}s  dtime=${DWALL}s  $LL"
    else
        pskip "PgRC2 not on PATH"
    fi

    rm -f "$IN"
    phase "${DS_IDX}.7" "[$DS] dataset complete"
    pdone "$DS finished — all tools done"
done

# ── Summary table ─────────────────────────────────────────────────────────────
echo ""
echo "=== SUMMARY ==="
printf "%-15s %-8s %-12s %-10s %-10s %-10s\n" "DATASET" "TOOL" "ARCHIVE_B" "CTIME_S" "CRAM_MB" "LOSSLESS"
for LOG in "$OUT_DIR"/*.log; do
    [ -f "$LOG" ] || continue
    DS=$(grep -oP 'DS=\S+' "$LOG" | head -1 | cut -d= -f2)
    TOOL=$(grep -oP 'TOOL=\S+' "$LOG" | head -1 | cut -d= -f2)
    ARCH=$(grep -oP 'archive=\S+' "$LOG" | head -1 | cut -d= -f2)
    CT=$(grep -oP 'ctime=\S+' "$LOG" | head -1 | cut -d= -f2)
    CRAMKB=$(grep -oP ' cram_kb=\S+' "$LOG" | head -1 | tr -d ' ' | cut -d= -f2)
    LL=$(grep -oP 'lossless=\S+' "$LOG" | head -1 | cut -d= -f2)
    CRAMMB=$(awk "BEGIN{printf \"%.0f\", ${CRAMKB:-0}/1024}")
    printf "%-15s %-8s %-12s %-10s %-10s %-10s\n" "$DS" "$TOOL" "$ARCH" "$CT" "${CRAMMB}MB" "$LL"
done

echo ""
echo "ALL DONE $(date +%T) — logs in $OUT_DIR/"
