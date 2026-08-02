#!/bin/bash
# Benchmark ARCS, SPRING, Genozip, and fqzcomp on up to 8 FASTQ datasets.
#
# Usage:
#   bash benchmark/run_linux.sh <DATA_DIR> <ARCS_BIN> [OUT_DIR]
#
# DATA_DIR : directory containing FASTQ files named DS1.fastq … DS7.fastq,
#            GIAB.fastq, NA12878.fastq (any subset is fine; missing files are skipped).
# ARCS_BIN : path to the compiled arcs binary (e.g. build/arcs).
# OUT_DIR  : where to write per-dataset log files (default: ./bench_results).
#
# External tools expected on PATH: spring, genozip, genounzip, fqzcomp.
# Data is copied to /tmp for fair I/O timing on native Linux ext4.
#
# Example:
#   cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)
#   bash benchmark/run_linux.sh ~/fastq_datasets ./build/arcs ./results
set -euo pipefail

DATA_DIR="${1:?Usage: $0 <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
ARCS_BIN="${2:?Usage: $0 <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
OUT_DIR="${3:-./bench_results}"
mkdir -p "$OUT_DIR"

WD=/tmp/arcs_bench_wd
mkdir -p "$WD"

DATASETS="DS1 DS2 DS4 DS5 DS6 DS7 GIAB NA12878"

# best-of-2 wall time in seconds
btime() {
    local best=99999 t
    for _ in 1 2; do
        /usr/bin/time -o /tmp/_tt -f "%e" "$@" >/dev/null 2>/dev/null
        t=$(cat /tmp/_tt)
        awk -v t="$t" -v b="$best" 'BEGIN{exit !(t+0 < b+0)}' && best=$t
    done
    echo "$best"
}

# byte-exact lossless check (order-free, handles CRLF)
losscmp() {
    paste - - - - < "$1" | tr -d '\r' | sort > /tmp/_orig
    paste - - - - < "$2" | tr -d '\r' | sort > /tmp/_dec
    cmp -s /tmp/_orig /tmp/_dec && echo LOSSLESS || echo LOSSY
}

for DS in $DATASETS; do
    SRC="$DATA_DIR/$DS.fastq"
    [ -f "$SRC" ] || { echo "Skipping $DS (not found: $SRC)"; continue; }

    IN="$WD/$DS.fq"
    cp "$SRC" "$IN"
    RAW=$(stat -c %s "$IN")
    echo ">>> $DS  raw=$RAW B  $(date +%T)"

    # ── ARCS ──────────────────────────────────────────────────────────────────
    A="$WD/$DS.arcs"; DEC="$WD/$DS.arcs.dec"
    CT=$(btime "$ARCS_BIN" compress "$IN" "$A")
    ARCH=$(stat -c %s "$A")
    DT=$(btime "$ARCS_BIN" decompress "$A" "$DEC")
    LL=$(losscmp "$IN" "$DEC")
    printf "TOOL=ARCS DS=%s raw=%s archive=%s ctime=%s dtime=%s lossless=%s\n" \
        "$DS" "$RAW" "$ARCH" "$CT" "$DT" "$LL" | tee "$OUT_DIR/${DS}__ARCS.log"
    rm -f "$A" "$DEC"

    # ── SPRING ────────────────────────────────────────────────────────────────
    if command -v spring &>/dev/null; then
        A="$WD/$DS.spring"; DEC="$WD/$DS.spring.dec"
        CT=$(btime spring -c -i "$IN" -o "$A" -w "$WD")
        ARCH=$(stat -c %s "$A")
        DT=$(btime spring -d -i "$A" -o "$DEC" -w "$WD")
        LL=$(losscmp "$IN" "$DEC")
        printf "TOOL=SPRING DS=%s raw=%s archive=%s ctime=%s dtime=%s lossless=%s\n" \
            "$DS" "$RAW" "$ARCH" "$CT" "$DT" "$LL" | tee "$OUT_DIR/${DS}__SPRING.log"
        rm -f "$A" "$DEC"
    fi

    # ── Genozip ───────────────────────────────────────────────────────────────
    if command -v genozip &>/dev/null; then
        A="$WD/$DS.genozip"; DEC="$WD/$DS.gz.dec"
        CT=$(btime genozip --force -o "$A" "$IN")
        ARCH=$(stat -c %s "$A")
        DT=$(btime genounzip --force -o "$DEC" "$A")
        LL=$(losscmp "$IN" "$DEC")
        printf "TOOL=Genozip DS=%s raw=%s archive=%s ctime=%s dtime=%s lossless=%s\n" \
            "$DS" "$RAW" "$ARCH" "$CT" "$DT" "$LL" | tee "$OUT_DIR/${DS}__Genozip.log"
        rm -f "$A" "$DEC"
    fi

    # ── fqzcomp ───────────────────────────────────────────────────────────────
    if command -v fqzcomp &>/dev/null; then
        A="$WD/$DS.fqz"; DEC="$WD/$DS.fqz.dec"
        CT=$(btime fqzcomp "$IN" "$A")
        ARCH=$(stat -c %s "$A")
        DT=$(btime fqzcomp -d "$A" "$DEC")
        LL=$(losscmp "$IN" "$DEC")
        printf "TOOL=fqzcomp DS=%s raw=%s archive=%s ctime=%s dtime=%s lossless=%s\n" \
            "$DS" "$RAW" "$ARCH" "$CT" "$DT" "$LL" | tee "$OUT_DIR/${DS}__fqzcomp.log"
        rm -f "$A" "$DEC"
    fi

    rm -f "$IN"
done

echo "ALL DONE $(date +%T) — results in $OUT_DIR/"
