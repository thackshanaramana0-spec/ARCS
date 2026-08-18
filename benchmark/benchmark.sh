#!/bin/bash
# ARCS Benchmark Suite — dispatcher
#
# Usage:
#   bash benchmark/benchmark.sh <claim1|claim2|claim3|all> <DATA_DIR> <ARCS_BIN> [OUT_DIR]
#
# claim1 : Compression ratio + speed + RAM (T1+T2) — any .fq files in DATA_DIR
# claim2 : Variant calling SNV/indel (T3+T4+T5)  — needs GIAB truth + chr20 FASTQ
# claim3 : Archive-derived analysis (T6)          — needs assembled datasets
# all    : Runs claim1 → claim2 → claim3 in order
#
# DATA_DIR : directory containing input FASTQ files (.fq or .fastq)
# ARCS_BIN : path to compiled arcs binary
# OUT_DIR  : results output directory (default: ./results)

set -euo pipefail

CLAIM="${1:?Usage: $0 <claim1|claim2|claim3|all> <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
DATA_DIR="${2:?Usage: $0 <claim1|claim2|claim3|all> <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
ARCS_BIN="${3:?Usage: $0 <claim1|claim2|claim3|all> <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
OUT_DIR="${4:-./results}"

BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Validate inputs ───────────────────────────────────────────────────────────
[ -d "$DATA_DIR" ]  || { echo "ERROR: DATA_DIR not found: $DATA_DIR"; exit 1; }
[ -x "$ARCS_BIN" ] || { echo "ERROR: ARCS binary not found or not executable: $ARCS_BIN"; exit 1; }

mkdir -p "$OUT_DIR"

echo "=========================================="
echo "  ARCS Benchmark Suite"
echo "  claim    : $CLAIM"
echo "  data_dir : $DATA_DIR"
echo "  arcs_bin : $ARCS_BIN"
echo "  out_dir  : $OUT_DIR"
echo "  started  : $(date)"
echo "=========================================="
echo ""

# ── Dispatcher ────────────────────────────────────────────────────────────────
run_claim1() {
    echo ">>> CLAIM 1: Compression ratio + speed + RAM"
    bash "$BENCH_DIR/run_block1.sh" "$DATA_DIR" "$ARCS_BIN" "$OUT_DIR/claim1"
    echo "<<< CLAIM 1 done"
}

run_claim2() {
    echo ">>> CLAIM 2: Variant calling (T3+T4+T5)"
    bash "$BENCH_DIR/run_claim2.sh" "$DATA_DIR" "$ARCS_BIN" "$OUT_DIR/claim2"
    echo "<<< CLAIM 2 done"
}

run_claim3() {
    echo ">>> CLAIM 3: Archive-derived analysis (T6)"
    bash "$BENCH_DIR/run_claim3.sh" "$DATA_DIR" "$ARCS_BIN" "$OUT_DIR/claim3"
    echo "<<< CLAIM 3 done"
}

if   [ "$CLAIM" = "claim1" ]; then run_claim1
elif [ "$CLAIM" = "claim2" ]; then run_claim2
elif [ "$CLAIM" = "claim3" ]; then run_claim3
elif [ "$CLAIM" = "all"    ]; then run_claim1; run_claim2; run_claim3
else
    echo "ERROR: unknown claim '$CLAIM' — must be claim1, claim2, claim3, or all"
    exit 1
fi

echo ""
echo "=========================================="
echo "  ALL DONE: $(date)"
echo "  Results in: $OUT_DIR"
echo "=========================================="
