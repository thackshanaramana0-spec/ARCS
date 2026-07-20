#!/bin/bash
# Isolate the pg-CODER time (not assembly/quality) and its response to model count.
# ARCS_ENC_NOPAR unoverlaps async so the pg coding phase is measured alone.
# Average 3 runs to fight laptop thermal noise. Print the raw [ENC] labels once.
set -u
ARCS=./build_wsl/arcs
IN=/mnt/c/Temp/arcs_test/d1.fastq
A=/tmp/tc.arcs

echo "--- raw ENC labels (baseline, NOPAR) ---"
ARCS_ENC_NOPAR=1 ARCS_ENC_TIMING=1 ARCSDNA_ORDERS=2,4,8,11,14,18,22,26 \
  $ARCS compress --chain-pg "$IN" "$A" 2>&1 | grep -E '\[ENC\]|pg_blob'

avg_pg() {  # $1=orders  -> mean of "pg+names" (or pg) phase over 3 runs
  local s=0 v
  for i in 1 2 3; do
    v=$(ARCS_ENC_NOPAR=1 ARCS_ENC_TIMING=1 ARCSDNA_ORDERS="$1" \
        $ARCS compress --chain-pg "$IN" "$A" 2>&1 | grep -oP '\[ENC\] pg[^:]*: \K[0-9.]+' | head -1)
    s=$(awk "BEGIN{print $s+$v}")
  done
  awk "BEGIN{printf \"%.2f\", $s/3}"
}

echo "--- mean pg-coder phase (3 runs each) ---"
printf "baseline  (8 models): %ss\n" "$(avg_pg 2,4,8,11,14,18,22,26)"
printf "prune-4   (7 models): %ss\n" "$(avg_pg 2,8,11,14,18,22,26)"
printf "prune-4-8 (6 models): %ss\n" "$(avg_pg 2,11,14,18,22,26)"
