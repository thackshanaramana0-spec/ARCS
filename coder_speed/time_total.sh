#!/bin/bash
# Coder speedup via TOTAL-time delta. Only the coder differs across configs
# (load/assembly/quality identical), so mean(total_baseline) - mean(total_prune)
# == coder time saved. Big pg (ds7_300k) so the delta clears thermal noise.
# 4 runs/config, report mean total wall + pg_blob.
set -u
ARCS=./build_wsl/arcs
IN=/mnt/c/Temp/arcs_test/ds7_100k.fastq   # arg1 overrides
[ $# -ge 1 ] && IN="$1"
A=/tmp/tt.arcs

meas() {  # $1=orders -> "mean_total pg_blob"
  local s=0 t0 t1 pg last_pg
  for i in 1 2 3 4; do
    t0=$(date +%s.%N)
    pg=$(ARCSDNA_ORDERS="$1" $ARCS compress --chain-pg "$IN" "$A" 2>&1 | grep -oP 'pg_blob=\K[0-9]+')
    t1=$(date +%s.%N)
    s=$(awk "BEGIN{print $s + ($t1-$t0)}")
    last_pg=$pg
  done
  awk "BEGIN{printf \"%.2fs  pg_blob=$last_pg\", $s/4}"
}

echo "== TOTAL-time coder speedup :: $IN =="
printf "baseline  (8 models): %s\n" "$(meas 2,4,8,11,14,18,22,26)"
printf "prune-4   (7 models): %s\n" "$(meas 2,8,11,14,18,22,26)"
printf "prune-4-8 (6 models): %s\n" "$(meas 2,11,14,18,22,26)"
