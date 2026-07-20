#!/bin/bash
# Per-model ablation of the pg context-mixing coder.
# For each config we record: pg_blob size (isolated coder ratio) + compress wall time.
# A model that costs time but contributes ~0 bytes = prune candidate (safe generalized win).
# Run from arcs_cpp/ so ./build_wsl/arcs resolves.
set -u
ARCS=./build_wsl/arcs
IN="${1:?usage: measure.sh <input.fastq> <tag>}"
TAG="${2:-run}"
OUT=/tmp/cs_${TAG}.arcs
BASE="2,4,8,11,14,18,22,26"

run() {   # $1=label  $2=orders
  local t0 t1 pg
  t0=$(date +%s.%N)
  pg=$(ARCSDNA_ORDERS="$2" $ARCS compress --chain-pg "$IN" "$OUT" 2>&1 | grep -oP 'pg_blob=\K[0-9]+')
  t1=$(date +%s.%N)
  printf "%-14s orders=%-28s pg_blob=%-9s time=%ss\n" "$1" "$2" "${pg:-ERR}" "$(awk "BEGIN{printf \"%.1f\", $t1-$t0}")"
}

echo "== $TAG :: $IN =="
run baseline "$BASE"
# leave-one-out: drop each order, see its marginal ratio contribution + time saved
for drop in 2 4 8 11 14 18 22 26; do
  ord=$(echo "$BASE" | tr ',' '\n' | grep -vx "$drop" | paste -sd,)
  run "drop-$drop" "$ord"
done
