#!/bin/bash
# Combined-prune confirmation: measure the dead-model prune as a real config change,
# check it stays ~0 ratio on the COMBINED drop (leave-one-out can miss redundancy),
# isolate coder time (NOPAR unoverlaps async), and verify LOSSLESS roundtrip
# (orders must apply to BOTH encode and decode → env set on both processes).
set -u
ARCS=./build_wsl/arcs
IN="${1:?usage: confirm.sh <input.fastq> <tag>}"
TAG="${2:-run}"
A=/tmp/cf_${TAG}.arcs
D=/tmp/cf_${TAG}.dec

run() {  # $1=label $2=orders
  local pg t0 t1 arch dec_ok
  t0=$(date +%s.%N)
  pg=$(ARCS_ENC_NOPAR=1 ARCS_ENC_TIMING=1 ARCSDNA_ORDERS="$2" $ARCS compress --chain-pg "$IN" "$A" 2>/tmp/enc.log; grep -oP 'pg_blob=\K[0-9]+' /tmp/enc.log)
  t1=$(date +%s.%N)
  local pgt=$(grep -oP 'pg\+names.*: \K[0-9.]+' /tmp/enc.log)
  arch=$(stat -c%s "$A")
  # lossless roundtrip: decode with SAME orders, compare to input
  ARCSDNA_ORDERS="$2" $ARCS decompress "$A" "$D" 2>/dev/null
  if cmp -s "$IN" "$D"; then dec_ok="LOSSLESS"; else dec_ok="DIFF"; fi
  printf "%-12s orders=%-26s pg_blob=%-8s archive=%-9s pg_time=%-5ss enc=%.1fs %s\n" \
    "$1" "$2" "${pg:-ERR}" "$arch" "${pgt:-?}" "$(awk "BEGIN{printf \"%.1f\",$t1-$t0}")" "$dec_ok"
}

echo "== CONFIRM $TAG :: $IN =="
run baseline   "2,4,8,11,14,18,22,26"
run prune-4    "2,8,11,14,18,22,26"
run prune-4-8  "2,11,14,18,22,26"
