#!/usr/bin/env bash
# run_fused_overhead_all.sh
# Measure "arcs compress" vs "arcs compress --call" wall time on every benchmark
# dataset to report calling overhead per dataset (not just GIAB_HG002).
#
# Usage:
#   bash run_fused_overhead_all.sh <arcs_exe> <datasets_dir> [runs=3]
#
# datasets_dir must contain one .fq/.fastq file per dataset, named after the
# dataset (e.g. ECOLI_35x.fq, GIAB_HG002.fq, …).
# Output: fused_overhead.csv
set -e
ARCS="$1"; DDIR="$2"; RUNS="${3:-3}"

DATASETS=(ECOLI_35x HUMAN_127bp MTB_51bp SARS2_AMP GIAB_HG002 GIAB_HG001 ECOLI_30x HUMAN_30x)

OUT=fused_overhead.csv
echo "Dataset,Run,Compress_s,Fused_s,Overhead_s,Overhead_pct,SNV_calls,Indel_calls" > "$OUT"

time_cmd(){ /usr/bin/time -f '%e' "$@" 2>&1 | tail -1; }

for DS in "${DATASETS[@]}"; do
    FQ=$(ls "$DDIR/${DS}".{fq,fastq,fq.gz,fastq.gz} 2>/dev/null | head -1)
    if [ -z "$FQ" ]; then
        echo "SKIP $DS — no FASTQ found in $DDIR"
        continue
    fi
    echo ""
    echo "=== $DS  ($FQ) ==="
    for i in $(seq 1 "$RUNS"); do
        # plain compress
        T_C=$( { /usr/bin/time -f '%e' "$ARCS" compress "$FQ" /tmp/arcs_bench.arcs; } 2>&1 | tail -1 )
        rm -f /tmp/arcs_bench.arcs

        # fused compress --call
        T_F=$( { /usr/bin/time -f '%e' "$ARCS" compress --call /tmp/arcs_call.vcf "$FQ" /tmp/arcs_bench.arcs; } 2>&1 | tail -1 )
        NSNV=$(grep -vc '^#' /tmp/arcs_call.vcf 2>/dev/null | tr -d ' ' || echo 0)
        # count SVTYPE=INDEL records
        NIND=$(grep -c 'SVTYPE=INDEL' /tmp/arcs_call.vcf 2>/dev/null || echo 0)
        NSNV=$(( NSNV - NIND ))
        rm -f /tmp/arcs_bench.arcs /tmp/arcs_call.vcf

        OVH=$(awk "BEGIN{printf \"%.3f\", $T_F - $T_C}")
        PCT=$(awk "BEGIN{printf \"%.1f\", 100*($T_F - $T_C)/($T_C+0.001)}")

        echo "  run $i: compress=${T_C}s  fused=${T_F}s  overhead=${OVH}s (${PCT}%)  SNV=${NSNV}  indel=${NIND}"
        echo "$DS,$i,$T_C,$T_F,$OVH,$PCT,$NSNV,$NIND" >> "$OUT"
    done
done

echo ""
echo "Written: $OUT"
echo ""
# Summary: mean overhead per dataset
echo "=== Mean overhead summary ==="
printf "%-14s  %8s  %8s  %8s  %7s\n" "Dataset" "Compress" "Fused" "Overhead" "Ovhd%"
awk -F',' 'NR==1{next}
{
  ds=$1; tc+=$3; tf+=$4; ov+=$5; pct+=$6; n++
  prev=ds
  if(NR==last || ds!=prev_ds) {
    printf "%-14s  %8.2f  %8.2f  %8.3f  %6.1f%%\n", prev_ds, tc/n, tf/n, ov/n, pct/n
    tc=0;tf=0;ov=0;pct=0;n=0
  }
  prev_ds=ds
}' "$OUT"
