#!/usr/bin/env bash
# Full reference-free indel benchmark on synthetic diploid truth, scored by rtg vcfeval
# (the GA4GH engine hap.py uses) — same pipeline as the GIAB het-SNV numbers, with BWA
# placing ARCS's contigs to the reference for coordinate lift. Reports SNV and INDEL
# precision/recall/F1 separately.
#   usage: run_indel_bench.sh <workdir> <arcs_exe> <scripts_dir>
set -e
WD="$1"; ARCS="$2"; SC="$3"
CHROM=chrSIM
cd "$WD"

# 1. ARCS reference-free call (indels default on) + dump its contigs.
# NB: relative paths only — $ARCS is a Windows .exe under WSL interop and cannot resolve
# /mnt/c env-var paths (its CWD is already $WD via interop, so relative names work).
# WSL only forwards env vars named in WSLENV to Windows processes, so declare it there.
export ARCS_DUMP_CONTIGS="contigs.tsv"
export WSLENV="${WSLENV}:ARCS_DUMP_CONTIGS"
"$ARCS" call reads.fq calls.vcf 2> arcs_call.log
grep -E "CALL\] H" arcs_call.log || true
awk -F'\t' '{print ">"$1"\n"$2}' contigs.tsv > contigs.fa

# 2. Place contigs on the reference with BWA (evaluation-only coordinate lift)
bwa index ref.fa 2> bwa_index.log
bwa mem -t 4 ref.fa contigs.fa 2> bwa_mem.log > c2r.sam

# 3. Lift contig-coordinate calls to genome coordinates (contigs.fa enables exact
#    haplotype-flank indel lift: resolves ins/del polarity + strand repeat-robustly)
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam ref.fa $CHROM lifted.vcf contigs.fa

# 4. Split truth + calls into SNV-only and INDEL-only, score each with rtg vcfeval
snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }

rtg format -o sdf ref.fa > /dev/null 2>&1 || true
NORM=${NORM:-1}    # left-align + canonicalize both truth and calls before scoring
prep() { # <in.vcf> <out_prefix> <filter>
  $3 < "$1" > "$2.raw.vcf"
  if [ "$NORM" = 1 ] && command -v bcftools >/dev/null 2>&1; then
    bcftools norm -f ref.fa -c s "$2.raw.vcf" 2>/dev/null | bcftools sort 2>/dev/null > "$2.vcf" || cp "$2.raw.vcf" "$2.vcf"
  else cp "$2.raw.vcf" "$2.vcf"; fi
  bgzip -f "$2.vcf"; tabix -f -p vcf "$2.vcf.gz"; }

score() { # <name> <truth.vcf.gz> <calls.vcf.gz>
  rm -rf "eval_$1"
  rtg vcfeval -b "$2" -c "$3" -t sdf --squash-ploidy --bed-regions regions.bed \
      -o "eval_$1" > /dev/null 2>&1 || true
  if [ -f "eval_$1/summary.txt" ]; then
    # rtg summary: last row "None" has TP-baseline FP FN Precision Recall F-measure
    awk 'NR>2{tp=$2;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{printf "%-6s TP=%s FP=%s FN=%s  P=%.3f R=%.3f F1=%.3f\n","'$1'",tp,fp,fn,p,r,f}' "eval_$1/summary.txt"
  else echo "$1: vcfeval produced no summary (see eval_$1)"; fi
}

prep truth.vcf  truth_snv   snv
prep truth.vcf  truth_indel indel
prep lifted.vcf call_snv    snv
prep lifted.vcf call_indel  indel

echo "=================== rtg vcfeval (gold-standard) ==================="
score SNV   truth_snv.vcf.gz   call_snv.vcf.gz
score INDEL truth_indel.vcf.gz call_indel.vcf.gz
echo "==================================================================="
