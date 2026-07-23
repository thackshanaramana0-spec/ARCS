#!/usr/bin/env bash
# Polyploid reference-free SNV benchmark: arcs call (ARCS_PLOIDY=k) -> BWA place contigs ->
# multi-allelic lift -> rtg vcfeval (--squash-ploidy = allele-presence match). Scores how
# many planted multi-allelic/biallelic het SNVs the k-ploid caller recovers.
#   usage: run_polyploid_bench.sh <workdir> <arcs_exe> <scripts_dir> <ploidy>
set -e
export PATH=~/miniconda3/bin:$PATH
WD="$1"; ARCS="$2"; SC="$3"; K="${4:-3}"; CHROM=chrPLD
cd "$WD"
export ARCS_DUMP_CONTIGS=contigs.tsv ARCS_PLOIDY=$K WSLENV=ARCS_DUMP_CONTIGS:ARCS_PLOIDY
"$ARCS" call reads.fq calls.vcf 2>&1 | grep -E "CALL\] H" || true
awk -F'\t' '{print ">"$1"\n"$2}' contigs.tsv > contigs.fa
[ -s ref.fa.bwt ] || bwa index ref.fa 2>/dev/null
bwa mem -t4 ref.fa contigs.fa 2>/dev/null > c2r.sam
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam ref.fa $CHROM lifted.vcf contigs.fa
rm -rf sdf; rtg format -o sdf ref.fa >/dev/null 2>&1
prep(){ (grep '^#' "$1"; grep -v '^#' "$1"|sort -k2,2n) | bgzip > "$2.vcf.gz"; tabix -f -p vcf "$2.vcf.gz"; }
prep truth.vcf  truth
prep lifted.vcf calls
rm -rf e_pld
rtg vcfeval -b truth.vcf.gz -c calls.vcf.gz -t sdf --squash-ploidy --bed-regions regions.bed -o e_pld >/dev/null 2>&1 || true
echo "======== POLYPLOID (k=$K) reference-free SNV — rtg vcfeval (--squash-ploidy) ========"
awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{printf "k=%s  TP=%s FP=%s FN=%s  P=%.3f R=%.3f F1=%.3f\n","'$K'",tp,fp,fn,p,r,f}' e_pld/summary.txt
echo "truth het SNV: $(zcat truth.vcf.gz|grep -vc '^#')   multiallelic truth: $(zcat truth.vcf.gz|grep -v '^#'|awk -F'\t' '$5~/,/'|wc -l)   calls: $(zcat calls.vcf.gz|grep -vc '^#')"
