#!/usr/bin/env bash
# REAL-GIAB reference-free indel benchmark, scored by rtg vcfeval (the GA4GH engine),
# same methodology as the het-SNV numbers: HG002 chr20:2.0-2.4Mb Illumina ~30x reads,
# GIAB v4.2.1 truth restricted to the read region ∩ the GIAB confident BED. ARCS calls
# are reference-free; BWA places its contigs to the reference for coordinate lift only.
# Truth is fetched for just the region via remote tabix (no multi-GB download); the
# confident BED is downloaded once and subset. SNV and INDEL are scored separately.
#   usage: run_giab_indel.sh <arcs_exe> <scripts_dir> <reads.fq> <ref.fa>
set -e
export PATH=~/miniconda3/bin:$PATH        # bcftools (for GA4GH left-alignment norm)
ARCS="$1"; SC="$2"; READS="$3"; REF="$4"
CHROM=20; RLO=2000000; RHI=2400000; REGION="$CHROM:$RLO-$RHI"
WD=~/giab_indel; mkdir -p "$WD"; cd "$WD"
FTP=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37
VCFURL="$FTP/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BEDURL="$FTP/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"

# 1. ARCS reference-free call + dump contigs (relative paths only for the Win exe)
export ARCS_DUMP_CONTIGS=contigs.tsv WSLENV=ARCS_DUMP_CONTIGS
[ -s calls.vcf ] || "$ARCS" call "$READS" calls.vcf 2>&1 | grep -E "CALL\] H" || true
awk -F'\t' '{print ">"$1"\n"$2}' contigs.tsv > contigs.fa

# 2. Place contigs (eval-only coordinate lift)
[ -s "$REF.bwt" ] || bwa index "$REF" 2>/dev/null
bwa mem -t4 "$REF" contigs.fa 2>/dev/null > c2r.sam

# 3. Lift contig-coordinate calls to genome coords (exact haplotype-flank indel lift)
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam "$REF" $CHROM lifted.vcf contigs.fa

# 4. Truth for the region (remote tabix — a few KB), header patched with contig length
tabix -h "$VCFURL" "$REGION" 2>/dev/null | awk -v OFS='\t' '
  /^##contig/{next} /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
  /^#/{print; next} {print}' > truth.vcf

# 5. Confident BED subset to the eval region (download once)
[ -s giab_conf.bed ] || curl -s "$BEDURL" -o giab_conf.bed
awk -v c=$CHROM -v lo=$RLO -v hi=$RHI '$1==c && $3>lo && $2<hi{
  s=($2>lo?$2:lo); e=($3<hi?$3:hi); if(e>s) print c"\t"s"\t"e}' giab_conf.bed > regions.bed
echo "confident bp in region: $(awk '{n+=$3-$2}END{print n}' regions.bed)"

# 6. rtg SDF + split SNV/indel + score each restricted to the confident region
rm -rf sdf; rtg format -o sdf "$REF" > /dev/null 2>&1
snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
# heterozygous-only: reference-free bubble calling sees ONLY het variants (a homozygous
# variant forms no bubble), so — as in the het-SNV methodology and Kmer2SNP's scope —
# the truth is restricted to heterozygous sites. Single-ALT het GT = 0/1|1/0|0|1|1|0.
het()   { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1]; if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
# prep: filter -> sort -> left-align + split-multiallelic against the reference (GA4GH
# standard, so homopolymer/STR indels have a canonical representation) -> bgzip+index.
prep(){ eval "$3" < "$1" > "$2.v.vcf"
  (grep '^#' "$2.v.vcf"; grep -v '^#' "$2.v.vcf"|sort -k2,2n) > "$2.s.vcf"
  bcftools norm -f "$REF" -m -any "$2.s.vcf" 2>/dev/null | bcftools sort 2>/dev/null | bgzip > "$2.vcf.gz" \
    || { bgzip -c "$2.s.vcf" > "$2.vcf.gz"; }
  tabix -f -p vcf "$2.vcf.gz"; }
prep truth.vcf  t_snv   'het | snv'
prep truth.vcf  t_ind   'het | indel'
prep lifted.vcf c_snv   snv
prep lifted.vcf c_ind   indel
score(){ rm -rf "e_$1"; rtg vcfeval -b "$2" -c "$3" -t sdf --squash-ploidy \
    --bed-regions regions.bed -o "e_$1" >/dev/null 2>&1 || true
  if [ -f "e_$1/summary.txt" ]; then
    awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8} END{printf "%-6s TP=%s FP=%s FN=%s  P=%.3f R=%.3f F1=%.3f\n","'$1'",tp,fp,fn,p,r,f}' "e_$1/summary.txt"
  else echo "$1: no summary"; fi; }
echo "======== REAL-GIAB HG002 chr20:2.0-2.4M — rtg vcfeval (gold-standard) ========"
score SNV   t_snv.vcf.gz c_snv.vcf.gz
score INDEL t_ind.vcf.gz c_ind.vcf.gz
echo "=============================================================================="
echo "truth het-indel in region: $(zcat t_ind.vcf.gz|grep -vc '^#')   arcs indel calls: $(zcat c_ind.vcf.gz|grep -vc '^#')"
