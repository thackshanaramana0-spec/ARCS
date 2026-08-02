#!/usr/bin/env bash
# run_indel_strata.sh
# Re-score an existing ARCS + DiscoSNP++ indel vcfeval run with GIAB stratification
# BEDs to produce Table S4: TP/FP/FN/P/R/F1 broken down by:
#   (A) event type: insertion vs deletion
#   (B) indel length bin: 1bp / 2-5bp / 6+bp
#   (C) repeat context: homopolymer / STR / other
#
# Depends on GIAB stratification BEDs from:
#   https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/
# Usage:
#   bash run_indel_strata.sh <lifted_arcs.vcf> <lifted_disco.vcf> \
#        <truth_het_indel.vcf.gz> <regions_conf.bed> <ref.fa> <strat_dir> [outdir]
#
# strat_dir must contain (GRCh37 v3.3):
#   GRCh37_All_Homopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz
#   GRCh37_SimpleRepeat_diTR_51to200_slop5.bed.gz  (proxy for STR)
set -e
export PATH=~/miniconda3/bin:$PATH

ARCS_VCF="$1"; DISCO_VCF="$2"; TRUTH="$3"; CONF_BED="$4"
REF="$5"; STRAT="$6"; WD="${7:-~/indel_strata}"
mkdir -p "$WD"; cd "$WD"

CHROM=20
HOMOPOL="$STRAT/GRCh37_All_Homopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz"
STR_BED="$STRAT/GRCh37_SimpleRepeat_diTR_51to200_slop5.bed.gz"

# rtg SDF
[ -d sdf ] || rtg format -o sdf "$REF" >/dev/null 2>&1

indel(){ awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
ins(){   awk -F'\t' '/^#/{print;next} length($5)>length($4)'; }
del(){   awk -F'\t' '/^#/{print;next} length($4)>length($5)'; }
len1(){  awk -F'\t' '/^#/{print;next} {d=length($4)-length($5); if(d<0)d=-d; if(d==1) print}'; }
len25(){ awk -F'\t' '/^#/{print;next} {d=length($4)-length($5); if(d<0)d=-d; if(d>=2&&d<=5) print}'; }
len6p(){ awk -F'\t' '/^#/{print;next} {d=length($4)-length($5); if(d<0)d=-d; if(d>=6) print}'; }

prep(){ eval "$3" < "$1" > "$2.v.vcf"
  { grep '^#' "$2.v.vcf"; grep -v '^#' "$2.v.vcf"|sort -k2,2n; } > "$2.s.vcf"
  bcftools norm -f "$REF" -m -any "$2.s.vcf" 2>/dev/null|bcftools sort 2>/dev/null \
      |bgzip>"$2.vcf.gz" || bgzip -c "$2.s.vcf">"$2.vcf.gz"
  tabix -f -p vcf "$2.vcf.gz"; }

score(){
    local label="$1" truth="$2" calls="$3" bed="$4" outd="e_${label}"
    rm -rf "$outd"
    rtg vcfeval -b "$truth" -c "$calls" -t sdf --squash-ploidy \
        --bed-regions "$bed" -o "$outd" >/dev/null 2>&1 || true
    if [ -f "$outd/summary.txt" ]; then
        awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8}
             END{printf "%-32s TP=%-4s FP=%-4s FN=%-4s P=%.3f R=%.3f F1=%.3f\n",
                        "'"$label"'",tp,fp,fn,p,r,f}' "$outd/summary.txt"
    else echo "$label: no summary"; fi
}

# intersect a BED with conf_bed to get the strat-within-conf region
strat_bed(){
    local name="$1" sbf="$2"
    bedtools intersect -a "$CONF_BED" -b "$sbf" -u > "${name}.bed" 2>/dev/null || true
    echo "${name}.bed"
}

echo "Building stratification BEDs..."
HOMO_B=$(strat_bed homopolymer "$HOMOPOL")
STR_B=$(strat_bed str "$STR_BED")
# "other" = conf minus (homopol union str)
bedtools intersect -a "$CONF_BED" -b "$HOMOPOL" "$STR_BED" -v > other.bed 2>/dev/null || cp "$CONF_BED" other.bed

# Prep truth once per filter
prep "$TRUTH"  t_all    indel
prep "$TRUTH"  t_ins    ins
prep "$TRUTH"  t_del    del
prep "$TRUTH"  t_len1   len1
prep "$TRUTH"  t_len25  len25
prep "$TRUTH"  t_len6p  len6p

# Prep caller VCFs
prep "$ARCS_VCF"  c_arcs  indel
prep "$DISCO_VCF" c_disco indel

echo ""
echo "============================================================"
echo "  INDEL STRATA TABLE (rtg vcfeval, HG002 chr20:2.0-2.4Mb)"
echo "============================================================"
echo ""
echo "--- A. Event type ---"
for TOOL in arcs disco; do
    score "${TOOL}_ins"  t_ins.vcf.gz  "c_${TOOL}.vcf.gz"  "$CONF_BED"
    score "${TOOL}_del"  t_del.vcf.gz  "c_${TOOL}.vcf.gz"  "$CONF_BED"
done
echo ""
echo "--- B. Indel length ---"
prep "$ARCS_VCF"  c_arcs_l1   len1
prep "$ARCS_VCF"  c_arcs_l25  len25
prep "$ARCS_VCF"  c_arcs_l6p  len6p
prep "$DISCO_VCF" c_disco_l1  len1
prep "$DISCO_VCF" c_disco_l25 len25
prep "$DISCO_VCF" c_disco_l6p len6p
for TOOL in arcs disco; do
    for LEN in l1 l25 l6p; do
        score "${TOOL}_${LEN}" "t_${LEN}.vcf.gz" "c_${TOOL}_${LEN}.vcf.gz" "$CONF_BED"
    done
done
echo ""
echo "--- C. Repeat context ---"
for TOOL in arcs disco; do
    for CTX in homopolymer str other; do
        BED_F="${CTX}.bed"
        [ -s "$BED_F" ] || { echo "  $CTX bed empty; skipping"; continue; }
        score "${TOOL}_${CTX}"  t_all.vcf.gz  "c_${TOOL}.vcf.gz"  "$BED_F"
    done
done
echo ""
echo "Done. Results in: $WD"
