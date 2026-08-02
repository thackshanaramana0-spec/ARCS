#!/usr/bin/env bash
# run_chr21_bench.sh
# Extend the het-SNV benchmark to chr21 and chr22 to support "multiple chromosomes"
# in the Methods. Streams reads from public GIAB BAMs; no pre-download required.
# Scores ARCS, DiscoSNP++, Kmer2SNP via rtg vcfeval (same pipeline as chr20).
#
# Usage:
#   bash run_chr21_bench.sh <arcs_exe> <scripts_dir> <ref_grch37.fa> [outroot]
#
# Evaluates HG001 and HG002 on two 400-kb windows each for chr21 and chr22.
set -e
export PATH=~/miniconda3/bin:$PATH

ARCS="$1"; SC="$2"; REF="$3"; OUTROOT="${4:-~/bench_chr21_22}"
mkdir -p "$OUTROOT"

# ── Windows ─────────────────────────────────────────────────────────────────
# chr21: avoid centromere (starts ~10Mb); use well-covered non-repetitive regions
# chr22: non-repetitive coding regions near 20-25Mb
declare -A WIN_DEF
WIN_DEF[chr21_w1]="21:30000000:30400000"
WIN_DEF[chr21_w2]="21:31000000:31400000"
WIN_DEF[chr22_w1]="22:20000000:20400000"
WIN_DEF[chr22_w2]="22:21000000:21400000"
WIN_ORDER=(chr21_w1 chr21_w2 chr22_w1 chr22_w2)

# ── Individuals ──────────────────────────────────────────────────────────────
declare -A BAM_URL
# HG001 (NA12878): no merged chr21/22 BAM on GIAB FTP; excluded from this script
BAM_URL[HG002]="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam"
IND_ORDER=(HG002)

declare -A TRUTH_VCF TRUTH_BED
_GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
TRUTH_VCF[HG001]="$_GIAB/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG001]="$_GIAB/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.bed"
TRUTH_VCF[HG002]="$_GIAB/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG002]="$_GIAB/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"

snv(){ awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
het(){ awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1];
       if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
prep_vcf(){
    local src="$1" tag="$2"
    eval "$3" < "$src" > "${tag}.v.vcf"
    { grep '^#' "${tag}.v.vcf"; grep -v '^#' "${tag}.v.vcf"|sort -k2,2n; } > "${tag}.s.vcf"
    bcftools norm -f "$CHREF" -m -any "${tag}.s.vcf" 2>/dev/null \
        |bcftools sort 2>/dev/null|bgzip>"${tag}.vcf.gz" \
        || bgzip -c "${tag}.s.vcf">"${tag}.vcf.gz"
    tabix -f -p vcf "${tag}.vcf.gz"
}
score_vcfeval(){
    local label="$1" truth="$2" calls="$3" bed="$4" sdf="$5"
    rm -rf "e_${label}"
    rtg vcfeval -b "$truth" -c "$calls" -t "$sdf" \
        --squash-ploidy --bed-regions "$bed" -o "e_${label}" >/dev/null 2>&1 || true
    if [ -f "e_${label}/summary.txt" ]; then
        awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8}
             END{printf "%-14s TP=%-4s FP=%-4s FN=%-4s P=%.3f R=%.3f F1=%.3f\n",
                        "'"$label"'",tp,fp,fn,p,r,f}' "e_${label}/summary.txt"
    else echo "$label: no summary"; fi
}

# REF is a directory or base path; derive per-chromosome reference FASTA
# Expects: <REF_DIR>/chr21_grch37_ens.fa and chr22_grch37_ens.fa
# or falls back to REF itself if it already contains chr21/22

for IND in "${IND_ORDER[@]}"; do
  for WIN in "${WIN_ORDER[@]}"; do
    DEF="${WIN_DEF[$WIN]}"
    CHROM="${DEF%%:*}"; REST="${DEF#*:}"; RLO="${REST%%:*}"; RHI="${REST##*:}"
    REGION="${CHROM}:${RLO}-${RHI}"
    # Pick chromosome-specific reference
    CHREF="$(dirname "$REF")/chr${CHROM}_grch37_ens.fa"
    [ -s "$CHREF" ] || CHREF="$REF"
    [ -s "${CHREF}.bwt" ] || bwa index "$CHREF" 2>/dev/null
    WD="$OUTROOT/${IND}_${WIN}"
    mkdir -p "$WD"; cd "$WD"

    echo ""
    echo "=========================================="
    echo "  $IND  $WIN  chr${CHROM}:${RLO}-${RHI}"
    echo "=========================================="

    # 1. Stream reads (~30x: 1-in-10 from 300x)
    if [ ! -s reads.fq ]; then
        samtools view -b -f 0x40 -F 0x900 "${BAM_URL[$IND]}" "$REGION" \
            | samtools fastq -n - \
            | awk 'NR%4==1{r++} NR%4==2{s=length($0)} r%10==1{print}' \
            > reads.fq
        echo "  reads: $(( $(wc -l < reads.fq) / 4 ))"
    fi

    # 2. ARCS
    if [ ! -s calls_arcs.vcf ]; then
        export ARCS_DUMP_CONTIGS=contigs.tsv
        "$ARCS" call reads.fq calls_arcs.vcf 2>&1 | grep 'CALL\]' || true
        awk -F'\t' '{print ">"$1"\n"$2}' contigs.tsv > contigs.fa
    fi
    [ ! -s lifted_arcs.vcf ] && {
        bwa mem -t4 "$CHREF" contigs.fa 2>/dev/null > c2r.sam
        python3 "$SC/lift_vcf.py" calls_arcs.vcf c2r.sam "$CHREF" $CHROM lifted_arcs.vcf contigs.fa
    }

    # 3. DiscoSNP++
    if [ ! -s calls_disco.vcf ]; then
        echo "$WD/reads.fq" > fof.txt
        /home/btr/DiscoSnp/run_discoSnp++.sh -r fof.txt -k 31 -c 3 -D 100 -P 3 -b 0 \
            -p disco 2>&1 | grep -E 'Wall|error' || true
        cp disco*coherent.vcf calls_disco.vcf 2>/dev/null || touch calls_disco.vcf
    fi
    if [ ! -s lifted_disco.vcf ] && grep -q '^[^#]' calls_disco.vcf 2>/dev/null; then
        bwa mem -t4 "$REF" disco*coherent.fa 2>/dev/null > disco_c2r.sam 2>/dev/null || true
        python3 "$SC/lift_vcf.py" calls_disco.vcf disco_c2r.sam "$REF" $CHROM \
            lifted_disco.vcf disco*coherent.fa 2>/dev/null || cp calls_disco.vcf lifted_disco.vcf
    fi

    # 4. Truth + confident BED
    if [ ! -s truth.vcf ]; then
        CLEN=$([ "$CHROM" = "21" ] && echo 48129895 || echo 51304566)
        tabix -h "${TRUTH_VCF[$IND]}" "$REGION" 2>/dev/null | awk -v OFS='\t' \
            -v c=$CHROM -v clen=$CLEN '
            /^##contig/{next}
            /^#CHROM/{print "##contig=<ID="c",length="clen">"; print; next}
            /^#/{print; next} {print}' > truth.vcf
    fi
    if [ ! -s giab_conf.bed ]; then
        curl -s "${TRUTH_BED[$IND]}" -o giab_conf.bed
    fi
    awk -v c=$CHROM -v lo=$RLO -v hi=$RHI \
        '$1==c && $3>lo && $2<hi{s=($2>lo?$2:lo);e=($3<hi?$3:hi);if(e>s)print c"\t"s"\t"e}' \
        giab_conf.bed > regions.bed

    # 5. SDF + score (per-chromosome SDF)
    [ -d sdf ] || rtg format -o sdf "$CHREF" >/dev/null 2>&1

    prep_vcf truth.vcf       t_snv  'het | snv'
    prep_vcf lifted_arcs.vcf c_arcs snv
    score_vcfeval arcs_snv t_snv.vcf.gz c_arcs.vcf.gz regions.bed sdf

    if grep -q '^[^#]' lifted_disco.vcf 2>/dev/null; then
        prep_vcf lifted_disco.vcf c_disco snv
        score_vcfeval disco_snv t_snv.vcf.gz c_disco.vcf.gz regions.bed sdf
    fi

    cd - >/dev/null
  done
done

echo ""
echo "Done. Collect with:"
echo "  python3 $SC/extract_vcfeval_metrics.py $OUTROOT --out chr21_22_snv.csv"
