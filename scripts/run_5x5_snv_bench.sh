#!/usr/bin/env bash
# run_5x5_snv_bench.sh — 5 chr20 windows × 5 GIAB individuals × 3 callers
# Streams reads from public GIAB remote BAMs; no pre-download required.
# Outputs one rtg vcfeval result directory per (individual, region, tool).
# Collect results afterward with:
#   python3 scripts/extract_vcfeval_metrics.py <OUTROOT> --out snv_table3.csv
#
# Usage:
#   bash run_5x5_snv_bench.sh <arcs_exe> <scripts_dir> <ref_grch37.fa> [outroot]
#
# Prerequisites (WSL): bwa, samtools, bcftools, bgzip, tabix, rtg-tools, DiscoSNP++
# Kmer2SNP and its VCF converter (scripts/kmer2snp_to_vcf.py) must be on PATH.
set -e
export PATH=~/miniconda3/bin:$PATH

ARCS="$1"; SC="$2"; REF="$3"; OUTROOT="${4:-~/bench_5x5}"
mkdir -p "$OUTROOT"

# ── Five chr20 windows (r2 = tuning; r3–r5+na = held-out) ──────────────────
declare -A WINDOWS
WINDOWS[r2]="3000000:3400000"   # tuning window — included for completeness
WINDOWS[r3]="4000000:4400000"
WINDOWS[na]="2000000:2400000"
WINDOWS[r4]="5000000:5400000"
WINDOWS[r5]="6000000:6400000"
WIN_ORDER=(r2 r3 na r4 r5)

# ── Five GIAB individuals with their public 300x BAM URLs (GRCh37) ──────────
declare -A BAM_URL
# HG001 (NA12878): no merged chr20 BAM on GIAB FTP; provide LOCAL_HG001_FQ env var to include it
# e.g. LOCAL_HG001_FQ=~/hg001_chr20_2_4M.fq bash run_5x5_snv_bench.sh ...
BAM_URL[HG002]="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x_chr20.bam"
BAM_URL[HG003]="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam"
BAM_URL[HG004]="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam"
BAM_URL[HG005]="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.hs37d5.300x.bam"
IND_ORDER=(HG002 HG003 HG004 HG005)

# ── GIAB truth VCF + BED per individual (GRCh37 v4.2.1) ─────────────────────
declare -A TRUTH_VCF TRUTH_BED
_GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
TRUTH_VCF[HG001]="$_GIAB/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG001]="$_GIAB/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.bed"
TRUTH_VCF[HG002]="$_GIAB/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG002]="$_GIAB/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH_VCF[HG003]="$_GIAB/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG003]="$_GIAB/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH_VCF[HG004]="$_GIAB/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG004]="$_GIAB/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH_VCF[HG005]="$_GIAB/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh37/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG005]="$_GIAB/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh37/HG005_GRCh37_1_22_v4.2.1_benchmark.bed"

# ── Shared helpers ───────────────────────────────────────────────────────────
CHROM=20

snv() { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
het() { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1];
        if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
prep_vcf(){
    local src="$1" tag="$2"
    eval "$3" < "$src" > "${tag}.v.vcf"
    { grep '^#' "${tag}.v.vcf"; grep -v '^#' "${tag}.v.vcf" | sort -k2,2n; } > "${tag}.s.vcf"
    bcftools norm -f "$REF" -m -any "${tag}.s.vcf" 2>/dev/null \
        | bcftools sort 2>/dev/null | bgzip > "${tag}.vcf.gz" \
        || bgzip -c "${tag}.s.vcf" > "${tag}.vcf.gz"
    tabix -f -p vcf "${tag}.vcf.gz"
}
score_vcfeval(){
    local label="$1" truth="$2" calls="$3" bed="$4" sdf="$5"
    local outd="e_${label}"   # separate local so label is already set
    rm -rf "$outd"
    rtg vcfeval -b "$truth" -c "$calls" -t "$sdf" \
        --squash-ploidy --bed-regions "$bed" -o "$outd" >/dev/null 2>&1 || true
    if [ -f "$outd/summary.txt" ]; then
        awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8}
             END{printf "%-14s TP=%-5s FP=%-5s FN=%-5s  P=%.3f R=%.3f F1=%.3f\n",
                        "'"$label"'",tp,fp,fn,p,r,f}' "$outd/summary.txt"
    else
        echo "$label: vcfeval produced no summary"
    fi
}

# ── BWA index once ───────────────────────────────────────────────────────────
[ -s "$REF.bwt" ] || { echo "[INDEX] bwa index $REF"; bwa index "$REF" 2>/dev/null; }

# ── Main loop ────────────────────────────────────────────────────────────────
for IND in "${IND_ORDER[@]}"; do
  for WIN in "${WIN_ORDER[@]}"; do
    RNG="${WINDOWS[$WIN]}"
    RLO="${RNG%%:*}"; RHI="${RNG##*:}"
    REGION="${CHROM}:${RLO}-${RHI}"
    WD="$OUTROOT/${IND}_${WIN}"
    mkdir -p "$WD"; cd "$WD"

    echo ""
    echo "========================================================"
    echo "  $IND  $WIN  chr${CHROM}:${RLO}-${RHI}"
    echo "========================================================"

    # 1. Stream reads (R1 only, 30x downsampled from 300x)
    # For HG002 na-window use pre-extracted local file if present (saves ~5min download)
    LOCAL_CACHE=""
    [ "$IND" = "HG002" ] && [ "$WIN" = "na" ] && LOCAL_CACHE="${LOCAL_HG002_NA:-}"
    if [ ! -s reads.fq ]; then
        if [ -n "$LOCAL_CACHE" ] && [ -s "$LOCAL_CACHE" ]; then
            echo "[READS] using local cache: $LOCAL_CACHE"
            cp "$LOCAL_CACHE" reads.fq
        else
            echo "[READS] streaming from remote BAM…"
            # -F 0x900: exclude secondary+supplementary; R1+R2 both included for ~30x
            # 1-in-10 from 300x paired (R1+R2) gives ~30x effective coverage
            samtools view -b -F 0x900 "${BAM_URL[$IND]}" "$REGION" \
                | samtools fastq -n - \
                | awk 'NR%4==1{r++} NR%4==2{s=length($0)} r%10==1{print}' \
                > reads.fq  # ~30x: keep 1-in-10 reads (R1+R2 interleaved) from 300x BAM
        fi
        echo "  reads: $(( $(wc -l < reads.fq) / 4 ))"
    fi

    # 2. ARCS call + dump contigs
    if [ ! -s calls_arcs.vcf ]; then
        export ARCS_DUMP_CONTIGS=contigs.tsv
        { /usr/bin/time -v "$ARCS" call reads.fq calls_arcs.vcf; } 2>arcs_time.txt || true
        awk -F'\t' '{print ">"$1"\n"$2}' contigs.tsv > contigs.fa
    fi

    # 3. BWA lift
    if [ ! -s lifted_arcs.vcf ]; then
        bwa mem -t4 "$REF" contigs.fa 2>/dev/null > c2r.sam
        python3 "$SC/lift_vcf.py" calls_arcs.vcf c2r.sam "$REF" $CHROM lifted_arcs.vcf contigs.fa
    fi

    # 4. DiscoSNP++
    if [ ! -s calls_disco.vcf ]; then
        echo "$WD/reads.fq" > fof.txt
        { /usr/bin/time -v \
          /home/btr/DiscoSnp/run_discoSnp++.sh -r fof.txt -k 31 -c 3 -D 100 -P 3 -b 0 -p disco; } \
          2>discosnp_time.txt || true
        cp disco*coherent.vcf calls_disco.vcf 2>/dev/null || touch calls_disco.vcf
    fi
    if [ ! -s lifted_disco.vcf ] && grep -q '^[^#]' calls_disco.vcf 2>/dev/null; then
        bwa mem -t4 "$REF" disco*coherent.fa 2>/dev/null > disco_c2r.sam 2>/dev/null || true
        python3 "$SC/lift_vcf.py" calls_disco.vcf disco_c2r.sam "$REF" $CHROM lifted_disco.vcf \
            disco*coherent.fa 2>/dev/null || cp calls_disco.vcf lifted_disco.vcf
    fi

    # 5. Kmer2SNP
    if [ ! -s calls_k2s.vcf ] && command -v kmer2snp >/dev/null 2>&1; then
        { /usr/bin/time -v kmer2snp reads.fq -o kmer2snp_out; } 2>k2s_time.txt || true
        python3 "$SC/kmer2snp_to_vcf.py" kmer2snp_out "$REF" $CHROM > calls_k2s.vcf 2>/dev/null || touch calls_k2s.vcf
    fi
    [ -s calls_k2s.vcf ] || touch calls_k2s.vcf

    # 6. Truth for this region
    if [ ! -s truth.vcf ]; then
        tabix -h "${TRUTH_VCF[$IND]}" "$REGION" 2>/dev/null | awk -v OFS='\t' '
          /^##contig/{next}
          /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
          /^#/{print; next} {print}' > truth.vcf
    fi

    # 7. Confident BED subset
    if [ ! -s giab_conf.bed ]; then
        curl -s "${TRUTH_BED[$IND]}" -o giab_conf.bed
    fi
    awk -v c=$CHROM -v lo=$RLO -v hi=$RHI \
        '$1==c && $3>lo && $2<hi{s=($2>lo?$2:lo); e=($3<hi?$3:hi); if(e>s) print c"\t"s"\t"e}' \
        giab_conf.bed > regions.bed

    # 8. rtg SDF
    [ -d sdf ] || rtg format -o sdf "$REF" >/dev/null 2>&1

    # 9. Prep and score
    prep_vcf truth.vcf       t_snv   'het | snv'
    prep_vcf lifted_arcs.vcf c_arcs  snv
    score_vcfeval arcs_snv  t_snv.vcf.gz c_arcs.vcf.gz regions.bed sdf

    if grep -q '^[^#]' lifted_disco.vcf 2>/dev/null; then
        prep_vcf lifted_disco.vcf c_disco snv
        score_vcfeval disco_snv t_snv.vcf.gz c_disco.vcf.gz regions.bed sdf
    fi

    if grep -q '^[^#]' calls_k2s.vcf 2>/dev/null; then
        prep_vcf calls_k2s.vcf c_k2s snv
        score_vcfeval k2s_snv t_snv.vcf.gz c_k2s.vcf.gz regions.bed sdf
    fi

    cd - >/dev/null
  done
done

echo ""
echo "All done. Collect results:"
echo "  python3 $SC/extract_vcfeval_metrics.py $OUTROOT --out snv_table3.csv"
