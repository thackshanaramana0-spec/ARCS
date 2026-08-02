#!/usr/bin/env bash
# Generalised indel + SNV benchmark for any HG002 chr20 region.
# Streams reads directly from the public GIAB remote BAM — no pre-download needed.
# Scores both ARCS and DiscoSNP++ via rtg vcfeval (GA4GH engine).
#
# Usage:
#   run_indel_region.sh <arcs_exe> <scripts_dir> <ref.fa> <RLO> <RHI> [<outdir>]
#
# Example (region 2):
#   run_indel_region.sh /path/arcs /path/scripts chr20_grch37.fa 10000000 10400000 ~/indel_r2
set -e
export PATH=~/miniconda3/bin:$PATH

ARCS="$1"; SC="$2"; REF="$3"; RLO="$4"; RHI="$5"
WD="${6:-~/indel_r${RLO}}"
CHROM=20
REGION="${CHROM}:${RLO}-${RHI}"
mkdir -p "$WD"; cd "$WD"

HG002_BAM="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x_chr20.bam"
FTP=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37
VCFURL="$FTP/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
BEDURL="$FTP/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"

echo "[REGION] chr${CHROM}:${RLO}-${RHI}  WD=$WD"

# 1. Stream reads from remote BAM → FASTQ (R1 only, single-end mode)
if [ ! -s reads.fq ]; then
    echo "[STEP 1] Streaming reads from remote GIAB BAM..."
    # R1-only: mixing R1+R2 confuses SNV bubble calling (R2 is reverse-complement).
    # Pre-filter to READ1 (0x40) before fastq conversion so no pair-matching runs —
    # without this, mates outside the queried window become "singletons" dropped to
    # /dev/null, yielding <2% of expected reads.
    samtools view -b -f 0x40 -F 0x900 "$HG002_BAM" "${CHROM}:${RLO}-${RHI}" \
        | samtools fastq -n - > reads.fq
    NREADS=$(( $(wc -l < reads.fq) / 4 ))
    echo "reads=${NREADS}"
    if [ "$NREADS" -lt 5000 ]; then
        echo "WARNING: only ${NREADS} reads — region may be repetitive; try a different window"
    fi
fi

# 2. ARCS reference-free call
export ARCS_DUMP_CONTIGS=contigs.tsv
[ -s calls.vcf ] || "$ARCS" call reads.fq calls.vcf 2>&1 | grep -E 'CALL\]|error' || true
awk -F'\t' '{print ">"$1"\n"$2}' contigs.tsv > contigs.fa
echo "[STEP 2] ARCS: $(grep -c '^[^#]' calls.vcf 2>/dev/null || echo 0) raw calls"

# 3. BWA place contigs to reference for coordinate lift
[ -s "$REF.bwt" ] || bwa index "$REF" 2>/dev/null
bwa mem -t4 "$REF" contigs.fa 2>/dev/null > c2r.sam

# 4. Lift contig calls to genome coords
python3 "$SC/lift_vcf.py" calls.vcf c2r.sam "$REF" $CHROM lifted.vcf contigs.fa

# 5. Truth for the region (remote tabix)
tabix -h "$VCFURL" "$REGION" 2>/dev/null | awk -v OFS='\t' '
  /^##contig/{next} /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
  /^#/{print; next} {print}' > truth.vcf

# 6. Confident BED subset
[ -s giab_conf.bed ] || curl -s "$BEDURL" -o giab_conf.bed
awk -v c=$CHROM -v lo=$RLO -v hi=$RHI '$1==c && $3>lo && $2<hi{
  s=($2>lo?$2:lo); e=($3<hi?$3:hi); if(e>s) print c"\t"s"\t"e}' giab_conf.bed > regions.bed
echo "confident bp: $(awk '{n+=$3-$2}END{print n}' regions.bed)"

# 7. DiscoSNP++ on same reads
if [ ! -s disco_coherent.vcf ]; then
    echo "[STEP 7] Running DiscoSNP++..."
    echo "$WD/reads.fq" > fof.txt
    /home/btr/DiscoSnp/run_discoSnp++.sh -r fof.txt -k 31 -c 3 -D 100 -P 3 -b 0 \
        -p disco 2>&1 | grep -E 'Wall|error|Error' || true
    cp disco*coherent.vcf disco_coherent.vcf 2>/dev/null || touch disco_coherent.vcf
fi

# 8. Score both via rtg vcfeval
rm -rf sdf; rtg format -o sdf "$REF" > /dev/null 2>&1

snv()   { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }
het()   { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1];
          if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }
prep(){ eval "$3" < "$1" > "$2.v.vcf"
  (grep '^#' "$2.v.vcf"; grep -v '^#' "$2.v.vcf"|sort -k2,2n) > "$2.s.vcf"
  bcftools norm -f "$REF" -m -any "$2.s.vcf" 2>/dev/null | bcftools sort 2>/dev/null \
      | bgzip > "$2.vcf.gz" || bgzip -c "$2.s.vcf" > "$2.vcf.gz"
  tabix -f -p vcf "$2.vcf.gz"; }

prep truth.vcf  t_snv 'het | snv'
prep truth.vcf  t_ind 'het | indel'
prep lifted.vcf c_snv snv
prep lifted.vcf c_ind indel

score(){ rm -rf "e_$1"
  rtg vcfeval -b "$2" -c "$3" -t sdf --squash-ploidy \
      --evaluation-regions regions.bed -o "e_$1" 2>&1 | grep -E 'Threshold|None' | tail -2; }

echo ""
echo "=== ARCS SNV ==="
score arcs_snv t_snv.vcf.gz c_snv.vcf.gz
echo "=== ARCS INDEL ==="
score arcs_ind t_ind.vcf.gz c_ind.vcf.gz

# DiscoSNP++ scoring (if vcf has content)
if [ -s disco_coherent.vcf ] && grep -q '^[^#]' disco_coherent.vcf 2>/dev/null; then
    # Lift DiscoSNP++ output: it uses its own coordinates, need BWA placement too
    bwa mem -t4 "$REF" disco*coherent.fa 2>/dev/null > disco_c2r.sam 2>/dev/null || true
    python3 "$SC/lift_vcf.py" disco_coherent.vcf disco_c2r.sam "$REF" $CHROM disco_lifted.vcf disco*coherent.fa 2>/dev/null || cp disco_coherent.vcf disco_lifted.vcf
    prep disco_lifted.vcf d_snv snv
    prep disco_lifted.vcf d_ind indel
    echo "=== DiscoSNP++ SNV ==="
    score disco_snv t_snv.vcf.gz d_snv.vcf.gz
    echo "=== DiscoSNP++ INDEL ==="
    score disco_ind t_ind.vcf.gz d_ind.vcf.gz
fi

echo ""
echo "=== DONE: chr${CHROM}:${RLO}-${RHI} ==="
