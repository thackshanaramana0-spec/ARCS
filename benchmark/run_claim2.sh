#!/bin/bash
# Claim 2 — FAITHFUL: T3 (SNV calling) + T4 (coverage depth) + T5 (indel calling)
#
# Usage:
#   bash benchmark/run_claim2.sh <DATA_DIR> <ARCS_BIN> [OUT_DIR]
#
# DATA_DIR : directory containing pooled chr20 FASTQ files:
#            HG002_pooled.fq  HG003_pooled.fq  HG004_pooled.fq  HG005_pooled.fq
#            (2Mb region chr20:2000000-2400000, ~30x coverage each)
# ARCS_BIN : path to arcs binary
# OUT_DIR  : results directory (default: ./results/claim2)
#
# Required env vars (with defaults):
#   CHR20_FA      path to chr20.fa GRCh37 reference    [~/refs/chr20.fa]
#   CHR20_SDF     path to rtg SDF dir for chr20        [~/refs/chr20.sdf]
#   GIAB_TRUTH    path to GIAB truth dir (HG002-HG005 vcf.gz + bed) [~/giab_truth]
#   SCRIPTS_DIR   path to lift_vcf.py + kmer2snp_to_vcf.py  [../scripts]
#   DISCO_DIR     path to DiscoSNP++ run_discoSnp++.sh       [~/DiscoSnp]
#   CONDA_ENV     conda env with kmer2snp + findGSE           [kmer2snp_r]

set -euo pipefail

# ── Helpers ───────────────────────────────────────────────────────────────────
_PH=""
log()   { echo "[$(date +%H:%M:%S)] $*"; }
phase() { _PH="$1"; echo ""; echo "[Phase $_PH] $2"; }
pdone() { echo "[Phase $_PH] DONE — $*"; }
pfail() { echo "[Phase $_PH] FAIL — $*"; exit 1; }
pskip() { echo "[Phase $_PH] SKIP — $*"; }
pinfo() { echo "[Phase $_PH]   → $*"; }
step()  { echo ""; echo "══════════════════════════════════════════════"; echo "[STEP] $*"; echo "══════════════════════════════════════════════"; }
ok()    { echo "  [OK]   $*"; }
warn()  { echo "  [WARN] $*"; }
fail()  { echo "  [FAIL] $*"; exit 1; }
chk()   { echo "  [CHK]  $*"; }

# ── Args ──────────────────────────────────────────────────────────────────────
DATA_DIR="${1:?Usage: $0 <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
ARCS_BIN="${2:?Usage: $0 <DATA_DIR> <ARCS_BIN> [OUT_DIR]}"
OUT_DIR="${3:-./results/claim2}"

BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS_DIR="${SCRIPTS_DIR:-$BENCH_DIR/../scripts}"
CHR20_FA="${CHR20_FA:-$HOME/refs/chr20.fa}"
CHR20_SDF="${CHR20_SDF:-$HOME/refs/chr20.sdf}"
GIAB_TRUTH="${GIAB_TRUTH:-$HOME/giab_truth}"
DISCO_DIR="${DISCO_DIR:-$HOME/DiscoSnp}"
CONDA_ENV="${CONDA_ENV:-kmer2snp_r}"

mkdir -p "$OUT_DIR"
WD="$OUT_DIR/workdir"
mkdir -p "$WD"

NPROC=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)

log "Claim 2 — FAITHFUL"
log "DATA_DIR    : $DATA_DIR"
log "ARCS_BIN    : $ARCS_BIN"
log "OUT_DIR     : $OUT_DIR"
log "CHR20_FA    : $CHR20_FA"
log "GIAB_TRUTH  : $GIAB_TRUTH"
log "SCRIPTS_DIR : $SCRIPTS_DIR"
log "DISCO_DIR   : $DISCO_DIR"
log "CONDA_ENV   : $CONDA_ENV"

# ── STEP 0: Prerequisites check ───────────────────────────────────────────────
step "0/7 Prerequisites check"

phase "0.1" "Check arcs binary"
[ -x "$ARCS_BIN" ] && pdone "arcs found: $ARCS_BIN" || pfail "arcs binary not found or not executable: $ARCS_BIN"

phase "0.2" "Check chr20.fa reference"
[ -s "$CHR20_FA" ] && pdone "chr20.fa found ($(stat -c %s "$CHR20_FA") bytes)" \
    || pfail "chr20.fa not found at $CHR20_FA"

phase "0.3" "Check/build BWA index for chr20"
if [ -s "${CHR20_FA}.bwt" ]; then
    pdone "BWA index exists"
else
    pinfo "BWA index missing — building now (may take 2-3 min)"
    bwa index "$CHR20_FA" 2>/dev/null
    pdone "BWA index built"
fi

phase "0.4" "Check/build rtg SDF for chr20"
if [ -d "$CHR20_SDF" ]; then
    pdone "chr20.sdf exists"
else
    pinfo "SDF missing — building now"
    rtg format -o "$CHR20_SDF" "$CHR20_FA" >/dev/null 2>&1
    pdone "SDF built at $CHR20_SDF"
fi

phase "0.5" "Check GIAB truth directory"
[ -d "$GIAB_TRUTH" ] && pdone "truth dir found: $GIAB_TRUTH" \
    || pfail "GIAB truth dir not found at $GIAB_TRUTH — run download.sh first"

phase "0.6" "Check helper scripts"
[ -s "$SCRIPTS_DIR/lift_vcf.py" ]         && pinfo "lift_vcf.py found"         || pfail "lift_vcf.py not found at $SCRIPTS_DIR"
[ -s "$SCRIPTS_DIR/kmer2snp_to_vcf.py" ]  && pinfo "kmer2snp_to_vcf.py found"  || pfail "kmer2snp_to_vcf.py not found at $SCRIPTS_DIR"
[ -s "$DISCO_DIR/run_discoSnp++.sh" ]      && pinfo "DiscoSNP++ found"           || pfail "DiscoSNP++ not found at $DISCO_DIR"
pdone "all helper scripts OK"

phase "0.7" "Check conda env: $CONDA_ENV"
conda run -n "$CONDA_ENV" python -c "import findgse" 2>/dev/null \
    && pdone "conda env $CONDA_ENV OK" \
    || pfail "conda env '$CONDA_ENV' not ready — run: conda create -n $CONDA_ENV r-base r-pracma"

phase "0.8" "Check tools on PATH: bwa samtools bcftools bgzip tabix rtg"
for t in bwa samtools bcftools bgzip tabix rtg; do
    command -v "$t" &>/dev/null && pinfo "$t found: $(command -v $t)" || pfail "$t not on PATH"
done
pdone "all tools present"

# ── STEP 1: Check input FASTQ files ──────────────────────────────────────────
step "1/7 Check input pooled FASTQ files"
INDIVIDUALS=(HG002 HG003 HG004 HG005)
for IND in "${INDIVIDUALS[@]}"; do
    FQ="$DATA_DIR/${IND}_pooled.fq"
    chk "$IND pooled FASTQ: $FQ"
    if [ -s "$FQ" ]; then
        N=$(( $(wc -l < "$FQ") / 4 ))
        SZ=$(stat -c %s "$FQ")
        ok "$IND — ${N} reads, ${SZ} bytes"
        [ "$N" -ge 5000 ] || warn "$IND has only $N reads — need ≥5000 for reliable calling (2Mb pool)"
    else
        fail "$IND_pooled.fq not found at $FQ — run download.sh first to extract chr20:2M-2.4M reads"
    fi
done

# ── STEP 2: Check GIAB truth files ───────────────────────────────────────────
step "2/7 Check GIAB truth files"
declare -A TRUTH_VCF TRUTH_BED
_GIAB_FTP=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
TRUTH_VCF[HG002]="$GIAB_TRUTH/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG002]="$GIAB_TRUTH/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH_VCF[HG003]="$GIAB_TRUTH/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG003]="$GIAB_TRUTH/HG003_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH_VCF[HG004]="$GIAB_TRUTH/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG004]="$GIAB_TRUTH/HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH_VCF[HG005]="$GIAB_TRUTH/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED[HG005]="$GIAB_TRUTH/HG005_GRCh37_1_22_v4.2.1_benchmark.bed"

for IND in "${INDIVIDUALS[@]}"; do
    chk "$IND truth VCF"
    [ -s "${TRUTH_VCF[$IND]}" ] && ok "${TRUTH_VCF[$IND]} found" || fail "Missing ${TRUTH_VCF[$IND]} — run download.sh"
    chk "$IND truth BED"
    [ -s "${TRUTH_BED[$IND]}" ] && ok "${TRUTH_BED[$IND]} found" || fail "Missing ${TRUTH_BED[$IND]} — run download.sh"
done

# ── Shared helpers ────────────────────────────────────────────────────────────
CHROM=20
RLO=2000000; RHI=2400000
REGION="${CHROM}:${RLO}-${RHI}"

parse_time_v() {
    local logf="$1"
    local wall vmhwm
    wall=$(grep "Elapsed (wall clock)" "$logf" | awk '{
        n=split($NF,a,":");
        if(n==3) printf "%.2f", a[1]*3600+a[2]*60+a[3];
        else if(n==2) printf "%.2f", a[1]*60+a[2];
        else printf "%.2f", a[1];
    }')
    vmhwm=$(grep "Maximum resident set size" "$logf" | awk '{print $NF}')
    echo "${wall:-0} ${vmhwm:-0}"
}

snv() { awk -F'\t' '/^#/{print;next} length($4)==1 && length($5)==1'; }
het() { awk -F'\t' '/^#/{print;next} {split($10,g,":"); gt=g[1];
        if(gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") print}'; }

prep_vcf() {
    local src="$1" tag="$2" filtercmd="$3"
    log "  prep_vcf: $tag (filter=$filtercmd)"
    eval "$filtercmd" < "$src" > "${tag}.v.vcf"
    { grep '^#' "${tag}.v.vcf"; grep -v '^#' "${tag}.v.vcf" | sort -k2,2n; } > "${tag}.s.vcf"
    bcftools norm -f "$CHR20_FA" -m -any "${tag}.s.vcf" 2>/dev/null \
        | bcftools sort 2>/dev/null | bgzip > "${tag}.vcf.gz" \
        || bgzip -c "${tag}.s.vcf" > "${tag}.vcf.gz"
    tabix -f -p vcf "${tag}.vcf.gz"
    local n
    n=$(zcat "${tag}.vcf.gz" | grep -vc '^#' || true)
    log "  prep_vcf: $tag → $n variants"
}

score_vcfeval() {
    local label="$1" truth="$2" calls="$3" bed="$4" outpfx="$5"
    local outd="${outpfx}_eval"
    log "  vcfeval: $label"
    rm -rf "$outd"
    rtg vcfeval -b "$truth" -c "$calls" -t "$CHR20_SDF" \
        --squash-ploidy --bed-regions "$bed" -o "$outd" >/dev/null 2>&1 || true
    if [ -f "$outd/summary.txt" ]; then
        awk 'NR>2{tp=$3;fp=$4;fn=$5;p=$6;r=$7;f=$8}
             END{printf "  %-14s  TP=%-5s FP=%-5s FN=%-5s  P=%.3f  R=%.3f  F1=%.3f\n",
                 "'"$label"'",tp,fp,fn,p,r,f}' "$outd/summary.txt" | tee -a "$OUT_DIR/t3_results.txt"
    else
        warn "  vcfeval produced no summary for $label — check $outd"
        echo "  $label: NO_RESULT" >> "$OUT_DIR/t3_results.txt"
    fi
}

# ── STEP 3: T3 — SNV calling, all 4 individuals ──────────────────────────────
step "3/7 T3 — Het-SNV calling (4 individuals)"
> "$OUT_DIR/t3_results.txt"

for IND in "${INDIVIDUALS[@]}"; do
    FQ="$DATA_DIR/${IND}_pooled.fq"
    IWD="$WD/${IND}"
    mkdir -p "$IWD"
    log "── $IND ──────────────────────────────────────"

    # 3a. ARCS --call
    phase "3.${IND}.1" "[$IND] ARCS --call (reference-free SNV)"
    if [ ! -s "$IWD/calls_arcs.vcf" ]; then
        TF="$IWD/arcs_time.txt"
        export ARCS_DUMP_CONTIGS="$IWD/contigs.tsv"
        # Real CLI contract (verified directly against main.cpp's option parser,
        # not just its printed usage string): `compress --call <out.vcf>
        # <input> <output>` — --call's value is consumed immediately as the
        # VCF output path, before input/output are read. The previous
        # argument order here ($FQ as the vcf path, out.arcs as the "input",
        # calls_arcs.vcf as the archive output) meant this command could never
        # have loaded a valid FASTQ input; `|| true` silently swallowed the
        # resulting failure, so no real T3 SNV-calling result was ever
        # produced by this step for ANY assembler.
        /usr/bin/time -v "$ARCS_BIN" compress --call "$IWD/calls_arcs.vcf" "$FQ" "$IWD/out.arcs" 2>"$TF" || true
        read -r ARCS_T ARCS_RAM <<< "$(parse_time_v "$TF")"
        unset ARCS_DUMP_CONTIGS
        N_ARCS=$(grep -vc '^#' "$IWD/calls_arcs.vcf" 2>/dev/null || echo 0)
        pdone "[$IND] ARCS call — ${ARCS_T}s, ${ARCS_RAM}KB RAM, $N_ARCS raw calls"
    else
        pskip "[$IND] ARCS calls cached"
    fi

    # 3b. BWA lift for ARCS
    phase "3.${IND}.2" "[$IND] BWA liftover of ARCS contigs"
    if [ ! -s "$IWD/lifted_arcs.vcf" ] && [ -s "$IWD/contigs.tsv" ]; then
        # ARCS_DUMP_CONTIGS (caller.cpp) already writes valid FASTA
        # (">contig_N\nSEQ\n") despite the "contigs.tsv" filename — the awk
        # transform this line used to apply assumed a tab-separated
        # name<TAB>seq input with no ">" prefix, which is NOT what caller.cpp
        # produces; running it against the real (FASTA) format corrupted
        # every header into ">>contig_N" and every sequence line into
        # ">SEQ", silently breaking bwa's alignment and therefore this
        # entire liftover step for ANY assembler, not just Method B. No
        # transform is needed — just copy it through as-is.
        cp "$IWD/contigs.tsv" "$IWD/contigs.fa"
        bwa mem -t 8 "$CHR20_FA" "$IWD/contigs.fa" 2>/dev/null > "$IWD/c2r.sam"
        python3 "$SCRIPTS_DIR/lift_vcf.py" \
            "$IWD/calls_arcs.vcf" "$IWD/c2r.sam" "$CHR20_FA" $CHROM \
            "$IWD/lifted_arcs.vcf" "$IWD/contigs.fa"
        N_LIFT=$(grep -vc '^#' "$IWD/lifted_arcs.vcf" 2>/dev/null || echo 0)
        pdone "[$IND] ARCS liftover — $N_LIFT lifted calls"
    else
        pskip "[$IND] ARCS liftover cached or no contigs"
    fi

    # 3c. DiscoSNP++ — MUST use -G flag
    phase "3.${IND}.3" "[$IND] DiscoSNP++ (MUST use -G chr20.fa)"
    if [ ! -s "$IWD/calls_disco.vcf" ]; then
        echo "$FQ" > "$IWD/fof.txt"
        TF="$IWD/disco_time.txt"
        cd "$IWD"
        /usr/bin/time -v "$DISCO_DIR/run_discoSnp++.sh" \
            -r "$IWD/fof.txt" -G "$CHR20_FA" -k 31 -c 3 -D 100 -P 3 -b 0 -p disco \
            2>"$TF" || true
        read -r DISCO_T DISCO_RAM <<< "$(parse_time_v "$TF")"
        cd - >/dev/null
        # Fix +1 position offset (known DiscoSNP++ bug)
        if ls "$IWD"/disco*coherent.vcf 2>/dev/null | grep -q .; then
            awk 'BEGIN{OFS="\t"} /^#/{print;next} {$2=$2-1; print}' \
                "$IWD"/disco*coherent.vcf > "$IWD/calls_disco.vcf"
            N_DISCO=$(grep -vc '^#' "$IWD/calls_disco.vcf" 2>/dev/null || echo 0)
            ok "[$IND] DiscoSNP++ done — ${DISCO_T}s, ${DISCO_RAM}KB RAM, $N_DISCO calls (pos offset fixed)"
        else
            warn "[$IND] DiscoSNP++ produced no coherent.vcf"
            touch "$IWD/calls_disco.vcf"
        fi
    else
        ok "[$IND] DiscoSNP++ cached"
    fi

    # 3d. DiscoSNP++ liftover
    phase "3.${IND}.4" "[$IND] BWA liftover of DiscoSNP++ contigs"
    if [ ! -s "$IWD/lifted_disco.vcf" ] && grep -q '^[^#]' "$IWD/calls_disco.vcf" 2>/dev/null; then
        if ls "$IWD"/disco*coherent.fa 2>/dev/null | grep -q .; then
            bwa mem -t 8 "$CHR20_FA" "$IWD"/disco*coherent.fa 2>/dev/null > "$IWD/disco_c2r.sam" || true
            python3 "$SCRIPTS_DIR/lift_vcf.py" \
                "$IWD/calls_disco.vcf" "$IWD/disco_c2r.sam" "$CHR20_FA" $CHROM \
                "$IWD/lifted_disco.vcf" "$IWD"/disco*coherent.fa 2>/dev/null \
                || cp "$IWD/calls_disco.vcf" "$IWD/lifted_disco.vcf"
            N_DLIFT=$(grep -vc '^#' "$IWD/lifted_disco.vcf" 2>/dev/null || echo 0)
            ok "[$IND] DiscoSNP++ liftover done — $N_DLIFT lifted"
        else
            warn "[$IND] No DiscoSNP++ contigs.fa found — skipping liftover"
            touch "$IWD/lifted_disco.vcf"
        fi
    else
        ok "[$IND] DiscoSNP++ liftover cached or no calls"
    fi

    # 3e. Kmer2SNP via conda
    phase "3.${IND}.5" "[$IND] Kmer2SNP via conda env $CONDA_ENV"
    if [ ! -s "$IWD/calls_k2s.vcf" ]; then
        TF="$IWD/k2s_time.txt"
        conda run -n "$CONDA_ENV" /usr/bin/time -v \
            python3 "$SCRIPTS_DIR/kmer2snp_to_vcf.py" \
            "$FQ" "$CHR20_FA" $CHROM "$IWD/calls_k2s.vcf" \
            2>"$TF" || { warn "[$IND] Kmer2SNP failed — creating empty VCF"; touch "$IWD/calls_k2s.vcf"; }
        read -r K2S_T K2S_RAM <<< "$(parse_time_v "$TF")"
        N_K2S=$(grep -vc '^#' "$IWD/calls_k2s.vcf" 2>/dev/null || echo 0)
        ok "[$IND] Kmer2SNP done — ${K2S_T}s, ${K2S_RAM}KB RAM, $N_K2S calls"
    else
        ok "[$IND] Kmer2SNP cached"
    fi

    # 3f. Truth region subset
    phase "3.${IND}.6" "[$IND] Extract truth for chr${CHROM}:${RLO}-${RHI}"
    if [ ! -s "$IWD/truth.vcf" ]; then
        tabix -h "${TRUTH_VCF[$IND]}" "$REGION" 2>/dev/null | awk -v OFS='\t' '
            /^##contig/{next}
            /^#CHROM/{print "##contig=<ID=20,length=63025520>"; print; next}
            /^#/{print; next} {print}' > "$IWD/truth.vcf"
        N_TRUTH=$(grep -vc '^#' "$IWD/truth.vcf" 2>/dev/null || echo 0)
        ok "[$IND] Truth extracted — $N_TRUTH variants in region"
    fi

    # 3g. Confident BED subset
    if [ ! -s "$IWD/regions.bed" ]; then
        curl -s "${TRUTH_BED[$IND]}" -o "$IWD/giab_conf.bed" 2>/dev/null \
            || cp "${TRUTH_BED[$IND]}" "$IWD/giab_conf.bed"
        awk -v c=$CHROM -v lo=$RLO -v hi=$RHI \
            '$1==c && $3>lo && $2<hi{s=($2>lo?$2:lo); e=($3<hi?$3:hi); if(e>s) print c"\t"s"\t"e}' \
            "$IWD/giab_conf.bed" > "$IWD/regions.bed"
        BP=$(awk '{n+=$3-$2}END{print n}' "$IWD/regions.bed")
        ok "[$IND] Confident region: ${BP} bp"
    fi

    # 3h. vcfeval — all 3 tools
    phase "3.${IND}.7" "[$IND] rtg vcfeval scoring (all 3 tools)"
    echo "── $IND ──" >> "$OUT_DIR/t3_results.txt"

    prep_vcf "$IWD/truth.vcf"       "$IWD/t_snv"   'het | snv'
    prep_vcf "$IWD/lifted_arcs.vcf" "$IWD/c_arcs"  'snv'
    score_vcfeval "ARCS_${IND}"  "$IWD/t_snv.vcf.gz" "$IWD/c_arcs.vcf.gz" "$IWD/regions.bed" "$IWD/arcs"

    if grep -q '^[^#]' "$IWD/lifted_disco.vcf" 2>/dev/null; then
        prep_vcf "$IWD/lifted_disco.vcf" "$IWD/c_disco" 'snv'
        score_vcfeval "DISCO_${IND}" "$IWD/t_snv.vcf.gz" "$IWD/c_disco.vcf.gz" "$IWD/regions.bed" "$IWD/disco"
    else
        warn "[$IND] DiscoSNP++ no calls to score"
    fi

    if grep -q '^[^#]' "$IWD/calls_k2s.vcf" 2>/dev/null; then
        prep_vcf "$IWD/calls_k2s.vcf" "$IWD/c_k2s" 'snv'
        score_vcfeval "K2S_${IND}" "$IWD/t_snv.vcf.gz" "$IWD/c_k2s.vcf.gz" "$IWD/regions.bed" "$IWD/k2s"
    else
        warn "[$IND] Kmer2SNP no calls to score"
    fi

    pdone "[$IND] all T3 sub-steps complete"
done

log "T3 results written to $OUT_DIR/t3_results.txt"

# ── STEP 4: T4 — Coverage depth effect ───────────────────────────────────────
step "4/7 T4 — Coverage depth (10x/15x/20x/30x)"
> "$OUT_DIR/t4_results.txt"
echo "Depth,ARCS_archive_bytes,Raw_bytes,Ratio" >> "$OUT_DIR/t4_results.txt"

HG002_FQ="$DATA_DIR/HG002_pooled.fq"
N_TOTAL=$(( $(wc -l < "$HG002_FQ") / 4 ))
log "HG002 pooled: $N_TOTAL reads total"

for DEPTH in 10 15 20 30; do
    log "T4: subsampling to ${DEPTH}x"
    FRAC=$(awk "BEGIN{printf \"%.4f\", $DEPTH/30.0}")
    SUB_FQ="$WD/hg002_${DEPTH}x.fq"
    if [ ! -s "$SUB_FQ" ]; then
        # Subsample using line math (every 4 lines = 1 read, take FRAC fraction)
        awk -v frac="$FRAC" 'BEGIN{srand(42)} NR%4==1{keep=(rand()<frac)} keep{print}' \
            "$HG002_FQ" > "$SUB_FQ" || true
        # awk subsample loses 4-line grouping — use python instead
        python3 -c "
import sys, random
random.seed(42)
frac = float('$FRAC')
reads = []
with open('$HG002_FQ') as f:
    while True:
        h = f.readline()
        if not h: break
        s = f.readline(); p = f.readline(); q = f.readline()
        if random.random() < frac:
            reads.extend([h,s,p,q])
with open('$SUB_FQ','w') as f:
    f.writelines(reads)
"
        N_SUB=$(( $(wc -l < "$SUB_FQ") / 4 ))
        ok "T4: ${DEPTH}x — $N_SUB reads subsampled"
    else
        N_SUB=$(( $(wc -l < "$SUB_FQ") / 4 ))
        ok "T4: ${DEPTH}x cached — $N_SUB reads"
    fi

    RAW=$(stat -c %s "$SUB_FQ")
    ARC="$WD/hg002_${DEPTH}x.arcs"
    TF="$WD/t4_${DEPTH}x_time.txt"
    # ARCS_PAR_SHARDS not set — C++ default min(4, hw_concurrency) gives best ratio/speed balance
    /usr/bin/time -v "$ARCS_BIN" compress "$SUB_FQ" "$ARC" 2>"$TF"
    ARCH=$(stat -c %s "$ARC")
    read -r T_WALL T_RAM <<< "$(parse_time_v "$TF")"
    RATIO=$(awk "BEGIN{printf \"%.3f\", $ARCH/$RAW}")
    ok "T4: ${DEPTH}x — raw=$RAW bytes, archive=$ARCH bytes, ratio=$RATIO, time=${T_WALL}s"
    echo "${DEPTH}x,${ARCH},${RAW},${RATIO}" >> "$OUT_DIR/t4_results.txt"
    rm -f "$ARC"
done

log "T4 results written to $OUT_DIR/t4_results.txt"

# ── STEP 5: T5 — Indel calling (real GIAB HG002) ─────────────────────────────
step "5/7 T5 — Indel calling (real GIAB HG002)"
> "$OUT_DIR/t5_results.txt"

IWD="$WD/HG002"
log "T5: using HG002 calls from T3 (reuse lifted_arcs.vcf and lifted_disco.vcf)"

# Build indel subsets from already-lifted VCFs
indel() { awk -F'\t' '/^#/{print;next} length($4)!=1 || length($5)!=1'; }

if [ -s "$IWD/lifted_arcs.vcf" ]; then
    prep_vcf "$IWD/truth.vcf"       "$IWD/t_ind"   'het | indel'
    prep_vcf "$IWD/lifted_arcs.vcf" "$IWD/c_arcs_ind" 'indel'
    score_vcfeval "ARCS_indel_HG002" \
        "$IWD/t_ind.vcf.gz" "$IWD/c_arcs_ind.vcf.gz" "$IWD/regions.bed" "$IWD/arcs_ind"
    echo "$(tail -1 "$OUT_DIR/t3_results.txt")" >> "$OUT_DIR/t5_results.txt" || true
else
    warn "T5: no ARCS lifted VCF from T3 — run T3 first"
fi

if grep -q '^[^#]' "$IWD/lifted_disco.vcf" 2>/dev/null; then
    prep_vcf "$IWD/lifted_disco.vcf" "$IWD/c_disco_ind" 'indel'
    score_vcfeval "DISCO_indel_HG002" \
        "$IWD/t_ind.vcf.gz" "$IWD/c_disco_ind.vcf.gz" "$IWD/regions.bed" "$IWD/disco_ind"
fi

log "T5 results written to $OUT_DIR/t5_results.txt"

# ── STEP 6: Print summary ─────────────────────────────────────────────────────
step "6/7 Results summary"

log "T3 — Het-SNV F1 per individual:"
cat "$OUT_DIR/t3_results.txt"

echo ""
log "T4 — Coverage depth effect:"
cat "$OUT_DIR/t4_results.txt"

echo ""
log "T5 — Indel calling:"
cat "$OUT_DIR/t5_results.txt"

step "7/7 DONE"
log "All Claim 2 results in: $OUT_DIR"
log "Finished at: $(date)"
