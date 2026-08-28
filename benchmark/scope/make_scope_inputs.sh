#!/usr/bin/env bash
# Build "scope-matched" FASTQ from a real one: take the first N reads, then
# replace names with @1..@N and set every quality to 'I'.
#
# WHY: PgRC2 stores DNA sequences only -- no names, no qualities. Comparing its
# whole archive against ARCS's whole archive is therefore not a comparison of
# the same thing. Flattening names and qualities makes them cost ~nothing for
# every tool, so what each archive measures is sequence coding.
#
# READ THE CAVEAT IN ../PGRC2_COMPARISON.md BEFORE USING THESE NUMBERS: this
# transformation is NOT neutral for ARCS. It destroys the numeric name field
# that ARCS derives read order from for free, forcing it to store an explicit
# permutation (2.57 MB on 1M yeast reads) that it does not pay for on real
# FASTQ. These inputs are for the sequence-stream comparison, not for a
# whole-archive one.
#
# THE THREE INPUTS BEHIND EVERY NUMBER IN ../PGRC2_COMPARISON.md AND
# ../reimpl/README.md. Subsetting is a plain head -- no sampling, no seed, so
# these are exactly reproducible:
#
#   yeast_sub.fq  1000000 reads  <- /data/fastq/DRR976266_1.fq
#   ecoli_sub.fq  1000000 reads  <- /data/fastq/SRR2584863_1.fq
#   pf_sub.fq      500000 reads  <- /data/fastq/SRR37283774_1.fq
#
#   bash make_scope_inputs.sh /data/fastq/DRR976266_1.fq   scope/yeast_sub.fq 1000000
#   bash make_scope_inputs.sh /data/fastq/SRR2584863_1.fq  scope/ecoli_sub.fq 1000000
#   bash make_scope_inputs.sh /data/fastq/SRR37283774_1.fq scope/pf_sub.fq     500000
#
# Verified 2026-08-28: the retained /tmp/scope copies match these commands to
# the read (sequence of records 1 and 2 identical to the source in all three).
set -eu
[ $# -ge 2 ] || { echo "usage: $0 <in.fq> <out.fq> [n_reads]" >&2; exit 1; }
IN=$1; OUT=$2; NREADS=${3:-0}

# Whole file when n_reads is omitted or 0, else the first n_reads records.
if [ "$NREADS" -gt 0 ] 2>/dev/null; then
    SRC=$(mktemp); trap 'rm -f "$SRC"' EXIT
    head -n $((NREADS * 4)) "$IN" > "$SRC"
    got=$(( $(wc -l < "$SRC") / 4 ))
    [ "$got" -eq "$NREADS" ] || { echo "FAIL: wanted $NREADS reads, got $got" >&2; exit 1; }
else
    SRC=$IN
fi

awk 'NR%4==1{n++; print "@"n}
     NR%4==2{print; L=length($0)}
     NR%4==3{print "+"}
     NR%4==0{q=""; for(i=0;i<L;i++) q=q "I"; print q}' "$SRC" > "$OUT"
