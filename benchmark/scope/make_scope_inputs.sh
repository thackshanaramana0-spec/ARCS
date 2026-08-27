#!/usr/bin/env bash
# Build "scope-matched" FASTQ from a real one: same sequences, but names
# replaced by @1..@N and every quality set to 'I'.
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
set -eu
[ $# -eq 2 ] || { echo "usage: $0 <in.fq> <out.fq>" >&2; exit 1; }
awk 'NR%4==1{n++; print "@"n}
     NR%4==2{print; L=length($0)}
     NR%4==3{print "+"}
     NR%4==0{q=""; for(i=0;i<L;i++) q=q "I"; print q}' "$1" > "$2"
