#!/bin/bash
# Proper losslessness check. The confirm run's "DIFF" hit baseline too → it's the
# known CRLF normalization, not a prune regression. Verify rigorously:
#   (a) is the input CRLF? (explains benign DIFF vs raw input)
#   (b) decode(baseline) == input after \r\n normalization?  (true losslessness)
#   (c) decode(prune-4)  == decode(baseline)?                (prune is byte-safe)
set -u
ARCS=./build_wsl/arcs
IN="${1:?usage}"
norm=/tmp/lc_in_norm.fastq
tr -d '\r' < "$IN" > "$norm"
crlf=$(grep -c $'\r' "$IN" 2>/dev/null || echo 0)
echo "input CRLF lines: $crlf"

ARCSDNA_ORDERS="2,4,8,11,14,18,22,26" $ARCS compress --chain-pg "$IN" /tmp/lc_base.arcs 2>/dev/null
ARCSDNA_ORDERS="2,4,8,11,14,18,22,26" $ARCS decompress /tmp/lc_base.arcs /tmp/lc_base.dec 2>/dev/null
ARCSDNA_ORDERS="2,8,11,14,18,22,26"   $ARCS compress --chain-pg "$IN" /tmp/lc_p4.arcs 2>/dev/null
ARCSDNA_ORDERS="2,8,11,14,18,22,26"   $ARCS decompress /tmp/lc_p4.arcs /tmp/lc_p4.dec 2>/dev/null

echo -n "(b) baseline decode == input(normalized): "; cmp -s "$norm" /tmp/lc_base.dec && echo LOSSLESS || echo DIFF
echo -n "(c) prune-4 decode == baseline decode:    "; cmp -s /tmp/lc_base.dec /tmp/lc_p4.dec && echo IDENTICAL || echo DIFF
echo -n "(c') prune-4 decode == input(normalized): "; cmp -s "$norm" /tmp/lc_p4.dec && echo LOSSLESS || echo DIFF
