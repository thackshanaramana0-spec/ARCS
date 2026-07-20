#!/bin/bash
IN=/mnt/c/Temp/arcs_test/d1.fastq
echo "--- sizes (input / decode) ---"; wc -c < "$IN"; wc -c < /tmp/lc_base.dec
echo "--- first diff ---"; cmp "$IN" /tmp/lc_base.dec | head
echo "--- input line1 ---"; head -1 "$IN" | cat -A
echo "--- decode line1 ---"; head -1 /tmp/lc_base.dec | cat -A
echo "--- input tail 40B ---"; tail -c 40 "$IN" | xxd | tail -3
echo "--- decode tail 40B ---"; tail -c 40 /tmp/lc_base.dec | xxd | tail -3
