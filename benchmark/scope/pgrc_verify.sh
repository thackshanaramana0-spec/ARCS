set -u
P=/root/arcs-clean/method_c/build/PgRC
W=/tmp/scope/pv; rm -rf $W; mkdir -p $W
OUT=/tmp/scope/pgrc_verify.csv
echo "dataset,mode,run,archive_bytes,reads_out,reads_expected,verdict" > $OUT
for d in yeast_sub ecoli_sub pf_sub; do
  IN=/tmp/scope/$d.fq
  awk 'NR%4==2' $IN > $W/ref            # reference sequences, file order
  sort $W/ref > $W/ref_sorted           # sorted once, reused by every no-order run
  EXP=$(wc -l < $W/ref)
  for run in 1 2 3; do
    # ── order-preserving (-o): must match exactly, read for read ──
    rm -f $W/a.pgrc $W/a.pgrc_out
    $P -o -t 12 -i $IN $W/a.pgrc > /dev/null 2>&1
    $P -d -t 12 $W/a.pgrc > /dev/null 2>&1
    if [ -s $W/a.pgrc_out ]; then
      N=$(wc -l < $W/a.pgrc_out)
      cmp -s $W/ref $W/a.pgrc_out && V=EXACT || V=DIFFERS
      echo "$d,order,$run,$(stat -c%s $W/a.pgrc),$N,$EXP,$V" >> $OUT
    else echo "$d,order,$run,,,,FAILED" >> $OUT; fi
    # ── default (no -o): order is discarded, so compare sorted content ──
    rm -f $W/b.pgrc $W/b.pgrc_out
    $P -t 12 -i $IN $W/b.pgrc > /dev/null 2>&1
    $P -d -t 12 $W/b.pgrc > /dev/null 2>&1
    if [ -s $W/b.pgrc_out ]; then
      N=$(wc -l < $W/b.pgrc_out)
      sort $W/b.pgrc_out > $W/b_sorted
      cmp -s $W/ref_sorted $W/b_sorted && V=CONTENT_OK || V=CONTENT_LOST
      echo "$d,default,$run,$(stat -c%s $W/b.pgrc),$N,$EXP,$V" >> $OUT
    else echo "$d,default,$run,,,,FAILED" >> $OUT; fi
  done
  rm -f $W/ref $W/ref_sorted $W/b_sorted
done
echo "=== PGRC VERIFY DONE ==="
