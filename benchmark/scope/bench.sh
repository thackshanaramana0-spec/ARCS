set -u
ARCS_NEW=/root/arcs-clean/build/arcs
ARCS_OLD=/tmp/arcs_old/build/arcs
PGRC=/root/arcs-clean/method_c/build/PgRC
W=/tmp/scope/work
mkdir -p $W
OUT=/tmp/scope/results.csv
echo "dataset,tool,bytes,comp_s,decomp_s,peak_rss_kb,lossless" > $OUT

# wall+RSS around one command; prints "secs rsskb"
run() { /usr/bin/time -f "%e %M" -o $W/t.txt "$@" > $W/log 2>&1; echo $? > $W/rc; cat $W/t.txt; }

seqs() { awk 'NR%4==2' "$1"; }

for d in yeast_sub ecoli_sub pf_sub; do
  IN=/tmp/scope/$d.fq
  seqs $IN > $W/ref.seq

  # ── ARCS new (zero flags) ────────────────────────────────────────────────
  for pair in "ARCS_new:$ARCS_NEW" "ARCS_old:$ARCS_OLD"; do
    NAME=${pair%%:*}; BIN=${pair#*:}
    rm -f $W/a.arcs $W/a.fq
    read C R <<< "$(run $BIN compress $IN $W/a.arcs)"; rc=$(cat $W/rc)
    if [ "$rc" != 0 ]; then echo "$d,$NAME,,,,,COMPRESS_FAIL" >> $OUT; continue; fi
    read D R2 <<< "$(run $BIN decompress $W/a.arcs $W/a.fq)"; rc=$(cat $W/rc)
    if [ "$rc" != 0 ]; then echo "$d,$NAME,$(stat -c%s $W/a.arcs),$C,,$R,DECOMPRESS_FAIL" >> $OUT; continue; fi
    if seqs $W/a.fq > $W/d.seq && cmp -s $W/ref.seq $W/d.seq; then L=LOSSLESS; else L=LOSSY; fi
    echo "$d,$NAME,$(stat -c%s $W/a.arcs),$C,$D,$R,$L" >> $OUT
    rm -f $W/a.arcs $W/a.fq $W/d.seq
  done

  # ── PgRC2 (order-preserving, 12 threads) ─────────────────────────────────
  rm -f $W/p.pgrc $W/p.pgrc_out
  read C R <<< "$(run $PGRC -o -t 12 -i $IN $W/p.pgrc)"; rc=$(cat $W/rc)
  if [ "$rc" != 0 ]; then echo "$d,PgRC2,,,,,COMPRESS_FAIL" >> $OUT; else
    read D R2 <<< "$(run $PGRC -d -t 12 $W/p.pgrc)"; rc=$(cat $W/rc)
    if [ "$rc" != 0 ]; then echo "$d,PgRC2,$(stat -c%s $W/p.pgrc),$C,,$R,DECOMPRESS_FAIL" >> $OUT; else
      if cmp -s $W/ref.seq $W/p.pgrc_out; then L=LOSSLESS; else L=LOSSY; fi
      echo "$d,PgRC2,$(stat -c%s $W/p.pgrc),$C,$D,$R,$L" >> $OUT
    fi
  fi
  rm -f $W/p.pgrc $W/p.pgrc_out

  # ── SPRING (sequences only: --no-quality --no-ids) ───────────────────────
  rm -f $W/s.spring $W/s.fq
  read C R <<< "$(run spring -c -i $IN -o $W/s.spring -t 12 -g --no-quality --no-ids)"; rc=$(cat $W/rc)
  if [ "$rc" != 0 ]; then echo "$d,SPRING,,,,,COMPRESS_FAIL" >> $OUT; else
    read D R2 <<< "$(run spring -d -i $W/s.spring -o $W/s.fq -t 12 -g)"; rc=$(cat $W/rc)
    if [ "$rc" != 0 ]; then echo "$d,SPRING,$(stat -c%s $W/s.spring),$C,,$R,DECOMPRESS_FAIL" >> $OUT; else
      if seqs $W/s.fq > $W/d.seq && cmp -s $W/ref.seq $W/d.seq; then L=LOSSLESS; else L=LOSSY; fi
      echo "$d,SPRING,$(stat -c%s $W/s.spring),$C,$D,$R,$L" >> $OUT
    fi
  fi
  rm -f $W/s.spring $W/s.fq $W/d.seq

  # ── Genozip ──────────────────────────────────────────────────────────────
  rm -f $W/g.genozip $W/g.fq
  read C R <<< "$(run genozip --force -o $W/g.genozip $IN)"; rc=$(cat $W/rc)
  if [ "$rc" != 0 ]; then echo "$d,Genozip,,,,,COMPRESS_FAIL" >> $OUT; else
    read D R2 <<< "$(run genounzip --force -o $W/g.fq $W/g.genozip)"; rc=$(cat $W/rc)
    if [ "$rc" != 0 ]; then echo "$d,Genozip,$(stat -c%s $W/g.genozip),$C,,$R,DECOMPRESS_FAIL" >> $OUT; else
      if seqs $W/g.fq > $W/d.seq && cmp -s $W/ref.seq $W/d.seq; then L=LOSSLESS; else L=LOSSY; fi
      echo "$d,Genozip,$(stat -c%s $W/g.genozip),$C,$D,$R,$L" >> $OUT
    fi
  fi
  rm -f $W/g.genozip $W/g.fq $W/d.seq $W/ref.seq
  echo "--- $d done ---"
done
echo "=== BENCH DONE ==="
