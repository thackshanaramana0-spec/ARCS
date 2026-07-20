#!/bin/bash
# Master benchmark: ARCS, SPRING, Genozip, fqzcomp on 8 datasets.
# Data copied to ext4 (/tmp) for fair timing; logs to /mnt/c/Temp/bench8/logs (persist).
set -u
LOGD=/mnt/c/Temp/bench8/logs
WD=/tmp/bench_wd; mkdir -p $WD
ARCS=/tmp/arcs
GZ=/tmp/gz/src/genozip; GC=/tmp/gz/src/genounzip
FQ=/tmp/fqzcomp
export USER="Thackshanaramana B"
cp /mnt/c/Users/Lenovo/OneDrive/Desktop/nu256_upgrade/arcs_cpp/build_wsl/arcs $ARCS; chmod +x $ARCS
cp /mnt/c/Temp/arcs_test/fqzcomp $FQ; chmod +x $FQ 2>/dev/null

declare -A SRC=(
 [DS1]="/mnt/c/Users/Lenovo/OneDrive/Desktop/nu256_upgrade/trial7/data/inputs/orig_DS1.fastq"
 [DS2]="/mnt/c/Users/Lenovo/OneDrive/Desktop/nu256_upgrade/trial7/data/inputs/orig_DS2.fastq"
 [DS4]="/mnt/c/Users/Lenovo/OneDrive/Desktop/nu256_upgrade/trial7/data/inputs/orig_DS4.fastq"
 [DS5]="/mnt/c/Users/Lenovo/OneDrive/Desktop/nu256_upgrade/trial7/data/inputs/orig_DS5.fastq"
 [DS6]="/mnt/c/Users/Lenovo/OneDrive/Desktop/nu256_upgrade/trial7/data/inputs/orig_DS6.fastq"
 [DS7]="/mnt/c/Users/Lenovo/OneDrive/Desktop/nu256_upgrade/trial7/data/inputs/DS7_30x_chr22_5Mbp.fastq"
 [GIAB]="/mnt/c/Temp/arcs_test/giab150.fq"
 [NA12878]="/mnt/c/Temp/arcs_test/na12878_region.fq"
)
ORDER="DS1 DS2 DS4 DS5 DS6 DS7 GIAB NA12878"

# peak VmHWM (KB) of a command run in background
peak(){ "$@" >/dev/null 2>&1 & local P=$!; local m=0 h; while kill -0 $P 2>/dev/null; do h=$(awk "/VmHWM/{print \$2}" /proc/$P/status 2>/dev/null); [ -n "$h" ]&&[ "$h" -gt "$m" ]&&m=$h; sleep 0.03; done; wait $P; echo $m; }
# best-of-2 wall time (s)
btime(){ local b=99999 t; for i in 1 2; do /usr/bin/time -o /tmp/_tt -f "%e" "$@" >/dev/null 2>/dev/null; t=$(cat /tmp/_tt); awk -v t="$t" -v b="$b" "BEGIN{exit !(t+0<b+0)}" && b=$t; done; echo $b; }
losscmp(){ paste - - - - < "$1" | tr -d "\r" | sort > /tmp/_o; paste - - - - < "$2" | tr -d "\r" | sort > /tmp/_d; cmp -s /tmp/_o /tmp/_d && echo LOSSLESS || echo LOSSY; }

for DS in $ORDER; do
  IN=$WD/$DS.fq; cp "${SRC[$DS]}" "$IN"
  raw=$(stat -c %s "$IN")
  echo ">>> $DS  raw=$raw B  $(date +%T)"

  # ---------- ARCS ----------
  A=$WD/$DS.arcs; DEC=$WD/$DS.arcs.dec
  ct=$(btime $ARCS --chain-pg compress "$IN" "$A"); cr=$(peak $ARCS --chain-pg compress "$IN" "$A")
  arch=$(stat -c %s "$A")
  dt=$(btime $ARCS decompress "$A" "$DEC"); dr=$(peak $ARCS decompress "$A" "$DEC")
  ll=$(losscmp "$IN" "$DEC")
  printf "TOOL=ARCS DS=%s raw=%s archive=%s ctime=%s crss_kb=%s dtime=%s drss_kb=%s lossless=%s\n" \
    "$DS" "$raw" "$arch" "$ct" "$cr" "$dt" "$dr" "$ll" | tee $LOGD/${DS}__ARCS.log
  rm -f "$A" "$DEC"

  # ---------- SPRING (default 8 threads) ----------
  A=$WD/$DS.spring; DEC=$WD/$DS.spring.dec
  ct=$(btime spring -c -i "$IN" -o "$A" -w $WD); cr=$(peak spring -c -i "$IN" -o "$A" -w $WD)
  arch=$(stat -c %s "$A")
  dt=$(btime spring -d -i "$A" -o "$DEC" -w $WD); dr=$(peak spring -d -i "$A" -o "$DEC" -w $WD)
  ll=$(losscmp "$IN" "$DEC")
  printf "TOOL=SPRING DS=%s raw=%s archive=%s ctime=%s crss_kb=%s dtime=%s drss_kb=%s lossless=%s\n" \
    "$DS" "$raw" "$arch" "$ct" "$cr" "$dt" "$dr" "$ll" | tee $LOGD/${DS}__SPRING.log
  rm -f "$A" "$DEC"

  # ---------- Genozip ----------
  A=$WD/$DS.genozip; DEC=$WD/$DS.gz.dec
  ct=$(btime $GZ --force -o "$A" "$IN"); cr=$(peak $GZ --force -o "$A" "$IN")
  arch=$(stat -c %s "$A")
  dt=$(btime $GC --force -o "$DEC" "$A"); dr=$(peak $GC --force -o "$DEC" "$A")
  ll=$(losscmp "$IN" "$DEC")
  printf "TOOL=Genozip DS=%s raw=%s archive=%s ctime=%s crss_kb=%s dtime=%s drss_kb=%s lossless=%s\n" \
    "$DS" "$raw" "$arch" "$ct" "$cr" "$dt" "$dr" "$ll" | tee $LOGD/${DS}__Genozip.log
  rm -f "$A" "$DEC"

  # ---------- fqzcomp (NOT byte-lossless: quality quirk) ----------
  A=$WD/$DS.fqz; DEC=$WD/$DS.fqz.dec
  ct=$(btime $FQ "$IN" "$A"); cr=$(peak $FQ "$IN" "$A")
  arch=$(stat -c %s "$A")
  dt=$(btime $FQ -d "$A" "$DEC"); dr=$(peak $FQ -d "$A" "$DEC")
  ll=$(losscmp "$IN" "$DEC")
  printf "TOOL=fqzcomp DS=%s raw=%s archive=%s ctime=%s crss_kb=%s dtime=%s drss_kb=%s lossless=%s\n" \
    "$DS" "$raw" "$arch" "$ct" "$cr" "$dt" "$dr" "$ll" | tee $LOGD/${DS}__fqzcomp.log
  rm -f "$A" "$DEC"

  rm -f "$IN"
done
echo "ALL LINUX TOOLS DONE $(date +%T)" | tee $LOGD/_DONE.flag
