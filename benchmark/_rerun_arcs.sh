set -u
ROOT=/mnt/c/Users/Lenovo/OneDrive/Desktop/nu256_upgrade
LOGD=$ROOT/arcs-clean/benchmark/logs
WD=/tmp/arcs_rerun2; mkdir -p $WD
ARCS=/tmp/arcs_cur; cp $ROOT/arcs_cpp/build_wsl/arcs $ARCS; chmod +x $ARCS
declare -A SRC=(
 [DS1]="$ROOT/trial7/data/inputs/orig_DS1.fastq"
 [DS2]="$ROOT/trial7/data/inputs/orig_DS2.fastq"
 [DS4]="$ROOT/trial7/data/inputs/orig_DS4.fastq"
 [DS5]="$ROOT/trial7/data/inputs/orig_DS5.fastq"
 [DS6]="$ROOT/trial7/data/inputs/orig_DS6.fastq"
 [DS7]="$ROOT/trial7/data/inputs/DS7_30x_chr22_5Mbp.fastq"
 [GIAB]="/mnt/c/Temp/arcs_test/giab150.fq"
 [NA12878]="/mnt/c/Temp/arcs_test/na12878_region.fq"
)
peak(){ "$@" >/dev/null 2>&1 & local P=$!; local m=0 h; while kill -0 $P 2>/dev/null; do h=$(awk "/VmHWM/{print \$2}" /proc/$P/status 2>/dev/null); [ -n "$h" ]&&[ "$h" -gt "$m" ]&&m=$h; sleep 0.03; done; wait $P; echo $m; }
btime(){ local b=99999 t; for i in 1 2 3; do /usr/bin/time -o /tmp/_tt -f "%e" "$@" >/dev/null 2>/dev/null; t=$(cat /tmp/_tt); awk -v t="$t" -v b="$b" "BEGIN{exit !(t+0<b+0)}" && b=$t; done; echo $b; }
losscmp(){ paste - - - - < "$1" | tr -d "\r" | sort > /tmp/_o; paste - - - - < "$2" | tr -d "\r" | sort > /tmp/_d; cmp -s /tmp/_o /tmp/_d && echo LOSSLESS || echo LOSSY; }
for DS in DS1 DS2 DS4 DS5 DS6 DS7 GIAB NA12878; do
  IN=$WD/$DS.fq; cp "${SRC[$DS]}" "$IN"; raw=$(stat -c %s "$IN")
  A=$WD/$DS.arcs; DEC=$WD/$DS.dec
  ct=$(btime $ARCS --chain-pg compress "$IN" "$A"); cr=$(peak $ARCS --chain-pg compress "$IN" "$A")
  arch=$(stat -c %s "$A")
  dt=$(btime $ARCS decompress "$A" "$DEC"); dr=$(peak $ARCS decompress "$A" "$DEC")
  ll=$(losscmp "$IN" "$DEC")
  printf "TOOL=ARCS DS=%s raw=%s archive=%s ctime=%s crss_kb=%s dtime=%s drss_kb=%s lossless=%s\n" \
    "$DS" "$raw" "$arch" "$ct" "$cr" "$dt" "$dr" "$ll" | tee $LOGD/${DS}__ARCS.log
  rm -f "$IN" "$A" "$DEC"
done
echo "ALLDONE"
