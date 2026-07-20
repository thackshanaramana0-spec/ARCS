import re,sys
# Score an ARCS-caller VCF (contig coordinates) against a reference-coordinate truth
# set. The contig->genome mapping (C2R = contigs aligned to the reference) is used
# ONLY for evaluation — the caller itself is reference-free. Mirrors the c2g/inc/near
# logic of rf_frozen2.py so numbers are comparable.
VCF,C2R,TRUTH,CONF=sys.argv[1:5]
cmap={}
for line in open(C2R):
    if line[0]=='@':continue
    f=line.split('\t');flag=int(f[1]);mapq=int(f[4])
    if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20:continue
    if f[0] in cmap:continue
    cmap[f[0]]=(int(f[3])-1,f[5],bool(flag&0x10),len(f[9]))
def c2g(rn,cp):
    if rn not in cmap:return None
    rs,cig,rev,clen=cmap[rn];target=(clen-1-cp) if rev else cp
    refidx=rs;readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            if readidx<=target<readidx+n: return refidx+(target-readidx)+1
            refidx+=n;readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
    return None
cf=[]
for line in open(CONF):
    f=line.split()
    if f and f[0]=='20': cf.append((int(f[1]),int(f[2])))
cf.sort()
def inc(p1):
    x=p1-1
    for s,e in cf:
        if s<=x<e:return True
        if s>x:break
    return False
th=set();tv=set();Bs="ACGT"
for line in open(TRUTH):
    if line[0]=='#':continue
    f=line.split('\t');pos=int(f[1]);ref=f[3];alt=f[4];gt=f[9].split(':')[0]
    tv.add(pos)
    if len(ref)==1 and len(alt)==1 and alt in Bs and gt in ('0/1','0|1','1|0'): th.add(pos)
g=set()
for line in open(VCF):
    if line[0]=='#':continue
    f=line.split('\t');rn=f[0];cp=int(f[1])-1
    gg=c2g(rn,cp)
    if gg is not None: g.add(gg)
gc=set(p for p in g if inc(p));hc=set(p for p in th if inc(p))
def near(p,S):return any((p+d) in S for d in range(-3,4))
tp=sum(1 for p in hc if near(p,gc));good=sum(1 for p in gc if near(p,tv))
P=good/max(len(gc),1);R=tp/max(len(hc),1);F=2*P*R/max(P+R,1e-9)
print(f"ARCS-call(cpp) het={len(hc)} calls={len(gc)} TP={tp} | P={P:.3f} R={R:.3f} F1={F:.3f}")
