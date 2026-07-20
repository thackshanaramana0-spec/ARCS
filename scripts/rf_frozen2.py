import collections,math,re,statistics,sys
# FROZEN bubble-caller: hd<=1, maf>=0.20, dhi=2.5, khi=1.1, mc>=3, tri=0.12.
# Bubble cleanliness = Hamming distance between major & minor haplotype flank consensus.
SAM,C2R,TRUTH,CONF,KC,KH=sys.argv[1:7]
HDMAX,MAF,DHI,KHI,MC,TRI=2,0.20,2.5,1.1,3,0.12
HALF=15;B2I={'A':0,'C':1,'G':2,'T':3};comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s):return ''.join(comp.get(c,'N') for c in reversed(s))
def canon(s):
    r=rc(s);return s if s<=r else r
kc={}
for line in open(KC):
    p=line.split()
    if len(p)==2: kc[p[0]]=int(p[1])
hist={}
for line in open(KH):
    a,c=line.split();hist[int(a)]=int(c)
H=max((a for a in hist if a>=5),key=lambda a:hist[a])
def kcount(km):
    if 'N' in km or len(km)!=31:return 0
    return kc.get(canon(km),0)
cf=[]
for line in open(CONF):
    f=line.split()
    if f[0]=='20': cf.append((int(f[1]),int(f[2])))
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
col=collections.defaultdict(list);recs=[]
for line in open(SAM):
    if line[0]=='@':continue
    f=line.split('\t');flag=int(f[1]);mapq=int(f[4]);cig=f[5];seq=f[9];qual=f[10];pos=int(f[3]);rn=f[2]
    if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20:continue
    recs.append((rn,pos-1,cig,seq));refidx=pos-1;readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            for t in range(n):
                b=seq[readidx+t]
                if b in B2I: col[(rn,refidx+t)].append((B2I[b],ord(qual[readidx+t])-33))
            refidx+=n;readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
def loglik(rl,al):
    m=len(al);ll=0.0
    for b,q in rl:
        e=10**(-q/10.0);p=sum(((1-e) if b==a else e/3.0) for a in al);ll+=math.log(max(p/m,1e-300))
    return ll
C={}
for k,rl in col.items():
    cnt=[0,0,0,0]
    for b,q in rl: cnt[b]+=1
    d=sum(cnt)
    if d<6:continue
    o=sorted(range(4),key=lambda b:-cnt[b]);M,mn=o[0],o[1]
    if cnt[mn]<2:continue
    if statistics.median([q for b,q in rl if b==mn])<20:continue
    if 10.0*(loglik(rl,[M,mn])-loglik(rl,[M]))/math.log(10)<10:continue
    C[k]=(M,mn,cnt,d)
med=sorted(v[3] for v in C.values())[len(C)//2] if C else 30
want=set(C.keys())
FLmin=collections.defaultdict(list);FLmaj=collections.defaultdict(list)
majoff=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0,0,0]))
minoff=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0,0,0]))
for rn,rs,cig,seq in recs:
    refidx=rs;readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            for t in range(n):
                key=(rn,refidx+t)
                if key in want:
                    ri=readidx+t;b=seq[ri]
                    if b in B2I:
                        isM=B2I[b]==C[key][0]; ismn=B2I[b]==C[key][1]
                        if (isM or ismn) and ri-HALF>=0 and ri+HALF+1<=len(seq):
                            (FLmaj if isM else FLmin)[key].append(seq[ri-HALF:ri+HALF+1])
                        if isM or ismn:
                            grp=majoff if isM else minoff
                            for off in range(-HALF,HALF+1):
                                if off==0: continue
                                p=ri+off
                                if 0<=p<len(seq) and seq[p] in B2I: grp[key][off][B2I[seq[p]]]+=1
            refidx+=n;readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
keep=set()
for key,(M,mn,cnt,d) in C.items():
    # bubble cleanliness (haplotype-flank Hamming)
    diffs=0
    for off in range(-HALF,HALF+1):
        if off==0: continue
        ma=majoff[key].get(off); mi=minoff[key].get(off)
        if not ma or not mi or sum(ma)<3 or sum(mi)<3: continue
        if max(range(4),key=lambda x:ma[x])!=max(range(4),key=lambda x:mi[x]): diffs+=1
    if diffs>HDMAX: continue
    if d>DHI*med: continue
    vmin=FLmin.get(key,[]);vmaj=FLmaj.get(key,[])
    cmaj=kcount(collections.Counter(vmaj).most_common(1)[0][0]) if vmaj else 0
    cmin=kcount(collections.Counter(vmin).most_common(1)[0][0]) if vmin else 0
    if max(cmaj,cmin)>KHI*H: continue
    if cnt[mn]/d<MAF: continue
    o=sorted(range(4),key=lambda b:-cnt[b])
    if cnt[o[2]]>TRI*d: continue
    mcnt=collections.Counter(vmin).most_common(1)[0][1] if vmin else 0
    if mcnt<MC: continue
    keep.add(key)
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
g=set()
for (rn,cp) in keep:
    gg=c2g(rn,cp)
    if gg is not None: g.add(gg)
gc=set(p for p in g if inc(p));hc=set(p for p in th if inc(p))
def near(p,S):return any((p+d) in S for d in range(-3,4))
tp=sum(1 for p in hc if near(p,gc));good=sum(1 for p in gc if near(p,tv))
P=good/max(len(gc),1);R=tp/max(len(hc),1);F=2*P*R/max(P+R,1e-9)
print(f"ARCS-BUBBLE(frozen) het={len(hc)} calls={len(gc)} TP={tp} | P={P:.3f} R={R:.3f} F1={F:.3f}")
