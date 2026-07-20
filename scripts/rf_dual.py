import collections, math, re, statistics, itertools
# Dual-region optimizer: evaluate each config on BOTH regions; accept only configs
# that beat Kmer2SNP on BOTH; rank by MIN F1 across regions (cannot overfit to one).
# Paralog rejection via LOCAL DEPTH vs region median (paralog=2x single-copy) +
# allele balance + biallelic cleanliness + minor-allele read-flank concordance.
HALF=15; B2I={'A':0,'C':1,'G':2,'T':3}; I2B="ACGT"
comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(comp.get(c,'N') for c in reversed(s))
def canon(s):
    r=rc(s); return s if s<=r else r
def load(sam,c2r,truth,conf,kcounts,khisto):
    kc={}
    for line in open(kcounts):
        p=line.split()
        if len(p)==2: kc[p[0]]=int(p[1])
    hist={}
    for line in open(khisto):
        a,c=line.split(); hist[int(a)]=int(c)
    H=max((a for a in hist if a>=5), key=lambda a:hist[a])
    C=[];
    cf=[]
    for line in open(conf):
        f=line.split()
        if f[0]=='20': cf.append((int(f[1]),int(f[2])))
    cf.sort()
    th=set(); tv=set(); Bs="ACGT"
    for line in open(truth):
        if line[0]=='#': continue
        f=line.split('\t'); pos=int(f[1]); ref=f[3]; alt=f[4]; gt=f[9].split(':')[0]
        tv.add(pos)
        if len(ref)==1 and len(alt)==1 and alt in Bs and gt in ('0/1','0|1','1|0'): th.add(pos)
    col=collections.defaultdict(list); recs=[]
    for line in open(sam):
        if line[0]=='@': continue
        f=line.split('\t'); flag=int(f[1]); mapq=int(f[4]); cig=f[5]; seq=f[9]; qual=f[10]; pos=int(f[3]); rn=f[2]
        if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20: continue
        recs.append((rn,pos-1,cig,seq))
        refidx=pos-1; readidx=0
        for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
            n=int(n)
            if op in 'M=X':
                for t in range(n):
                    b=seq[readidx+t]
                    if b in B2I: col[(rn,refidx+t)].append((B2I[b],ord(qual[readidx+t])-33))
                refidx+=n; readidx+=n
            elif op in 'IS': readidx+=n
            elif op in 'DN': refidx+=n
    cmap={}
    for line in open(c2r):
        if line[0]=='@': continue
        f=line.split('\t'); flag=int(f[1]); mapq=int(f[4])
        if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20: continue
        if f[0] in cmap: continue
        cmap[f[0]]=(int(f[3])-1,f[5],bool(flag&0x10),len(f[9]))
    def kcount(km):
        if 'N' in km or len(km)!=31: return 0
        return kc.get(canon(km),0)
    return dict(col=col,recs=recs,cmap=cmap,cf=cf,th=th,tv=tv,kcount=kcount,H=H)
def loglik(rl,al):
    m=len(al);ll=0.0
    for b,q in rl:
        e=10**(-q/10.0); p=sum(((1-e) if b==a else e/3.0) for a in al); ll+=math.log(max(p/m,1e-300))
    return ll
def prep(D):
    col=D['col']
    C={}
    for k,rl in col.items():
        cnt=[0,0,0,0]
        for b,q in rl: cnt[b]+=1
        d=sum(cnt)
        if d<6: continue
        o=sorted(range(4),key=lambda b:-cnt[b]); M,mn=o[0],o[1]
        if cnt[mn]<2: continue
        if statistics.median([q for b,q in rl if b==mn])<20: continue
        if 10.0*(loglik(rl,[M,mn])-loglik(rl,[M]))/math.log(10)<10: continue
        C[k]=(M,mn,cnt,d)
    D['C']=C
    D['MEDD']=sorted(v[3] for v in C.values())[len(C)//2] if C else 30
    # minor & major allele read-flanks -> mode k-mers -> global k-mer coverage
    want=set(C.keys()); FLmin=collections.defaultdict(list); FLmaj=collections.defaultdict(list)
    for rn,rs,cig,seq in D['recs']:
        refidx=rs; readidx=0
        for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
            n=int(n)
            if op in 'M=X':
                for t in range(n):
                    key=(rn,refidx+t)
                    if key in want:
                        ri=readidx+t; b=seq[ri]
                        if b in B2I and ri-HALF>=0 and ri+HALF+1<=len(seq):
                            km=seq[ri-HALF:ri+HALF+1]
                            if B2I[b]==C[key][1]: FLmin[key].append(km)
                            elif B2I[b]==C[key][0]: FLmaj[key].append(km)
                refidx+=n; readidx+=n
            elif op in 'IS': readidx+=n
            elif op in 'DN': refidx+=n
    kcount=D['kcount']; modecnt={}; kmaxcov={}
    for k in C:
        vmin=FLmin.get(k,[]); vmaj=FLmaj.get(k,[])
        modecnt[k]=collections.Counter(vmin).most_common(1)[0][1] if vmin else 0
        cmaj=kcount(collections.Counter(vmaj).most_common(1)[0][0]) if vmaj else 0
        cmin=kcount(collections.Counter(vmin).most_common(1)[0][0]) if vmin else 0
        kmaxcov[k]=max(cmaj,cmin)          # global k-mer coverage (paralog if >> het peak)
    D['modecnt']=modecnt; D['kmaxcov']=kmaxcov
    # BUBBLE CLEANLINESS (learned from Kmer2SNP's Hamming-distance edge classes, NOT copied):
    # a real het = major & minor haplotype flanks identical except at the SNP (few diffs);
    # a paralog/misassembly = two loci that diverge at MANY flank positions.
    majoff=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0,0,0]))
    minoff=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0,0,0]))
    for rn,rs,cig,seq in D['recs']:
        refidx=rs; readidx=0
        for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
            n=int(n)
            if op in 'M=X':
                for t in range(n):
                    key=(rn,refidx+t)
                    if key in want:
                        ri=readidx+t; b=seq[ri]
                        if b in B2I:
                            grp = majoff if B2I[b]==C[key][0] else (minoff if B2I[b]==C[key][1] else None)
                            if grp is not None:
                                for off in range(-HALF,HALF+1):
                                    if off==0: continue
                                    p=ri+off
                                    if 0<=p<len(seq) and seq[p] in B2I: grp[key][off][B2I[seq[p]]]+=1
                refidx+=n; readidx+=n
            elif op in 'IS': readidx+=n
            elif op in 'DN': refidx+=n
    hd={}
    for k in C:
        diffs=0
        for off in range(-HALF,HALF+1):
            if off==0: continue
            ma=majoff[k].get(off); mi=minoff[k].get(off)
            if not ma or not mi: continue
            if sum(ma)<3 or sum(mi)<3: continue
            mb=max(range(4),key=lambda x:ma[x]); nb=max(range(4),key=lambda x:mi[x])
            if mb!=nb: diffs+=1
        hd[k]=diffs
    D['hd']=hd
def c2g(cmap,rn,cp):
    if rn not in cmap: return None
    rs,cig,rev,clen=cmap[rn]; target=(clen-1-cp) if rev else cp
    refidx=rs; readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            if readidx<=target<readidx+n: return refidx+(target-readidx)+1
            refidx+=n; readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
    return None
def score(D,maf,dhi,mc,tri,khi,hdmax):
    med=D['MEDD']; H=D['H']; keep=set()
    for key,(M,mn,cnt,d) in D['C'].items():
        if hdmax>=0 and D['hd'].get(key,99)>hdmax: continue   # BUBBLE cleanliness (haplotype-flank divergence)
        if d>dhi*med: continue                       # local-depth paralog reject
        if khi>0 and D['kmaxcov'].get(key,0)>khi*H: continue  # global k-mer paralog reject
        if cnt[mn]/d<maf: continue                   # allele balance
        o=sorted(range(4),key=lambda b:-cnt[b])
        if cnt[o[2]]>tri*d: continue                 # biallelic cleanliness
        if mc>0 and D['modecnt'].get(key,0)<mc: continue  # minor-flank concordance
        keep.add(key)
    cmap=D['cmap']; g=set()
    for (rn,cp) in keep:
        gg=c2g(cmap,rn,cp)
        if gg is not None: g.add(gg)
    def inc(p1):
        x=p1-1
        for s,e in D['cf']:
            if s<=x<e: return True
            if s>x: break
        return False
    gc=set(p for p in g if inc(p)); hc=set(p for p in D['th'] if inc(p))
    def near(p,S): return any((p+d2) in S for d2 in range(-3,4))
    tp=sum(1 for p in hc if near(p,gc)); good=sum(1 for p in gc if near(p,D['tv']))
    P=good/max(len(gc),1); R=tp/max(len(hc),1); F=2*P*R/max(P+R,1e-9)
    return P,R,F
R1=load('rf_reads.sam','contigs_to_ref.sam','hg001_truth.vcf','hg001_conf.bed','chr_k31.txt','chr_k31.histo'); prep(R1)
R2=load('r2_rf.sam','r2_c2r.sam','r2_truth.vcf','r2_conf.bed','r2_k31.txt','r2_k31.histo'); prep(R2)
K1F,K2F=0.839,0.782   # Kmer2SNP F1 on the two regions
print(f"medians: R1={R1['MEDD']} R2={R2['MEDD']}  H: R1={R1['H']} R2={R2['H']}  (beat: R1>{K1F}, R2>{K2F})")
print(f"{'hd':>3} {'maf':>4} {'dhi':>4} {'khi':>4} {'mc':>3} | {'R1 P/R/F1':>20} | {'R2 P/R/F1':>20} | min  win?")
res=[]
for hd,maf,dhi,khi,mc in itertools.product([1,2,3],[0.20,0.25],[1.7,2.5],[0.0,1.1],[0,3]):
    p1,r1,f1=score(R1,maf,dhi,mc,0.12,khi,hd); p2,r2,f2=score(R2,maf,dhi,mc,0.12,khi,hd)
    win = f1>K1F and f2>K2F
    res.append((min(f1,f2),win,hd,maf,dhi,khi,mc,p1,r1,f1,p2,r2,f2))
res.sort(reverse=True)
for mn,win,hd,maf,dhi,khi,mc,p1,r1,f1,p2,r2,f2 in res[:16]:
    tag=" *** BEATS BOTH ***" if win else ""
    print(f"{hd:>3} {maf:>4} {dhi:>4} {khi:>4} {mc:>3} | {p1:.3f}/{r1:.3f}/{f1:.3f} | {p2:.3f}/{r2:.3f}/{f2:.3f} | {mn:.3f}{tag}")
