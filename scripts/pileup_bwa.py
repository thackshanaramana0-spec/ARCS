import collections, math, re, statistics
L=151
B2I={'A':0,'C':1,'G':2,'T':3}
truth={}
for line in open('ecoli_sim.truth'):
    a=line.split(); truth[(int(a[0]),int(a[1]))]=a[2]
# parse SAM -> pileup per ref column
col=collections.defaultdict(list)   # refpos -> list of (base, q, rid, orig_readpos)
for line in open('aln.sam'):
    if line[0]=='@': continue
    f=line.split('\t')
    qname=f[0]; flag=int(f[1]); rname=f[2]; pos=int(f[3]); mapq=int(f[4]); cig=f[5]; seq=f[9]; qual=f[10]
    if flag & 0x4 or flag & 0x100 or flag & 0x800: continue
    if mapq<20: continue
    rid=int(qname.split('.')[1]); rev=bool(flag&0x10)
    refidx=pos-1; readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            for t in range(n):
                b=seq[readidx+t]
                if b in B2I:
                    q=ord(qual[readidx+t])-33
                    sp=readidx+t
                    orig= (L-1-sp) if rev else sp
                    col[(rname,refidx+t)].append((B2I[b],q,rid,orig))
            refidx+=n; readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
def loglik(rl,al):
    m=len(al);ll=0.0
    for b,q,_,_ in rl:
        e=10**(-q/10.0); p=sum(((1-e) if b==a else e/3.0) for a in al); ll+=math.log(max(p/m,1e-300))
    return ll
def call(QUAL,dmin,minalt,altq):
    va={}
    for c,rl in col.items():
        cnt=[0,0,0,0]
        for b,q,_,_ in rl: cnt[b]+=1
        d=sum(cnt)
        if d<dmin: continue
        o=sorted(range(4),key=lambda b:-cnt[b]); M,mn=o[0],o[1]
        if cnt[mn]<minalt: continue
        if altq and statistics.median([q for b,q,_,_ in rl if b==mn])<altq: continue
        if 10.0*(loglik(rl,[M,mn])-loglik(rl,[M]))/math.log(10)>=QUAL: va[c]=(M,{M,mn})
    return va
def evalc(va):
    tpV=fpV=fnV=tpE=0;nV=nE=0; tp6=tot6=0
    for c,rl in col.items():
        cnt=[0,0,0,0]
        for b,q,_,_ in rl: cnt[b]+=1
        M=max(range(4),key=lambda b:cnt[b])
        alleles=va[c][1] if c in va else None
        for b,q,rid,orig in rl:
            if b==M: continue          # matches consensus, not a mismatch
            t=truth.get((rid,orig)); cV=(alleles is not None and b in alleles)
            if t=='V':
                nV+=1
                if cV: tpV+=1
                else: fnV+=1
                if len(rl)>=6:
                    tot6+=1
                    if cV: tp6+=1
            elif t=='E':
                nE+=1
                if cV: fpV+=1
                else: tpE+=1
    pV=tpV/max(tpV+fpV,1); rV=tpV/max(nV,1); rV6=tp6/max(tot6,1)
    pE=tpE/max(tpE+fnV,1); rE=tpE/max(nE,1)
    return pV,rV,rV6,pE,rE
print("aligned pileup columns:",len(col))
print("QUAL dmin mA altq | VarP  VarR  VarR(d>=6) | ErrP  ErrR")
for cfg in [(15,3,2,28),(20,4,2,28),(30,6,2,0),(15,4,2,25)]:
    pV,rV,rV6,pE,rE=evalc(call(*cfg))
    flag="  <-- all99" if min(pV,rV,pE,rE)>=0.99 else ("  99@d6" if min(pV,rV6,pE,rE)>=0.99 else "")
    print(f" {cfg[0]:2d}  {cfg[1]}   {cfg[2]}  {cfg[3]:2d}  | {pV:.3f} {rV:.3f}  {rV6:.3f}     | {pE:.3f} {rE:.3f}{flag}")
