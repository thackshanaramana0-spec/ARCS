import collections, math, statistics
truth={}
for line in open('ecoli_sim.truth'):
    a=line.split(); truth[(int(a[0]),int(a[1]))]=a[2]
cnt=collections.defaultdict(lambda:[0,0,0,0]); reads_at=collections.defaultdict(list); mm_rows=[]
for line in open('sim_pileup.tsv'):
    if line[0]=='r': continue
    v=line.split()
    rid=int(v[0]);q=int(v[1]);mm=int(v[8]);pos=int(v[9]);col=int(v[10]);rb=int(v[12])
    if col<0 or rb>3: continue
    cnt[col][rb]+=1; reads_at[col].append((rb,q))
    if mm: mm_rows.append((rid,pos,col,rb))
def loglik(rl,al):
    m=len(al);ll=0.0
    for b,q in rl:
        e=10**(-q/10.0); p=sum(((1-e) if b==a else e/3.0) for a in al); ll+=math.log(max(p/m,1e-300))
    return ll
def call(QUAL,dmin,minalt,altq):
    va={}
    for col,rl in reads_at.items():
        c=cnt[col]; d=sum(c)
        if d<dmin: continue
        o=sorted(range(4),key=lambda b:-c[b]); M,mn=o[0],o[1]
        if c[mn]<minalt: continue
        aq=[q for b,q in rl if b==mn]
        if statistics.median(aq)<altq: continue
        qual=10.0*(loglik(rl,[M,mn])-loglik(rl,[M]))/math.log(10)
        if qual>=QUAL: va[col]={M,mn}
    return va
def evalv(va):
    tpV=fpV=fnV=tpE=0; nV=nE=0
    for rid,pos,col,rb in mm_rows:
        t=truth.get((rid,pos)); callV=(col in va and rb in va[col])
        if t=='V':
            nV+=1
            if callV: tpV+=1
            else: fnV+=1
        elif t=='E':
            nE+=1
            if callV: fpV+=1
            else: tpE+=1
    pV=tpV/max(tpV+fpV,1); rV=tpV/max(nV,1)
    pE=tpE/max(tpE+fnV,1); rE=tpE/max(nE,1)
    return pV,rV,pE,rE
print("QUAL dmin minAlt altqMin | VarP  VarR | ErrP  ErrR")
for QUAL,dmin,ma,aq in [(30,6,2,0),(20,4,2,20),(15,4,2,25),(12,3,2,28),(15,3,2,30),(20,3,2,30)]:
    pV,rV,pE,rE=evalv(call(QUAL,dmin,ma,aq))
    flag=" <-- all99+" if min(pV,rV,pE,rE)>=0.99 else ""
    print(f" {QUAL:3d}  {dmin}    {ma}     {aq:3d}   | {pV:.3f} {rV:.3f} | {pE:.3f} {rE:.3f}{flag}")
