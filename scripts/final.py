import collections, math, statistics
truth={}
for line in open('ecoli_sim.truth'):
    a=line.split(); truth[(int(a[0]),int(a[1]))]=a[2]
cnt=collections.defaultdict(lambda:[0,0,0,0]); reads_at=collections.defaultdict(list); mm_rows=[]; pd={}
for line in open('sim_pileup.tsv'):
    if line[0]=='r': continue
    v=line.split()
    rid=int(v[0]);q=int(v[1]);dep=int(v[7]);mm=int(v[8]);pos=int(v[9]);col=int(v[10]);rb=int(v[12])
    if col<0 or rb>3: continue
    cnt[col][rb]+=1; reads_at[col].append((rb,q)); pd[(rid,pos)]=dep
    if mm: mm_rows.append((rid,pos,col,rb))
def loglik(rl,al):
    m=len(al);ll=0.0
    for b,q in rl:
        e=10**(-q/10.0); p=sum(((1-e) if b==a else e/3.0) for a in al); ll+=math.log(max(p/m,1e-300))
    return ll
def run(QUAL,dmin,minalt,altq):
    va={}
    for col,rl in reads_at.items():
        c=cnt[col]; d=sum(c)
        if d<dmin: continue
        o=sorted(range(4),key=lambda b:-c[b]); M,mn=o[0],o[1]
        if c[mn]<minalt: continue
        if altq and statistics.median([q for b,q in rl if b==mn])<altq: continue
        if 10.0*(loglik(rl,[M,mn])-loglik(rl,[M]))/math.log(10)>=QUAL: va[col]={M,mn}
    tpV=fpV=fnV=tpE=0;nV=nE=0; tp6=tot6=0
    for rid,pos,col,rb in mm_rows:
        t=truth.get((rid,pos)); cV=(col in va and rb in va[col])
        if t=='V':
            nV+=1; 
            if cV: tpV+=1
            else: fnV+=1
            if pd[(rid,pos)]>=6:
                tot6+=1
                if cV: tp6+=1
        elif t=='E':
            nE+=1
            if cV: fpV+=1
            else: tpE+=1
    pV=tpV/max(tpV+fpV,1); rV=tpV/max(nV,1); rV6=tp6/max(tot6,1)
    pE=tpE/max(tpE+fnV,1); rE=tpE/max(nE,1)
    return pV,rV,rV6,pE,rE
print("QUAL dmin mA altq | VarP  VarR  VarR(d>=6) | ErrP  ErrR")
for cfg in [(15,3,2,28),(12,3,2,20),(10,4,3,25),(15,5,3,25),(12,4,3,28)]:
    pV,rV,rV6,pE,rE=run(*cfg)
    f=" *ALL99@d>=6" if min(pV,rV6,pE,rE)>=0.99 else ""
    print(f" {cfg[0]:2d}  {cfg[1]}   {cfg[2]}  {cfg[3]:2d}  | {pV:.3f} {rV:.3f}  {rV6:.3f}     | {pE:.3f} {rE:.3f}{f}")
