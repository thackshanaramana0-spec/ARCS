import collections
L=151
meta={}
for line in open('ecoli_sim.meta'):
    i,s,st=map(int,line.split()); meta[i]=(s,st)
truthE=set()
for line in open('ecoli_sim.truth'):
    a=line.split()
    if a[2]=='E': truthE.add((int(a[0]),int(a[1])))
varloci=set(int(line.split()[0]) for line in open('ecoli_sim.varloci'))
def gpos(rid,pos):
    s,st=meta[rid]; return s+(L-1-pos if st else pos)
col_mm=collections.defaultdict(list); col_depth={}
for line in open('sim_pileup.tsv'):
    if line[0]=='r': continue
    v=line.split()
    rid=int(v[0]);depth=int(v[7]);mm=int(v[8]);pos=int(v[9]);col=int(v[10]);rb=int(v[12])
    if col<0: continue
    col_depth[col]=depth
    if mm: col_mm[col].append((rid,pos,rb))
def near(g,S): return any((g+d) in S for d in (-1,0,1))
print("MAF thr | minCnt | Variant P | Variant R | Error P | Error R")
for maf,minc in [(0.05,3),(0.10,4),(0.15,4),(0.20,5),(0.30,6)]:
    call={}; var_cols=set()
    for col,m in col_mm.items():
        d=col_depth[col]
        c=collections.Counter(rb for _,_,rb in m)
        top,topn=c.most_common(1)[0]
        is_var = (d>=8 and topn>=minc and topn/d>=maf and topn/len(m)>=0.8)
        for rid,pos,rb in m: call[(rid,pos)]=('V' if (is_var and rb==top) else 'E')
        if is_var: var_cols.add(col)
    callE={k for k,v in call.items() if v=='E'}
    tp=len(callE & truthE); fp=len(callE-truthE)
    pe=tp/max(tp+fp,1); re=tp/max(len(truthE),1)
    cv=set()
    for col in var_cols:
        gs=collections.Counter(gpos(r,p) for r,p,_ in col_mm[col]); cv.add(gs.most_common(1)[0][0])
    tpv=sum(1 for g in cv if near(g,varloci))
    pv=tpv/max(len(cv),1); rv=sum(1 for g in varloci if near(g,cv))/len(varloci)
    print(f"  {maf:.2f}  |   {minc}    |   {pv:.3f}   |   {rv:.3f}   |  {pe:.3f} |  {re:.3f}")
