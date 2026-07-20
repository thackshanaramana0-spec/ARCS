import math
from collections import defaultdict
names=[]; quals=[]
with open('giab150.fq') as f:
    while True:
        h=f.readline()
        if not h: break
        s=f.readline(); p=f.readline(); q=f.readline().rstrip('\n')
        names.append(h[1:].strip()); quals.append(q)
# parse tile:x:y  (fields: instr run flowcell lane tile x y  [/mate])
tiles=[];xs=[];ys=[]
for nm in names:
    core=nm.split()[0].split('/')[0]
    f=core.split(':')
    # last three colon fields = tile,x,y
    tiles.append(int(f[-3])); xs.append(int(f[-2])); ys.append(int(f[-1]))
xmn,xmx=min(xs),max(xs); ymn,ymx=min(ys),max(ys)
tset=sorted(set(tiles))
print(f"reads={len(names)}  X:[{xmn},{xmx}] Y:[{ymn},{ymx}]  tiles={len(tset)}")
NB=8
def xbin(x): return min(NB-1,(x-xmn)*NB//max(1,(xmx-xmn)))
def ybin(y): return min(NB-1,(y-ymn)*NB//max(1,(ymx-ymn)))

# train on even-indexed reads, test on odd (honest held-out cross-entropy)
def build_eval(ctxfn):
    cnt=defaultdict(lambda: defaultdict(int))
    for i in range(0,len(quals),2):       # train
        L=quals[i]; q1=q2=run=0
        for j,ch in enumerate(L):
            qv=ord(ch)-33; k=ctxfn(i,j,q1,q2,run,len(L))
            cnt[k][qv]+=1
            run=min(run+1,7) if qv==q1 else 0
            q2=q1; q1=qv
    tot={k:sum(d.values()) for k,d in cnt.items()}
    bits=0.0; nsym=0
    for i in range(1,len(quals),2):        # test
        L=quals[i]; q1=q2=run=0
        for j,ch in enumerate(L):
            qv=ord(ch)-33; k=ctxfn(i,j,q1,q2,run,len(L))
            d=cnt.get(k); s=tot.get(k,0)
            # Laplace-smoothed predictive prob over 48 symbols (unseen ctx -> uniform)
            p=((d.get(qv,0) if d else 0)+1)/((s if s else 0)+48)
            bits-=math.log2(p); nsym+=1
            run=min(run+1,7) if qv==q1 else 0
            q2=q1; q1=qv
    return bits/nsym
base=lambda i,j,q1,q2,run,L:(q1,q2,j>>3,run)
f2d =lambda i,j,q1,q2,run,L:(q1,q2,j>>3,run,xbin(xs[i]),ybin(ys[i]))
ftile=lambda i,j,q1,q2,run,L:(q1,q2,j>>3,run,xbin(xs[i]),ybin(ys[i]),tiles[i])
print(f"baseline  H(q|q1,q2,pos,run)            = {build_eval(base):.4f} bpq")
print(f"+2D       H(...,xbin8,ybin8)            = {build_eval(f2d):.4f} bpq")
print(f"+2D+tile  H(...,xbin8,ybin8,tile)       = {build_eval(ftile):.4f} bpq")

# coarser spatial bins — does ANY spatial signal survive the sparsity cost?
for nbb in (4,3,2):
    def xb(x,n=nbb): return min(n-1,(x-xmn)*n//max(1,(xmx-xmn)))
    def yb(y,n=nbb): return min(n-1,(y-ymn)*n//max(1,(ymx-ymn)))
    fc=lambda i,j,q1,q2,run,L,n=nbb:(q1,q2,j>>3,run,xb(xs[i]),yb(ys[i]))
    print(f"+2D {nbb}x{nbb}   = {build_eval(fc):.4f} bpq  (baseline 1.5941)")
