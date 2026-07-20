import random
random.seed(7)
seq=[]
with open('ecoli_ref.fa') as f:
    next(f)
    for line in f: seq.append(line.strip())
G=''.join(seq)
REGION=300000; L=151; COV=30
region=G[:REGION]; NREADS=REGION*COV//L
B="ACGT"; comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(comp[c] for c in reversed(s))
V=300; varloci={}
for p in random.sample(range(1,REGION-1), V):
    ref=region[p]
    if ref in B: varloci[p]=random.choice([b for b in B if b!=ref])
def perr_at(k): return 0.002 + 0.03*((k/L)**3)
fq=open('ecoli_sim.fastq','w'); truth=open('ecoli_sim.truth','w')
meta=open('ecoli_sim.meta','w'); vf=open('ecoli_sim.varloci','w')
nerr=0; nvar=0
for i in range(NREADS):
    start=random.randint(0,REGION-L); strand=random.random()<0.5
    bases=list(region[start:start+L]); qual=[0]*L; events=[]
    for k in range(L):
        gpos=start+k
        if gpos in varloci and random.random()<0.5:
            bases[k]=varloci[gpos]; events.append((k,'V')); nvar+=1
        if random.random()<perr_at(k):
            ob=bases[k]; bases[k]=random.choice([b for b in B if b!=ob])
            events.append((k,'E')); nerr+=1; q=max(2,int(random.gauss(10,4)))
        else: q=min(40,max(2,int(random.gauss(37,3))))
        qual[k]=q
    rs=''.join(bases); rq=''.join(chr(q+33) for q in qual)
    if strand:
        rs=rc(rs); rq=rq[::-1]; events=[(L-1-k,t) for (k,t) in events]
    fq.write(f"@sim.{i}\n{rs}\n+\n{rq}\n")
    meta.write(f"{i}\t{start}\t{int(strand)}\n")
    for (k,t) in events: truth.write(f"{i}\t{k}\t{t}\n")
for p,a in varloci.items(): vf.write(f"{p}\t{a}\n")
for h in (fq,truth,meta,vf): h.close()
print(f"reads={NREADS} varloci={len(varloci)} errors={nerr} variant-bases={nvar}")
