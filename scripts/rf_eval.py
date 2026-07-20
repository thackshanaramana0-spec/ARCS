import collections, math, re, statistics
OFF=1999000
conf=[]
for line in open('conf.bed'):
    f=line.split(); conf.append((int(f[1]),int(f[2])))
conf.sort()
def in_conf(p1):
    x=p1-1
    for s,e in conf:
        if s<=x<e: return True
        if s>x: break
    return False
truth_het=set(); truth_var=set(); B="ACGT"; B2I={'A':0,'C':1,'G':2,'T':3}
for line in open('truth.vcf'):
    if line[0]=='#': continue
    f=line.split('\t'); pos=int(f[1]); ref=f[3]; alt=f[4]; gt=f[9].split(':')[0]
    truth_var.add(pos)
    if len(ref)==1 and len(alt)==1 and alt in B and gt in ('0/1','0|1','1|0'): truth_het.add(pos)
# our pileup on our consensus (contig,pos)
col=collections.defaultdict(list)
for line in open('rf_reads.sam'):
    if line[0]=='@': continue
    f=line.split('\t'); flag=int(f[1]); mapq=int(f[4]); cig=f[5]; seq=f[9]; qual=f[10]; pos=int(f[3]); rn=f[2]
    if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20: continue
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
def loglik(rl,al):
    m=len(al);ll=0.0
    for b,q in rl:
        e=10**(-q/10.0); p=sum(((1-e) if b==a else e/3.0) for a in al); ll+=math.log(max(p/m,1e-300))
    return ll
def callvars(QUAL=20,dmin=10,altq=28):
    S=set()
    for k,rl in col.items():
        cnt=[0,0,0,0]
        for b,q in rl: cnt[b]+=1
        d=sum(cnt)
        if d<dmin: continue
        o=sorted(range(4),key=lambda b:-cnt[b]); M,mn=o[0],o[1]
        if cnt[mn]<2: continue
        if statistics.median([q for b,q in rl if b==mn])<altq: continue
        if 10.0*(loglik(rl,[M,mn])-loglik(rl,[M]))/math.log(10)>=QUAL: S.add(k)
    return S
# contig->genome coordinate map from contigs_to_ref.sam
cmap={}  # contig -> (ref_start0, cigar, rev)
for line in open('contigs_to_ref.sam'):
    if line[0]=='@': continue
    f=line.split('\t'); flag=int(f[1]); mapq=int(f[4])
    if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20: continue
    if f[0] in cmap: continue
    cmap[f[0]]=(int(f[3])-1, f[5], bool(flag&0x10), len(f[9]))
def contig_pos_to_genome(rn, cp):
    if rn not in cmap: return None
    rs, cig, rev, clen = cmap[rn]
    # if contig aligned reverse, position on contig maps to RC; convert
    target = (clen-1-cp) if rev else cp
    refidx=rs; readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            if readidx<=target<readidx+n: return OFF + refidx + (target-readidx)
            refidx+=n; readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
    return None
def near(p,S,t=3): return any((p+d) in S for d in range(-t,t+1))
S=callvars()
gcalls=set()
for (rn,cp) in S:
    g=contig_pos_to_genome(rn,cp)
    if g is not None: gcalls.add(g)
gcalls_conf=set(p for p in gcalls if in_conf(p))
het_conf=set(p for p in truth_het if in_conf(p))
tp=sum(1 for p in het_conf if near(p,gcalls_conf))
good=sum(1 for p in gcalls_conf if near(p,truth_var))
prec=good/max(len(gcalls_conf),1); rec=tp/max(len(het_conf),1)
print(f"REFERENCE-FREE (own consensus) on REAL data:")
print(f"  reads self-placed ~64%; het SNVs in confident={len(het_conf)}")
print(f"  our genome-mapped calls in confident={len(gcalls_conf)}")
print(f"  het-SNV site precision={prec:.3f}  recall={rec:.3f}  (TP={tp})")
