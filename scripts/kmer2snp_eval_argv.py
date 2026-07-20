import re,sys,collections
REGFA,TRUTH,CONF,SNP,NONSEP = sys.argv[1:6]
comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(comp.get(c,'N') for c in reversed(s))
reg=[]; hdr=None
for line in open(REGFA):
    if line[0]=='>': hdr=line.strip()
    else: reg.append(line.strip())
REF=''.join(reg); RSTART=int(re.search(r':(\d+)-',hdr).group(1))
K=31; MID=15; sig={}
for i in range(len(REF)-K+1):
    kk=REF[i:i+K]
    if 'N' in kk: continue
    sig.setdefault(kk[:MID]+kk[MID+1:], RSTART+i+MID)
for i in range(len(REF)-K+1):
    kk=rc(REF[i:i+K])
    if 'N' in kk: continue
    sig.setdefault(kk[:MID]+kk[MID+1:], RSTART+i+MID)
conf=[]
for line in open(CONF):
    f=line.split()
    if f[0]=='20': conf.append((int(f[1]),int(f[2])))
conf.sort()
def in_conf(p1):
    x=p1-1
    for s,e in conf:
        if s<=x<e: return True
        if s>x: break
    return False
truth_het=set(); truth_var=set(); B="ACGT"
for line in open(TRUTH):
    if line[0]=='#': continue
    f=line.split('\t'); pos=int(f[1]); ref=f[3]; alt=f[4]; gt=f[9].split(':')[0]
    truth_var.add(pos)
    if len(ref)==1 and len(alt)==1 and alt in B and gt in ('0/1','0|1','1|0'): truth_het.add(pos)
calls=set()
for src in (SNP,NONSEP):
    try: lines=open(src)
    except: continue
    for line in lines:
        p=line.split()
        if len(p)<2: continue
        for kmer in (p[0],p[1]):
            if len(kmer)!=K: continue
            key=kmer[:MID]+kmer[MID+1:]
            if key in sig: calls.add(sig[key]); break
cc=set(p for p in calls if in_conf(p)); hc=set(p for p in truth_het if in_conf(p))
def near(p,S): return any((p+d) in S for d in range(-3,4))
tp=sum(1 for p in hc if near(p,cc)); good=sum(1 for p in cc if near(p,truth_var))
P=good/max(len(cc),1); R=tp/max(len(hc),1); F=2*P*R/max(P+R,1e-9)
print(f"KMER2SNP het conf={len(hc)} calls={len(cc)} TP={tp} | P={P:.3f} R={R:.3f} F1={F:.3f}")
