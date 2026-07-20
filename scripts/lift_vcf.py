import re,sys
# Lift an ARCS contig-coordinate call VCF to genome coordinates (allele-aware) so a
# standard tool (rtg vcfeval / hap.py) can score it. Uses the contig->reference
# alignment (C2R, eval-only) to place each call, fetches the genome REF base from the
# reference FASTA, and complements alleles for reverse-strand contigs. Emits a valid
# genome-coordinate VCF (GT=0/1, or 1/2 when neither observed allele is the ref).
#   usage: lift_vcf.py <calls.contig.vcf> <c2r.sam> <ref.fa> <chrom> <out.vcf>
CALLS,C2R,REF,CHROM,OUT=sys.argv[1:6]
comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}

# reference sequence for CHROM (1-based indexing via seq[pos-1])
seq=[]
cur=None
for line in open(REF):
    if line[0]=='>':
        cur=line[1:].split()[0]
        continue
    if cur==CHROM: seq.append(line.strip())
refseq=''.join(seq)
def gref(pos1): return refseq[pos1-1].upper() if 1<=pos1<=len(refseq) else 'N'

# contig -> (genome_start0, cigar, reverse, contig_len)
cmap={}
for line in open(C2R):
    if line[0]=='@':continue
    f=line.split('\t');flag=int(f[1]);mapq=int(f[4])
    if flag&0x4 or flag&0x100 or flag&0x800 or mapq<20:continue
    if f[0] in cmap:continue
    cmap[f[0]]=(int(f[3])-1,f[5],bool(flag&0x10),len(f[9]))
def c2g(rn,cp0):
    if rn not in cmap:return None
    rs,cig,rev,clen=cmap[rn];target=(clen-1-cp0) if rev else cp0
    refidx=rs;readidx=0
    for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
        n=int(n)
        if op in 'M=X':
            if readidx<=target<readidx+n: return (refidx+(target-readidx)+1, rev)  # 1-based genome pos
            refidx+=n;readidx+=n
        elif op in 'IS': readidx+=n
        elif op in 'DN': refidx+=n
    return None

rows=[]
for line in open(CALLS):
    if line[0]=='#':continue
    f=line.rstrip('\n').split('\t')
    rn=f[0];cp1=int(f[1]);cref=f[3].upper();calt=f[4].upper()
    if len(cref)!=1 or len(calt)!=1: continue          # SNVs only
    r=c2g(rn,cp1-1)
    if r is None: continue
    gpos,rev=r
    gR=gref(gpos)
    if gR=='N': continue
    # observed alleles in genome frame (major=cref, minor=calt), complement if reverse
    a1 = comp.get(cref,'N') if rev else cref
    a2 = comp.get(calt,'N') if rev else calt
    if a1=='N' or a2=='N': continue
    # standard VCF representation: REF=gR, ALT=non-ref observed allele(s)
    if   a1==gR and a2!=gR: alt=a2; gt="0/1"
    elif a2==gR and a1!=gR: alt=a1; gt="0/1"
    elif a1!=gR and a2!=gR:
        # neither observed allele is the reference -> both are ALT
        if a1==a2: alt=a1; gt="1/1"
        else:      alt=a1+","+a2; gt="1/2"
    else:
        continue                                        # both == ref, not a variant
    rows.append((gpos,gR,alt,gt))

rows.sort()
with open(OUT,'w') as o:
    o.write("##fileformat=VCFv4.2\n")
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write("##contig=<ID=%s>\n"%CHROM)
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    seen=set()
    for gpos,gR,alt,gt in rows:
        if gpos in seen: continue                       # dedup positions (one call per site)
        seen.add(gpos)
        o.write("%s\t%d\t.\t%s\t%s\t30\tPASS\t.\tGT\t%s\n"%(CHROM,gpos,gR,alt,gt))
print("lifted %d calls -> %s"%(len(seen),OUT))
