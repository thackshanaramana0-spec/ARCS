import lzma, re, math
names=[]
with open('giab150.fq') as f:
    while True:
        h=f.readline()
        if not h: break
        f.readline(); f.readline(); f.readline()
        names.append(h[1:].rstrip('\n'))
n=len(names)
raw=('\n'.join(names)+'\n').encode()
base=len(lzma.compress(raw, preset=9))
print(f"names={n}  LZMA-9={base} B = {8*base/n:.2f} bits/name")
# parse fields: instr:run:flowcell:lane:tile:x:y/mate
cols=[[] for _ in range(8)]
for nm in names:
    core,_,mate = nm.partition('/')
    f=core.split(':')
    for i in range(7): cols[i].append(f[i] if i<len(f) else '')
    cols[7].append(mate)
# per-field cardinality + information content
labels=['instr','run','flowcell','lane','tile','X','Y','mate']
total_bits=0
for i,lab in enumerate(labels):
    vals=cols[i]; card=len(set(vals))
    if all(v.isdigit() for v in vals if v):
        ints=[int(v) for v in vals if v]
        rng=max(ints)-min(ints)+1 if ints else 1
        bits=math.log2(rng) if rng>1 else 0
    else:
        bits=math.log2(card) if card>1 else 0
    total_bits+=bits
    print(f"  {lab:9s} cardinality={card:6d}  ~{bits:5.1f} bits/read")
floor_bytes = total_bits*n/8
print(f"INFORMATION FLOOR (independent fields): ~{floor_bytes:.0f} B = {total_bits:.1f} bits/name")
print(f"LZMA-9 is {100*(base/floor_bytes-1):+.1f}% vs this floor")
# what fraction is the irreducible X,Y coordinates?
xy_bits = 0
for i in (5,6):
    ints=[int(v) for v in cols[i] if v]; rng=max(ints)-min(ints)+1
    xy_bits+=math.log2(rng)
print(f"X+Y coords alone = {xy_bits:.1f} bits/read = {xy_bits*n/8:.0f} B ({100*xy_bits*n/8/base:.0f}% of the LZMA archive)")

print("\n=== VARIANT: binary-pack X,Y as raw ints + LZMA structural fields ===")
import struct
# structural = everything except X,Y, joined; XY = raw little-endian
struct_names=[]
xy=bytearray()
ok=True
for nm in names:
    core,_,mate=nm.partition('/')
    f=core.split(':')
    if len(f)!=7 or not f[5].isdigit() or not f[6].isdigit(): ok=False; break
    struct_names.append(f"{f[0]}:{f[1]}:{f[2]}:{f[3]}:{f[4]}/{mate}")
    xy+=struct.pack('<II', int(f[5]), int(f[6]))
if ok:
    sblob=('\n'.join(struct_names)+'\n').encode()
    s_lz=len(lzma.compress(sblob,preset=9))
    xy_lz=len(lzma.compress(bytes(xy),preset=9))
    xy_raw=len(xy)
    print(f"  structural LZMA-9 = {s_lz} B")
    print(f"  XY raw = {xy_raw} B, XY LZMA-9 = {xy_lz} B (use min)")
    tot = s_lz + min(xy_lz, xy_raw)
    print(f"  TOTAL = {tot} B  vs LZMA-9 baseline {base} B = {100*(1-tot/base):+.1f}%")
    # also try: XY delta-sorted? no, random. try 2-byte packing if ranges fit
