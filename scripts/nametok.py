import lzma, re
names=[]
with open('giab150.fq') as f:
    while True:
        h=f.readline()
        if not h: break
        f.readline(); f.readline(); f.readline()
        names.append(h[1:].rstrip('\n'))
raw=('\n'.join(names)+'\n').encode()
print(f"names={len(names)}  raw={len(raw)}")
def L(b,preset=9): return len(lzma.compress(b, preset=preset))
# ---- baseline: what we do now (LZMA newline-joined) ----
base = L(raw, 6)
base9 = L(raw, 9)
print(f"[baseline] LZMA-6 = {base}   LZMA-9 = {base9}")
# ---- tokenized column-major + numeric delta ----
# tokenize each name into maximal runs of digits / non-digits
tok_re = re.compile(r'(\d+|\D+)')
toks = [tok_re.findall(n) for n in names]
ncol = max(len(t) for t in toks)
# transpose into columns; for numeric columns delta-encode
cols_out = []
for c in range(ncol):
    colvals = [ (t[c] if c < len(t) else '') for t in toks ]
    # numeric if ALL nonempty entries are pure digits
    nonempty = [v for v in colvals if v!='']
    is_num = nonempty and all(v.isdigit() for v in nonempty)
    if is_num:
        # zigzag delta of integer value, keep leading-zero width via string fallback if needed
        prev=0; buf=bytearray()
        # detect fixed leading zeros: if any value has leading zero and len>1, treat as string
        haslz = any(len(v)>1 and v[0]=='0' for v in nonempty)
        if haslz:
            cols_out.append(('\n'.join(colvals)).encode())
        else:
            for v in colvals:
                iv = int(v) if v!='' else prev
                d = iv-prev; prev=iv
                zz=(d<<1)^(d>>63)
                while True:
                    b=zz&0x7f; zz>>=7
                    if zz: buf.append(b|0x80)
                    else: buf.append(b); break
            cols_out.append(bytes(buf))
    else:
        cols_out.append(('\n'.join(colvals)).encode())
blob = b'\x00'.join(cols_out)
tok = L(blob, 9)
print(f"[tokenized] column-major + delta, LZMA-9 = {tok}")
print(f"  vs LZMA-9 baseline: {100*(1-tok/base9):.1f}% smaller")
print(f"  vs LZMA-6 (current default): {100*(1-tok/base):.1f}% smaller")
