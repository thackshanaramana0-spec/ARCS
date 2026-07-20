#!/usr/bin/env python3
# Parse per-(dataset,tool) logs -> SUMMARY.md with 4 tables, every cell linking to its log.
import os, glob, re
LOGD="/mnt/c/Temp/bench8/logs"
OUT="/mnt/c/Temp/bench8/SUMMARY.md"
DATASETS=["DS1","DS2","DS4","DS5","DS6","DS7","GIAB","NA12878"]
TOOLS=["ARCS","SPRING","Genozip","fqzcomp","PgRC2"]
DESC={"DS1":"E.coli 151bp","DS2":"Human 127bp","DS4":"M.tb 51bp","DS5":"SARS 221bp",
      "DS6":"bacteria 37bp","DS7":"Human chr22 ~1M reads","GIAB":"HG002 150bp","NA12878":"HG001 148bp"}

def parse(ds,tool):
    p=f"{LOGD}/{ds}__{tool}.log"
    if not os.path.exists(p): return None
    d=dict(re.findall(r"(\w+)=(\S+)", open(p).read()))
    d["_log"]=f"logs/{ds}__{tool}.log"
    return d

data={(ds,t):parse(ds,t) for ds in DATASETS for t in TOOLS}

def link(val, logrel):
    return f"[{val}]({logrel})" if logrel else val

def cell(ds,tool,field,fmt=lambda x:x):
    r=data.get((ds,tool))
    if not r or field not in r: return "—"
    try: v=fmt(r[field])
    except: v=r[field]
    return link(v, r["_log"])

def mb(b): return f"{int(b)/1e6:.2f} MB"
def secs(s): return f"{float(s):.2f}s"
def ram(kb): return f"{int(kb)//1024} MB"

with open(OUT,"w") as f:
    f.write("# ARCS full benchmark — 8 datasets × 5 tools\n\n")
    f.write("Every value links to the raw log it was measured from. Native ext4, best-of-2 wall time, "
            "peak RAM via kernel VmHWM. Lossless verified by byte comparison.\n\n")
    f.write("Datasets: " + ", ".join(f"**{d}** ({DESC[d]})" for d in DATASETS) + "\n\n")

    f.write("## Table 1 — Archive size (compression ratio; smaller = better)\n\n")
    f.write("| Dataset | raw | " + " | ".join(TOOLS) + " |\n")
    f.write("|"+"---|"*(len(TOOLS)+2)+"\n")
    for ds in DATASETS:
        r=data.get((ds,"ARCS")); raw=mb(r["raw"]) if r and "raw" in r else "—"
        row=[ds, raw]+[cell(ds,t,"archive",mb) for t in TOOLS]
        f.write("| "+" | ".join(row)+" |\n")

    f.write("\n## Table 2a — Compression speed (wall time; lower = better)\n\n")
    f.write("| Dataset | " + " | ".join(TOOLS) + " |\n"); f.write("|"+"---|"*(len(TOOLS)+1)+"\n")
    for ds in DATASETS:
        f.write("| "+ds+" | "+" | ".join(cell(ds,t,"ctime",secs) for t in TOOLS)+" |\n")

    f.write("\n## Table 2b — Decompression speed (wall time; lower = better)\n\n")
    f.write("| Dataset | " + " | ".join(TOOLS) + " |\n"); f.write("|"+"---|"*(len(TOOLS)+1)+"\n")
    for ds in DATASETS:
        f.write("| "+ds+" | "+" | ".join(cell(ds,t,"dtime",secs) for t in TOOLS)+" |\n")

    f.write("\n## Table 3 — Peak RAM (compress / decompress, VmHWM)\n\n")
    f.write("| Dataset | " + " | ".join(TOOLS) + " |\n"); f.write("|"+"---|"*(len(TOOLS)+1)+"\n")
    for ds in DATASETS:
        cells=[]
        for t in TOOLS:
            r=data.get((ds,t))
            if r and "crss_kb" in r:
                cells.append(f"{link(ram(r['crss_kb']),r['_log'])} / {ram(r['drss_kb'])}")
            else: cells.append("—")
        f.write("| "+ds+" | "+" | ".join(cells)+" |\n")

    f.write("\n## Table 5 — Losslessness (byte-identical roundtrip)\n\n")
    f.write("| Dataset | " + " | ".join(TOOLS) + " |\n"); f.write("|"+"---|"*(len(TOOLS)+1)+"\n")
    for ds in DATASETS:
        f.write("| "+ds+" | "+" | ".join(cell(ds,t,"lossless") for t in TOOLS)+" |\n")

    f.write("\n## Table 4 — Reference-free variant calling (analysis-free byproduct)\n\n")
    f.write("Only datasets with a truth VCF (GIAB, NA12878, E.coli-sim). Values from the ref-free "
            "eval docs; see linked result files.\n\n")
    f.write("| Dataset | ARCS F1 | Kmer2SNP F1 | DiscoSNP++ | source |\n|---|---|---|---|---|\n")
    f.write("| NA12878 (held-out) | 0.816 | 0.833 | 0.71 | [REFFREE](../docs/REFFREE_NA12878_COMPARISON.md) |\n")
    f.write("| NA12878 (avg 3-region) | 0.840 | 0.818 | 0.71 | [REFFREE](../docs/REFFREE_NA12878_COMPARISON.md) |\n")
    f.write("| E.coli sim (haploid) | >0.99 P/R | — | — | [SHARED_LATENT](../docs/SHARED_LATENT_ERROR_STATE.md) |\n")

print("wrote", OUT)
