# Deep literature survey — "Has anyone in the world done our style?"

Survey date: 2026-07-16. Written to be read cold by a fresh session with no prior context.

## 0. What "our style" is (so we can check it honestly)
ARCS `--chain-pg` = a 4-part pipeline:
1. **Self-assembled pseudogenome** built from the reads (no external reference).
2. Reads **mapped onto it** and stored as (position, orientation, mismatch list, N list) — the read
   itself is not stored, only its differences vs the assembled sequence.
3. A **dependency-free context-mixing DNA entropy coder** (ARCS-DNA) compresses the assembled sequence.
4. (Planned, idea B) **cross-stream quality**: use the sequence coder's per-base surprise + mismatch
   flag as context for the quality coder, in one unified lossless codec.

The blunt question: is any of this world-first? Answer below, component by component. **Short version:
almost every individual component has strong prior art. The specific *combination* is not a single
existing tool, but we must NOT claim the pieces are new.**

---

## 1. Component-by-component prior art

### 1a. Assemble reads into a reference, store reads as edits vs it
**Prior art — well established.** This is the dominant high-ratio strategy.
- **PgRC** (Kowalski & Grabowski, Bioinformatics 2020) — the "pseudogenome": approximate shortest
  common superstring over high-quality reads, then map reads onto it. This is our direct ancestor
  and the thing we benchmark against. https://academic.oup.com/bioinformatics/article/36/7/2082/5670526
- **CURC** (2022) — CUDA pseudogenome compressor, same family. https://academic.oup.com/bioinformatics/article/38/12/3294/6586792
- **Assembltrie / light assembly** (Ginart et al., Nature Communications 2018) — proves that
  assembly-based representations give the shortest output; stores reads in a compact trie.
  https://www.nature.com/articles/s41467-017-02480-6

### 1b. Multi-contig assembly + consensus + encode reads vs consensus  ← our assembler
**Prior art — THIS IS THE CLOSEST MATCH, and it is important to be honest about it.**
- **NanoSpring** (Meng, Chandak et al., Scientific Reports 2023 / bioRxiv 2021) does *exactly* our
  assembler shape: "reads are first **assembled into contigs**, a **consensus sequence** is obtained
  for each contig, and the consensus sequence and **encoded reads with respect to the consensus** are
  stored." That is multi-contig + consensus polish + map-reads-with-differences — the same three moves
  as `build_multicontig_pg`. https://www.nature.com/articles/s41598-023-29267-8
  - **The one real distinction:** NanoSpring targets **nanopore long reads** (thousands of bp, ~10%+
    error), using minimizer-based approximate assembly + MSA consensus. Ours targets **short Illumina
    reads** (~150 bp, <1% error) competing with PgRC, using a global k-mer index + majority-vote
    consensus + emit-in-position-order. Same *idea*, different regime and implementation.
- So: **the "multi-contig + consensus + encode-vs-consensus" idea is NOT world-first.** Our assembler
  is a legitimate independent re-implementation for the short-read regime, not a new concept.

### 1c. Reorder reads to improve compression (order-not-preserving)
**Prior art — standard.** ORCOM, Mince, SCALCE, BEETL, SPRING, Minicom all reorder reads first.
Our multi-contig emits reads in pg-position order = the same reorder-allowed contract. Not novel.
- SPRING: https://github.com/shubhamchandak94/SPRING · overview of the family:
  https://academic.oup.com/bioinformatics/article/34/4/558/4386919 (hash-based reordering).

### 1d. Context-mixing entropy coder for DNA  ← ARCS-DNA
**Prior art — the whole field.** Multi-order finite-context models + inverted-repeat handling +
mixing + arithmetic coding is exactly GeCo/GeCo2/**GeCo3**, XM, cmix, PAQ, Jarvis, NAF.
- **GeCo3** (Silva, Pratas, Pinho, GigaScience 2020) — neural-net mixing of context + substitution-
  tolerant models, explicit inverted-repeat sub-programs. https://academic.oup.com/gigascience/article/9/11/giaa119/5974977
- Our ARCS-DNA is a from-scratch, dependency-free re-implementation of this well-known architecture
  (FCM + IR + logistic mix + APM + arithmetic). **The architecture is not novel**; the value is that
  it is self-contained, 3.4× faster / 3.5× lighter than GeCo3 (at ~18% worse ratio) — an
  engineering trade-off, not a new algorithm.

### 1e. Quality conditioned on sequence / mismatch context  ← idea B (planned)
**Core is prior art** (see `PRIOR_ART.md` for detail):
- Match-flag / mismatch-type / surrounding-reference bases as quality context: **US patent 12125562**
  (MPEG-G lineage); standard in CRAM. https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12125562
- Sequence-predictability → shrink quality: Quartz, GeneCodeq (both LOSSY).
- Context-mixing on quality: fqzcomp, ACO.

---

## 2. Closest existing tools — one-line each

| Tool | Overlap with us | Key difference from us |
|---|---|---|
| **PgRC** | pseudogenome + map reads (1a,1c) | no consensus multi-contig; PPMd streams; our benchmark target |
| **NanoSpring** | **contig + consensus + encode-vs-consensus (1b)** | **nanopore long reads, not short Illumina** — closest match |
| **SPRING** | reorder + assemble-ish + all streams | different reorder + BSC/entropy; no pseudogenome |
| **Assembltrie** | assembly-optimal representation (1a) | trie representation, not map+mismatch |
| **GeCo3** | CM DNA coder (1d) | our ARCS-DNA reimplements this dependency-free |
| **CRAM / MPEG-G** | mismatch-as-quality-context (1e) | reference-based/aligned; ours reference-free + surprise |

---

## 3. Honest verdict

**Has anyone done our exact end-to-end pipeline as one tool? No single tool combines
(self-assembled short-read pseudogenome via multi-contig+consensus) + (dependency-free CM DNA coder)
+ (cross-stream lossless quality).** But every *component* has clear prior art:
- assemble+map: PgRC, Assembltrie
- multi-contig + consensus + encode-vs-consensus: **NanoSpring (near-identical idea, different regime)**
- reorder: ORCOM/Mince/SPRING
- CM DNA coder: GeCo3/XM/cmix
- quality-from-mismatch: MPEG-G patent / CRAM

**Therefore the defensible novelty is a COMBINATION / integration claim, not a new-primitive claim:**
> A reference-free short-read compressor that unifies self-assembled multi-contig consensus with a
> self-contained context-mixing coder and (idea B) couples sequence-model surprise into a lossless
> quality coder in a single codec.

That is the *same kind* of novelty PgRC has (its contribution was the pseudogenome assembler +
integration, not new entropy math). It is real but modest. **Do NOT claim the multi-contig assembler,
the CM coder, or mismatch-as-quality-context are individually new — reviewers will cite NanoSpring,
GeCo3, and MPEG-G immediately.** The strongest genuinely-fresh angle left is idea B's specific form:
reference-free + real-valued model *surprise* (not just a binary flag) + fully lossless + unified.

## 4. Sources
- PgRC — https://academic.oup.com/bioinformatics/article/36/7/2082/5670526
- NanoSpring — https://www.nature.com/articles/s41598-023-29267-8
- Assembltrie (light assembly) — https://www.nature.com/articles/s41467-017-02480-6
- SPRING — https://github.com/shubhamchandak94/SPRING · reorder analysis https://academic.oup.com/bioinformatics/article/34/4/558/4386919
- CURC — https://academic.oup.com/bioinformatics/article/38/12/3294/6586792
- GeCo3 — https://academic.oup.com/gigascience/article/9/11/giaa119/5974977
- MPEG-G quality-context patent — https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12125562
- Quartz / smoothing — https://pmc.ncbi.nlm.nih.gov/articles/PMC6873394/
- GeneCodeq — https://academic.oup.com/bioinformatics/article/32/20/3124/2196578
- ACO quality CM — https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04712-z
