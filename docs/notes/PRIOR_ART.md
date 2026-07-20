# Prior-art check — cross-stream quality conditioning (idea B)

Searched 2026-07-16. Two ideas, separated:

- **(A) DNA compressor also codes quality** — reuse the CM engine as a *separate* quality model
  (context = position + previous quality). Independent streams. **Not novel** — fqzcomp / ACO
  already do context-mixing on quality; CM is the known-best.
- **(B) cross-stream coupling** — feed the SEQUENCE model's per-base state (prediction
  probability / surprise + mismatch-vs-pg flag) as CONTEXT into the quality coder. Couples the
  two into one coder. The real lever.

## Verdict: the CORE of B is prior art
- Match-flag / mismatch-type / surrounding-reference bases as quality context:
  **US patent 12125562** ("Quality value compression framework in aligned sequencing data
  based on novel contexts", MPEG-G lineage). Standard in CRAM/aligned pipelines.
- Sequence-predictability → shrink quality: **Quartz**, **GeneCodeq**, sequence-based smoothing
  (PMC6873394) — but LOSSY.
- Context-mixing on quality generally: fqzcomp, ACO (BMC 2022), ZPAQ.

## What is still defensible for us (combination novelty, do NOT overclaim)
Reference-**FREE** (self-assembled pseudogenome, no external genome)
+ model-internal **probability/surprise** (real-valued, not just a binary match flag)
+ fully **LOSSLESS**
+ **single unified integer codec** (one engine, seq state exported directly as quality context).
No searched tool does all four together. This is "new engineering combination" novelty —
same *kind* PgRC2's assembler has — not a new primitive.

## Sources
- US patent 12125562 — https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12125562
- ACO (BMC Bioinformatics 2022) — https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04712-z
- Quartz / sequence-based smoothing — https://pmc.ncbi.nlm.nih.gov/articles/PMC6873394/
- GeneCodeq — https://academic.oup.com/bioinformatics/article/32/20/3124/2196578
