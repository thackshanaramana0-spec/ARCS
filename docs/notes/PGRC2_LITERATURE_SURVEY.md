# Deep literature survey — PgRC / PgRC2's own novelty vs prior art

Survey date: 2026-07-16. Companion to `LITERATURE_SURVEY.md` (which does the same for OUR pipeline).
Purpose: apply the exact same honest, component-by-component prior-art test to the tool we benchmark
against, so we can see what PgRC2 actually invented vs borrowed — and whether our novelty is the
"same kind."

Note on naming: **PgRC** = Kowalski & Grabowski, Bioinformatics 2020. **PgRC2** = the v2/benchmark
binary used in this project (`C:\Temp\arcs_test\PgRC2.exe`); the GeCo3-upgrade experiments in this
repo (`PGRC2_GeCo3_UPGRADE_PAPER.md`) are OURS, not the original authors'. Component knowledge below
comes from the PgRC papers + our own line-by-line clone (Trial18, see memory).

## 0. PgRC's pipeline, component by component
1. Classify reads **HQ / LQ / N** using quality (an error-level threshold).
2. Build a **pseudogenome** = approximate **shortest common superstring (SCS)** over the HQ reads
   (greedy overlap assembly).
3. **Reverse-complement folding** — reads matched on either strand.
4. **Approximately map** the LQ reads onto the pseudogenome (store position + orientation + mismatches).
5. Serialize semantic streams (offsets, RC flags, mismatch counts/positions/bases) and entropy-code
   them (**PPMd / range coder**; PgRC2-in-this-project swaps in **GeCo3** for some streams).

---

## 1. Component-by-component prior art

### 1a. Build a reference FROM the reads, then encode reads against it  ← PgRC's headline
**CORRECTION (2026-07-16): PgRC was NOT first to "construct a reference from the reads."** That is a
well-populated subfield that predates PgRC (2020) by years:
- **Leon** (Benoit et al., BMC Bioinformatics 2015) — builds a de-novo **de Bruijn graph reference
  from the reads** (in a Bloom filter) and encodes each read as a path. Reference-built-from-reads,
  5 years before PgRC. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0709-7
- **BdBG** (Wang et al., 2018) — bucketed **dynamic de Bruijn graph**, reference-free, reads encoded
  as bifurcation lists.
- **ORCOM / Mince** (2015) — cluster overlapping reads (a local implicit reference).
- **Assembltrie / light assembly** (Ginart et al., Nature Comms 2018) — assembly-based representation
  proven near-optimal; reads in a compact trie. https://www.nature.com/articles/s41467-017-02480-6

So the *concept* — a reference synthesized from the reads + encode-reads-against-it — is prior art,
not a PgRC invention. **PgRC's actual, much narrower contribution** is the specific *form*: an explicit
**linear concatenated approximate shortest-common-superstring** ("pseudogenome") over **HQ reads only**
(1b), done in-memory at scale, engineered into a competitive lossless codec. That is a real but modest
**engineering + framing** delta over Leon/BdBG/ORCOM — NOT a foundational new idea, and NOT "first."

### 1b. HQ / LQ / N read separation by quality
**Mildly novel to PgRC, but not new as a concept.** Using quality to decide which reads are
trustworthy is standard (quality-guided variant calling, error correction). PgRC's specific use —
build the reference from HQ reads only, then map LQ/N onto it — is a neat design choice, not a new
primitive.

### 1c. Overlap discovery via maximal exact matches (copMEM)
**Borrowed — from the same group's own prior tool.** PgRC finds read overlaps with **copMEM**
(Grabowski & Bieniecki, Bioinformatics 2019 — sampling both sequences with coprime steps).
Overlap-by-minimizers for read compression is **ORCOM** (Grabowski et al. 2015). So the "find
overlaps to assemble/bin reads" machinery is prior art, not a PgRC invention.
- copMEM — https://dx.doi.org/10.1093/bioinformatics/bty670

### 1d. Reverse-complement folding
**Borrowed / standard.** Treating a read and its reverse complement as the same is universal in DNA
compression (GeCo family, ORCOM, SPRING, MFCompress). Not novel to PgRC.

### 1e. Approximate mapping of reads onto the reference (position + mismatches)
**Borrowed / standard.** This is exactly reference-based read alignment + edit encoding, the core of
reference-based compressors (CRAM, Quip-ref) — PgRC just uses a *self-built* reference instead of an
external genome. The reference-free twist ties back to 1a; the mapping mechanics themselves are old.

### 1f. Entropy coders (PPMd / range coder; GeCo3 in our PgRC2 variant)
**Borrowed entirely.** PPMd and range/arithmetic coding are decades old; ORCOM already used PPMd on
its streams. The GeCo3 swap in this project's PgRC2 is literally an external tool
(https://academic.oup.com/gigascience/article/9/11/giaa119/5974977). Zero entropy-coder novelty.

---

## 2. Honest verdict on PgRC's novelty
PgRC's real contribution is **1a — the pseudogenome: framing an approximate global SCS over HQ reads
as an explicit, reusable reference and building a clean lossless read codec around it (with the HQ/LQ
split, 1b).** Everything else — overlap discovery (copMEM/ORCOM), RC folding, map-and-edit encoding,
PPMd/range/GeCo3 coding — is borrowed or standard.

**So PgRC's novelty is exactly an ASSEMBLER / representation + integration novelty, NOT new entropy
math.** This confirms the user's read: *"PgRC2's novelty is not their layer but their assembler."*
A well-informed reviewer would grant PgRC the pseudogenome framing and cite copMEM, ORCOM, and
Assembltrie for the rest.

---

## 3. PgRC vs OURS — is our novelty the same kind?
| | PgRC's novelty | ARCS `--chain-pg` (ours) |
|---|---|---|
| Core idea | one global SCS-approx pseudogenome over HQ reads | **multi-contig assembly + per-contig consensus** |
| Closest prior art to the core | SCS / Assembltrie (but PgRC's framing was fresh for short reads) | **NanoSpring** (contig+consensus+encode-vs-consensus, nanopore) |
| Overlap engine | copMEM (borrowed, own group) | global k-mer index (standard) |
| Entropy coder | PPMd/range (borrowed); GeCo3 in our variant | **ARCS-DNA** (our dependency-free CM coder) + GeCo3 fallback |
| Quality coupling | none (streams independent) | **idea B: seq-surprise→quality (planned)** — the freshest angle |
| Kind of novelty | assembler + integration (not new primitives) | assembler + integration (not new primitives) |

**Conclusion.** Our novelty is the *same category* as PgRC's — an assembler/representation +
integration contribution, not a new entropy primitive. Two honest caveats keep us from overclaiming:
1. PgRC was **NOT first** to build a reference from reads (Leon 2015, BdBG 2018, ORCOM/Mince 2015,
   Assembltrie 2018 all predate it). Its contribution is a narrow engineering/framing delta — the
   linear HQ-only pseudogenome. **This actually helps us:** PgRC earned a strong paper standing on a
   crowded field of prior art, purely on a good *specific form + integration*. Our position is the
   same shape, so it is defensible on the same terms.
2. Our specific assembler (multi-contig + consensus) is **very close to NanoSpring's** (different
   read regime). So our cleanest *fresh* contribution is not the assembler alone — it's the
   **combination** (self-assembled short-read pseudogenome + our own CM coder + idea-B cross-stream
   quality). That combination, especially idea B, is where we can claim something PgRC does not do —
   see `IDEA_B_NOVELTY.md` for the honest ceiling on how novel idea B really is.

## 4. Sources
- PgRC — https://academic.oup.com/bioinformatics/article/36/7/2082/5670526 · preprint https://www.biorxiv.org/content/10.1101/710822v1
- copMEM — https://dx.doi.org/10.1093/bioinformatics/bty670 · copMEM2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10209524/
- ORCOM (minimizer read binning) — referenced via https://academic.oup.com/bioinformatics/article/34/4/558/4386919
- Assembltrie / light assembly — https://www.nature.com/articles/s41467-017-02480-6
- NanoSpring — https://www.nature.com/articles/s41598-023-29267-8
- GeCo3 — https://academic.oup.com/gigascience/article/9/11/giaa119/5974977
