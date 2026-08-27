#pragma once
#include "chain_encoder.h"

// ── Variable-order-style greedy overlap assembler (opt-in, ARCS_DBG_ASSEMBLY=1) ──
//
// Motivation (see repeat_elim.h / session notes for the full diagnosis): ARCS's
// default assembler (build_multicontig_pg) places each read in a single forward
// pass, checking only a handful of fixed seed positions against whatever's
// already been built — order-dependent, and each read is decided once and
// never reconsidered. That leaves far more redundancy in the assembled pg than
// a thorough overlap search would (measured: ~2.6x reduction from raw reads,
// vs ~7x from PgRC2's exhaustive greedy shortest-common-superstring approach).
//
// This assembler targets that gap directly, via a global greedy-overlap growth
// using a full k-mer occurrence index built once over ALL reads (not just
// "whatever's been placed so far"): to extend a contig's tail, every read
// containing that seed k-mer is checked, and whichever one has the LONGEST
// *verified* overlap wins — not just the first one found within a small cap.
//
// This is also where the "variable-order" idea (Boucher & Bowe, 2015 — using
// more context exactly where a fixed-length k-mer would be ambiguous, e.g. in
// a repeat) is realized *practically*: rather than building a formal succinct
// multi-order graph, each candidate read's own full sequence already IS the
// extra context — verifying the overlap length directly against the whole
// read naturally disambiguates repeat regions a plain fixed-k de Bruijn graph
// would collapse, with no extra data structure needed.
//
// Natively variable-length-safe (unlike PgRC2's fixed-length-only pipeline):
// nothing here assumes reads are the same size — k-mers are extracted from
// each read independently of its length.
//
// Produces the same ChainEncodeResult pg_* fields build_multicontig_pg does
// (via the same record_mapped/record_append bookkeeping), so serialization,
// the quality path, and the decoder are all unchanged — this is purely an
// alternative way to arrive at a pg + per-read placement, not a new format.
//
// call_out (default nullptr, zero overhead when omitted): same contract as
// build_multicontig_pg's CallData* parameter (see chain_encoder.h) — every
// original read gets a (contig id, contig-local position, RC flag) triple,
// covering all three populations this assembler ever produces a read from:
//   - grown reads:    place[i].cid/off directly (already global post-merge).
//   - fallback-mapped reads: inverted from their absolute pg position via
//     contig_base (binary search) back to (contig id, local offset) — these
//     land inside one of the same growth contigs, never a new one.
//   - reads that fail even the fallback pass (record_append'd): each becomes
//     its own new singleton contig, pos=0, rc=0 — same convention
//     build_multicontig_pg uses for a "new contig start" read (pl_new=1).
ChainEncodeResult build_vodbg_pg(const std::vector<Read>& reads, CallData* call_out = nullptr);
