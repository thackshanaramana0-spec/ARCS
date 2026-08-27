#pragma once
#include "chain_encoder.h"

// ── K-mer-anchor greedy overlap assembler ("Method A", opt-in ARCS_KA_ASSEMBLY=1) ──
//
// Reimplemented as its own standalone assembler so its assembly-only cost can
// be measured in isolation (ARCS_ENC_TIMING's "[ENC] assembly:" line wraps
// whichever assembler build_chain_pg_impl dispatches to, generically), for a
// fair comparison against build_vodbg_pg ("Method B", suffix-array/APSP-based)
// and the default assembler (build_multicontig_pg).
//
// Architecture: a global K0-mer occurrence index (both views of every read),
// used to find anchor candidates for a contig's current tail/head read, then
// an X-drop (BLAST-style seed-and-extend) scan past the anchor to determine
// the accepted overlap length with bounded mismatch tolerance — NOT exact
// match only. This is the error-tolerance mechanism (as opposed to Method B's
// HQ/LQ read-filtering approach): every read participates in growth directly,
// tolerant extension absorbs sequencing errors instead of excluding them.
//
// Same bidirectional reverse-complement-mirror growth and global greedy
// (priority-queue) contig-growth loop as build_vodbg_pg, and the same
// mismatch-tolerant fallback placement pass for reads growth never touched.
// Produces the same ChainEncodeResult pg_* fields, so serialization, quality,
// decoder, and Claims 2/3 tooling are unchanged.
//
// allow_leftover_pass (default true): after the fallback mismatch-tolerant
// mapping pass, reads that STILL don't place get reassembled into their own
// small pg (a recursive call to this same function, with the flag forced
// false to bound the recursion to one extra level) instead of being appended
// raw with zero deduplication against each other — confirmed via direct
// reading of PgRC's real source (pgrc-encoder.cpp's runLQPgGeneration) that
// this is exactly what it does with reads that fail to map onto its main pg.
// Internal/recursive use only; callers should not need to pass this.
ChainEncodeResult build_kmer_anchor_pg(const std::vector<Read>& reads, bool allow_leftover_pass = true);
