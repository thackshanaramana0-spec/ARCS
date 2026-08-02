#pragma once
#include <string>
#include <vector>
#include "common.h"        // Read

struct CallData;           // defined in chain_encoder.h

// Reference-free SNV caller. Builds a read→consensus pileup from the placements
// the assembler already computed (CallData), counts k-mers internally (no external
// counter), applies the frozen bubble-cleanliness filters, and writes a VCF in
// contig coordinates. No external aligner, no external k-mer tool. Returns the
// number of variants called.
int run_variant_call(const std::vector<Read>& reads, const CallData& cd,
                     const std::string& out_vcf);
