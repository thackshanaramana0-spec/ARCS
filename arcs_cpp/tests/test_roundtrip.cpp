#include "../src/encoder.h"
#include "../src/decoder.h"
#include "../src/fastq_io.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>

int main() {
    // Create a small test FASTQ
    const char* tmp_in  = "arcs_test_input.fastq";
    const char* tmp_out = "arcs_test_output.arcs";
    const char* tmp_dec = "arcs_test_decoded.fastq";

    // Write test FASTQ
    {
        FILE* f = fopen(tmp_in, "w");
        assert(f);
        // 20 reads of 20bp from a simple repeating sequence
        std::string bases = "ACGTACGTACGTACGTACGT";
        for (int i = 0; i < 20; ++i) {
            // Shift start position for variation
            std::string seq = (bases + bases).substr(i % 5, 20);
            fprintf(f, "@read_%d\n%s\n+\n%s\n",
                    i, seq.c_str(), "IIIIIIIIIIIIIIIIIIII");
        }
        fclose(f);
    }

    // Compress
    ARCSParams params;
    params.k     = 5;
    params.mink  = 5;
    params.w     = 3;
    params.pilot = 10;

    ARCSEncoder enc(params);
    auto prog = enc.compress(tmp_in, tmp_out);
    printf("Compressed: %zu reads, %zu mapped, genome=%zu bp\n",
           prog.total_reads, prog.mapped_reads, prog.genome_len);

    assert(prog.total_reads == 20);

    // Decompress
    ARCSDecoder dec;
    dec.decompress(tmp_out, tmp_dec);

    // Count decoded reads
    FASTQReader rdr(tmp_dec);
    size_t n_decoded = 0;
    Read r;
    while (rdr.next(r)) ++n_decoded;
    printf("Decoded: %zu reads\n", n_decoded);
    assert(n_decoded == 20);

    // Print info
    dec.print_info(tmp_out);

    // Cleanup
    remove(tmp_in);
    remove(tmp_out);
    remove(tmp_dec);

    printf("PASS: test_roundtrip\n");
    return 0;
}
