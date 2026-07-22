#include <stdexcept>
#include "encoder.h"
#include "decoder.h"
#include "akc.h"
#include "chain_encoder.h"
#include "caller.h"
#include "fastq_io.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <chrono>

#ifdef _WIN32
  #define SET_ENV(k, v) _putenv_s((k), (v))
#else
  #define SET_ENV(k, v) setenv((k), (v), 1)
#endif

static void print_usage(const char* prog) {
    fprintf(stderr,
        "ARCS-CPP: Dominant genomic FASTQ compressor\n"
        "  Novel algorithms: Eulerian dBG + Hamming-MST + context rANS\n\n"
        "Usage:\n"
        "  %s compress   <input.fastq[.gz]> <output.arcs>\n"
        "  %s decompress <input.arcs>       <output.fastq>\n"
        "  %s compress_pe <r1.fastq[.gz]> <r2.fastq[.gz]> <output.arcs>\n"
        "  %s info        <input.arcs>\n"
        "  %s akc         <input.fastq[.gz]>   (show AKC score only)\n"
        "  %s call        <reads.fastq[.gz]> <out.vcf>   (reference-free SNV calling)\n"
        "  %s compress --call <out.vcf> <reads.fastq[.gz]> <out.arcs>  (archive + call, one pass)\n\n"
        "Options (set before subcommand):\n"
        "  --fast           fast mode: lean FCM (5 orders) + sparse merge/place + static quality;\n"
        "                   speed/RAM optimized, always lossless; parallel encode phases preserved\n"
        "  --k <int>        k-mer size (default: auto)\n"
        "  --w <int>        minimizer window (default: 10)\n"
        "  --mink <int>     minimizer k (default: 15)\n"
        "  --knn <int>      Hamming-MST k-neighbors (default: 40)\n"
        "  --pilot <int>    pilot reads for model training (default: 5000)\n"
        "  --no-names       do not store read names\n"
        "  --gutensor       enable GU-TENSOR codon token stream\n"
        "  --lossy          lossy quality: Novaseq 4-bin (2,12,23,37)\n"
        "  --lossy8         lossy quality: uniform 8-bin\n"
        "  --lossy2         lossy quality: binary 2-bin (low/high)\n"
        "  --chunk-size N   chunked MST mode: compress N reads at a time\n"
        "                   (reduces peak RAM from O(n*L) to O(N*L); default: off)\n"
        "                   recommended: --chunk-size 5000000 for whole-genome files\n",
        prog, prog, prog, prog, prog, prog, prog
    );
}

int main(int argc, char** argv) {
    if (argc < 2) { print_usage(argv[0]); return 1; }
  try {
    ARCSParams params;
    int i = 1;

    // Parse options (before subcommand)
    while (i < argc && argv[i][0] == '-' && argv[i][1] == '-') {
        std::string opt = argv[i++];
        if (opt == "--k"      && i < argc) params.k        = std::stoi(argv[i++]);
        else if (opt == "--w" && i < argc) params.w        = std::stoi(argv[i++]);
        else if (opt == "--mink"  && i < argc) params.mink = std::stoi(argv[i++]);
        else if (opt == "--knn"   && i < argc) params.knn  = std::stoi(argv[i++]);
        else if (opt == "--pilot" && i < argc) params.pilot = std::stoi(argv[i++]);
        else if (opt == "--no-names")  params.store_names  = false;
        else if (opt == "--gutensor")  params.gutensor     = true;
        else if (opt == "--lossy")          params.quality_bins = 4;
        else if (opt == "--lossy8")         params.quality_bins = 8;
        else if (opt == "--lossy2")         params.quality_bins = 2;
        else if (opt == "--chunk-size" && i < argc) params.chunk_size = std::stoi(argv[i++]);
        else if (opt == "--no-mst")         params.use_mst      = false;
        else if (opt == "--chain")          { params.use_chain = true; params.use_mst = false; }
        else if (opt == "--chain-pg")       { params.use_chain_pg = true; params.use_chain = false; params.use_mst = false; }
        else if (opt == "--fast")           { params.use_chain_pg = true; params.use_chain = false; params.use_mst = false;
                                              // --fast: lean FCM (5 orders) + sparse merge/place index + static quality coder.
                                              // Each of these genuinely trades ratio for speed/RAM.
                                              // ARCS_ENC_NOPAR is intentionally NOT set: the three heavy phases (pg compress,
                                              // quality CM, names LZMA-9) run concurrently as async futures. On GIAB that
                                              // parallelism saves ~8-9s wall-time; serializing them would cost more than the
                                              // lean-FCM saves. NOPAR is only useful if the user needs to cap peak RAM below
                                              // the sum of all three phases — available via ARCS_ENC_NOPAR=1 env var.
                                              SET_ENV("ARCSDNA_ORDERS", "2,4,8,14,22");
                                              SET_ENV("ARCS_MERGE_SPARSE", "8"); SET_ENV("ARCS_PLACE_SPARSE", "4");
                                              SET_ENV("ARCS_QUAL_NOCM", "1"); }
        else { fprintf(stderr, "Unknown option: %s\n", opt.c_str()); return 1; }
    }

    if (i >= argc) { print_usage(argv[0]); return 1; }

    std::string cmd = argv[i++];

    if (cmd == "compress") {
        // Parse options that appear after the subcommand (e.g. ./arcs compress --chain in out)
        std::string call_vcf;   // set by --call: fused compress-and-call in one pass
        while (i < argc && argv[i][0] == '-' && argv[i][1] == '-') {
            std::string opt = argv[i++];
            if      (opt == "--chain")    { params.use_chain = true; params.use_mst = false; }
            else if (opt == "--chain-pg") { params.use_chain_pg = true; params.use_chain = false; params.use_mst = false; }
            else if (opt == "--fast")     { params.use_chain_pg = true; params.use_chain = false; params.use_mst = false;
                                            SET_ENV("ARCSDNA_ORDERS", "2,4,8,14,22");
                                            SET_ENV("ARCS_MERGE_SPARSE", "8"); SET_ENV("ARCS_PLACE_SPARSE", "4");
                                            SET_ENV("ARCS_QUAL_NOCM", "1"); }
            else if (opt == "--no-mst")   { params.use_mst   = false; }
            else if (opt == "--chunk-size" && i < argc) params.chunk_size = std::stoi(argv[i++]);
            else if (opt == "--call" && i < argc) {
                call_vcf = argv[i++];
                params.use_chain_pg = true; params.use_chain = false; params.use_mst = false;
                // LZMA-7 saves ~2.9s on names (serial pole) at +0.4% archive cost.
                // Only active in --call mode where total-time matters more than ratio.
                SET_ENV("ARCS_NAMES_LZMA_LEVEL", "7");
            }
            else { fprintf(stderr, "Unknown compress option: %s\n", opt.c_str()); return 1; }
        }
        if (i + 1 >= argc) { fprintf(stderr, "compress: need input and output\n"); return 1; }
        std::string inp = argv[i++];
        std::string out = argv[i++];

        ARCSEncoder enc(params);
        // Fused mode: capture the assembly's placements + reads so the same single
        // assembly that builds the archive also drives reference-free calling.
        CallData fused_cd; std::vector<Read> fused_reads;
        int call_result = 0; bool calling_done = false;
        if (!call_vcf.empty()) {
            enc.call_capture_ = &fused_cd; enc.call_reads_ = &fused_reads;
            // Hook runs as 3rd async task overlapping quality+names encoding (Win 1).
            // fused_cd and fused_reads are both ready before the hook is launched.
            enc.post_assembly_call_hook_ = [&]() {
                if (!fused_cd.valid) { call_result = -1; calling_done = true; return; }
                fprintf(stderr, "[ARCS] fused calling from the compression assembly...\n");
                call_result = run_variant_call(fused_reads, fused_cd, call_vcf);
                calling_done = true;
            };
        }
        enc.set_progress_callback([](const EncodeProgress& p) {
            fprintf(stderr,
                "[ARCS] Reads: %zu | Mapped: %.1f%% | G*: %.1f MB | k: %d | AKC: %.4f | Method: %s | Time: %.1fs\n",
                p.total_reads,
                p.total_reads > 0 ? 100.0 * p.mapped_reads / p.total_reads : 0.0,
                p.genome_len / 1e6,
                p.k_used,
                p.akc_score,
                p.pg_method.c_str(),
                p.elapsed_sec
            );
        });

        auto prog = enc.compress(inp, out);

        // Print compression ratio vs uncompressed FASTQ size
        // Uncompressed = n_reads × (read_len + name_len + 4 newlines + '+')
        // Approximate: prog.total_reads × avg_record_size
        // Better: use gzip to get uncompressed size
        {
            // Get compressed output size
            FILE* ouf = fopen(out.c_str(), "rb");
            long comp_sz = 0;
            if (ouf) { fseek(ouf, 0, SEEK_END); comp_sz = ftell(ouf); fclose(ouf); }

            // Estimate uncompressed size: for .gz input, uncompress to count
            // Simple approach: read first 1000 reads to get avg record size
            double avg_record = 0;
            {
                FASTQReader rdr(inp);
                Read r; size_t cnt = 0;
                while (cnt < 1000 && rdr.next(r)) {
                    avg_record += 4 + r.name.size() + r.seq.size() * 2 + 3;
                    ++cnt;
                }
                if (cnt > 0) avg_record /= cnt;
            }
            long raw_approx = (long)(prog.total_reads * avg_record);

            fprintf(stderr,
                "[ARCS] Uncompressed: ~%.3f MB | Archive: %.3f MB | Ratio: ~%.2fx\n",
                raw_approx/1e6, comp_sz/1e6,
                comp_sz > 0 ? (double)raw_approx/comp_sz : 0.0);
        }
        // Fused calling: if the hook ran (parallel path), result is already in
        // call_result. If ARCS_ENC_NOPAR was set the hook never ran — do it now.
        if (!call_vcf.empty()) {
            if (!calling_done) {
                if (!fused_cd.valid) { fprintf(stderr, "arcs: --call requires chain-pg assembly (auto-set); no placements captured\n"); return 1; }
                fprintf(stderr, "[ARCS] fused calling from the compression assembly...\n");
                call_result = run_variant_call(fused_reads, fused_cd, call_vcf);
            }
            if (call_result < 0) return 1;
        }
        return 0;

    } else if (cmd == "decompress") {
        if (i + 1 >= argc) { fprintf(stderr, "decompress: need input and output\n"); return 1; }
        std::string inp = argv[i++];
        std::string out = argv[i++];

        ARCSDecoder dec;
        auto t0 = std::chrono::steady_clock::now();
        dec.decompress(inp, out);
        auto t1 = std::chrono::steady_clock::now();
        double secs = std::chrono::duration<double>(t1-t0).count();
        fprintf(stderr, "[ARCS] Decompressed in %.2fs\n", secs);
        return 0;

    } else if (cmd == "compress_pe") {
        if (i + 2 >= argc) { fprintf(stderr, "compress_pe: need r1, r2, output\n"); return 1; }
        std::string r1  = argv[i++];
        std::string r2  = argv[i++];
        std::string out = argv[i++];
        params.paired_end = true;

        ARCSEncoder enc(params);
        enc.compress_pe(r1, r2, out);
        return 0;

    } else if (cmd == "info") {
        if (i >= argc) { fprintf(stderr, "info: need input file\n"); return 1; }
        ARCSDecoder dec;
        dec.print_info(argv[i++]);
        return 0;

    } else if (cmd == "akc") {
        if (i >= argc) { fprintf(stderr, "akc: need input file\n"); return 1; }
        AKCRouter router;
        auto result = router.compute_from_file(argv[i++], 1000);
        static const char* names[] = {"AMPLICON","TARGETED","WGS","METAGENOMIC"};
        printf("AKC score: %.6f\nRegime: %s\nUse BSC: %s\n",
               result.score,
               names[(int)result.regime],
               result.use_bsc ? "yes" : "no");
        return 0;

    } else if (cmd == "call") {
        // Reference-free SNV calling: assemble reads → build pileup from ARCS's own
        // placements → internal k-mer filter → VCF (contig coords). No external
        // aligner, no external k-mer tool.
        if (i + 1 >= argc) { fprintf(stderr, "call: need <reads.fastq[.gz]> <out.vcf>\n"); return 1; }
        std::string reads_path = argv[i++];
        std::string vcf_path   = argv[i++];
        // Use the COMPRESSION-default assembly (tight placement MAXMM=4/OVERR=0.06 +
        // merge-on). Measured to call BEST (F1 ~0.95 across 3 regions incl held-out,
        // robust to exact-position matching) — the loose config the old bwa-pileup
        // pipeline used force-places paralog/repeat reads and injects noisy columns.
        // This also means compression and calling share ONE assembly config (enables
        // the fused compress-and-call one-pass mode). Env vars still override.
        std::vector<Read> reads;
        { FASTQReader rdr(reads_path); Read r; while (rdr.next(r)) reads.push_back(std::move(r)); }
        fprintf(stderr, "[CALL] loaded %zu reads; assembling...\n", reads.size());
        CallData cd;
        (void)build_multicontig_pg(reads, &cd);
        int nc = run_variant_call(reads, cd, vcf_path);
        return nc < 0 ? 1 : 0;

    } else {
        fprintf(stderr, "Unknown command: %s\n", cmd.c_str());
        print_usage(argv[0]);
        return 1;
    }
  } catch (const std::exception& e) {
    fprintf(stderr, "arcs: error: %s\n", e.what());
    return 1;
  }
}
