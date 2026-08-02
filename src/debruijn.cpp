#include "debruijn.h"
#include <stack>
#include <set>
#include <algorithm>
#include <numeric>
#include <cmath>

// ── DeBruijnGraph::build ──────────────────────────────────────────────────────
void DeBruijnGraph::build(const KmerCounter& counts,
                          const KmerEncoder& enc,
                          count_t min_count) {
    k_ = enc.k();
    adj_.clear();
    n_edges_ = 0;

    // For each solid k-mer, add directed edge:
    // node_left  = k-mer[0..k-2]  (prefix, k-1 bases)
    // node_right = k-mer[1..k-1]  (suffix, k-1 bases)
    // base       = k-mer[k-1]     (last base)

    counts.foreach([&](kmer_t km, count_t cnt) {
        if (cnt < min_count) return;

        // Determine forward k-mer (vs reverse complement canonical)
        // Since we stored canonical k-mers, we need both orientations.
        // Add forward edge:
        kmer_t fwd   = km; // treat stored canonical as forward
        kmer_t rc    = enc.revcomp(km);

        // Process forward k-mer
        auto add_edge = [&](kmer_t kmer, bool is_rc) {
            // prefix = kmer >> 2 (drop last base) with k-1 bits
            kmer_t prefix = (kmer >> 2) & ((kmer_t(1) << (2*(k_-1))) - 1);
            uint8_t last  = kmer & 3;
            // suffix = kmer & ((1<<2(k-1))-1) (keep last k-1 bases)
            kmer_t suffix = kmer & ((kmer_t(1) << (2*(k_-1))) - 1);

            // Add out-edge from prefix to suffix
            auto& nd = adj_[prefix];
            nd.out.push_back({last, is_rc, cnt});
            adj_[suffix].in_deg++;
            ++n_edges_;
        };

        add_edge(fwd, false);
        if (rc != fwd) add_edge(rc, true);
    });
}

// ── DeBruijnGraph::remove_tips ────────────────────────────────────────────────
void DeBruijnGraph::remove_tips(int max_tip_len) {
    // A tip is a path starting from a node with in_deg=0 and ending at
    // a node before a branch, with total length < max_tip_len bases.
    // Iteratively remove such paths.
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<kmer_t> to_remove;

        for (auto& [node, nd] : adj_) {
            if (nd.in_deg == 0 && !nd.out.empty()) {
                // Follow this path
                kmer_t cur = node;
                int    len = 0;
                bool   is_tip = false;

                while (len < max_tip_len) {
                    auto it = adj_.find(cur);
                    if (it == adj_.end()) { is_tip = true; break; }
                    auto& cur_nd = it->second;
                    if (cur_nd.out.empty()) { is_tip = true; break; }
                    if (cur_nd.out.size() > 1) break; // joins main graph
                    // Follow single out-edge
                    uint8_t b = cur_nd.out[0].base;
                    kmer_t next = ((cur << 2) | b) & ((kmer_t(1) << (2*(k_-1)))-1);
                    to_remove.push_back(cur);
                    cur = next;
                    ++len;
                }

                if (is_tip && len < max_tip_len) {
                    for (kmer_t n : to_remove) {
                        auto it = adj_.find(n);
                        if (it != adj_.end()) {
                            n_edges_ -= it->second.out.size();
                            adj_.erase(it);
                            changed = true;
                        }
                    }
                }
                to_remove.clear();
            }
        }
    }
}

// ── Hierholzer's Eulerian path ─────────────────────────────────────────────────
// Finds an Eulerian path (or circuit) in a directed multigraph.
// Requires: at most 2 nodes with unequal in/out degree (semi-Eulerian).
// Algorithm: O(V + E) with stack-based DFS.
std::vector<kmer_t> DeBruijnGraph::hierholzer(kmer_t start) const {
    // Copy adjacency (we consume edges during traversal)
    std::unordered_map<kmer_t, std::vector<EdgeInfo>> adj_copy;
    for (const auto& [node, nd] : adj_)
        adj_copy[node] = nd.out;

    std::stack<kmer_t>  stack;
    std::vector<kmer_t> path;

    stack.push(start);

    while (!stack.empty()) {
        kmer_t v = stack.top();
        auto it = adj_copy.find(v);
        if (it != adj_copy.end() && !it->second.empty()) {
            // Take last edge (pop from back for efficiency)
            EdgeInfo e = it->second.back();
            it->second.pop_back();
            // Next node = last (k-1) bases of (v << 2 | e.base) with k-1 bits
            kmer_t next = (((v << 2) | e.base)) & ((kmer_t(1) << (2*(k_-1)))-1);
            stack.push(next);
        } else {
            path.push_back(v);
            stack.pop();
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}

// ── Path → string ──────────────────────────────────────────────────────────────
std::string DeBruijnGraph::path_to_string(const std::vector<kmer_t>& path,
                                           const KmerEncoder& enc) const {
    if (path.empty()) return "";
    // First node contributes k-1 characters, each subsequent adds 1.
    std::string s;
    s.reserve(path.size() + enc.k() - 1);

    // Decode first node (k-1 bases)
    kmer_t first = path[0];
    std::string first_str(enc.k() - 1, 'N');
    for (int i = enc.k()-2; i >= 0; --i) {
        first_str[i] = BASE_TO_CHAR[first & 3];
        first >>= 2;
    }
    s += first_str;

    // Each subsequent node adds 1 base (the last base of the node + next edge)
    // Actually: each edge in the path adds 1 base.
    // The last base of node v is determined by the edge taken from v.
    // But in our path encoding, each node is a (k-1)-mer.
    // To get the base between node[i] and node[i+1]:
    //   base = last base of k-mer = (path[i+1] & 3) but we need to check direction.
    // Simplest: last base of node[i+1] in the path order
    for (size_t i = 1; i < path.size(); ++i) {
        // The base appended = last 2 bits of path[i] in forward direction
        // Actually: path[i] is a (k-1)-mer. The overlap gives us 1 new base.
        // path[i] = last k-1 bases of the k-mer going from path[i-1] to path[i].
        // So: k-mer = (path[i-1] << 2) | last_base, last k-1 bases = path[i]
        // last_base = path[i] & ((1<<2)-1) ... wait this isn't quite right.
        // path[i] is a (k-1)-mer. Its encoding has the last base at bits[1:0].
        // The "new" base added at step i is bits [0:1] of path[i], i.e., path[i] & 3.
        // BUT: this only works if all nodes are in forward orientation.
        // Since we added both forward and RC edges, orientation may vary.
        // For simplicity: extract last base of each subsequent node.
        s += BASE_TO_CHAR[path[i] & 3]; // last base of (k-1)-mer
    }

    return s;
}

// ── find_start_nodes ──────────────────────────────────────────────────────────
std::vector<kmer_t> DeBruijnGraph::find_start_nodes() const {
    std::vector<kmer_t> starts;
    // Prefer nodes with out_degree > in_degree (Eulerian path start)
    for (const auto& [node, nd] : adj_) {
        int out = (int)nd.out.size();
        int in  = nd.in_deg;
        if (out > in) starts.push_back(node);
    }
    // If no unbalanced nodes (Eulerian circuit), start anywhere
    if (starts.empty() && !adj_.empty())
        starts.push_back(adj_.begin()->first);
    return starts;
}

// ── out_degree / in_degree ────────────────────────────────────────────────────
int DeBruijnGraph::out_degree(kmer_t node) const {
    auto it = adj_.find(node);
    return (it != adj_.end()) ? (int)it->second.out.size() : 0;
}

int DeBruijnGraph::in_degree(kmer_t node) const {
    auto it = adj_.find(node);
    return (it != adj_.end()) ? it->second.in_deg : 0;
}

// ── extract_contigs ───────────────────────────────────────────────────────────
std::vector<std::string> DeBruijnGraph::extract_contigs(
    const KmerEncoder& enc) const {

    std::vector<std::string> contigs;
    std::set<kmer_t> visited;

    // Find all nodes that start a contig:
    // - nodes with in_degree = 0
    // - nodes with in_degree > 1 (junction entry)
    // - nodes that are only reachable from themselves (in a cycle)

    for (const auto& [node, nd] : adj_) {
        bool is_start = (nd.in_deg == 0) ||
                        (nd.in_deg > 1)   ||
                        ((int)nd.out.size() > 1);
        if (!is_start) continue;
        if (visited.count(node)) continue;

        // Follow single-path forward until branch/end
        std::vector<kmer_t> path;
        path.push_back(node);
        visited.insert(node);

        kmer_t cur = node;
        while (true) {
            auto it = adj_.find(cur);
            if (it == adj_.end()) break;
            const auto& cur_nd = it->second;
            if (cur_nd.out.size() != 1) break;

            uint8_t b = cur_nd.out[0].base;
            kmer_t next = ((cur << 2) | b) & ((kmer_t(1) << (2*(k_-1)))-1);
            if (visited.count(next)) break;

            path.push_back(next);
            visited.insert(next);
            cur = next;
        }

        std::string contig = path_to_string(path, enc);
        if ((int)contig.size() >= enc.k())
            contigs.push_back(std::move(contig));
    }

    return contigs;
}

// ── eulerian_pseudogenome ──────────────────────────────────────────────────────
std::string DeBruijnGraph::eulerian_pseudogenome(const KmerEncoder& enc) const {
    if (adj_.empty()) return "";

    // For simple genomes (linear or circular), the pseudogenome is
    // the single Eulerian path through the graph.
    // For fragmented assemblies (low coverage), we concatenate contigs.

    // Strategy:
    // 1. Find all unbalanced nodes (Eulerian path start/end candidates)
    // 2. If 0 or 2 unbalanced nodes: Eulerian path exists → use Hierholzer
    // 3. If many unbalanced nodes (fragmented): use contig concatenation

    int n_unbalanced = 0;
    kmer_t start_node = adj_.begin()->first;

    for (const auto& [node, nd] : adj_) {
        int diff = (int)nd.out.size() - nd.in_deg;
        if (diff != 0) {
            ++n_unbalanced;
            if (diff > 0) start_node = node; // good start for Eulerian path
        }
    }

    if (n_unbalanced <= 2) {
        // Eulerian path exists — use Hierholzer
        auto path = hierholzer(start_node);
        return path_to_string(path, enc);
    }

    // Fragmented assembly — concatenate all contigs
    // Sort by length descending (longest first for best compression)
    auto contigs = extract_contigs(enc);
    std::sort(contigs.begin(), contigs.end(),
              [](const std::string& a, const std::string& b) {
                  return a.size() > b.size();
              });

    std::string result;
    result.reserve(n_edges_ + enc.k());
    for (const auto& c : contigs)
        result += c;  // concatenation; no separator needed for mapping
    return result;
}

// ── Greedy SCS builder (returns positions) ────────────────────────────────────
SCSResult build_greedy_scs_with_positions(const std::vector<std::string>& seqs, int k) {
    SCSResult res;
    if (seqs.empty()) return res;
    int n = (int)seqs.size();
    int rlen = (int)seqs[0].size();

    res.positions.resize(n, 0);
    res.rc.resize(n, false);
    res.n_mm.resize(n, 0);
    res.genome.reserve((size_t)n * rlen);

    // Build prefix hash: first k bases → list of read indices
    std::unordered_map<std::string, std::vector<int>> prefix_hash;
    prefix_hash.reserve(n * 2);
    for (int i = 0; i < n; ++i)
        if ((int)seqs[i].size() >= k)
            prefix_hash[seqs[i].substr(0, k)].push_back(i);

    std::vector<bool> used(n, false);

    for (int start = 0; start < n; ++start) {
        if (used[start]) continue;

        // Start a new chain
        pos_t chain_start = (pos_t)res.genome.size();
        res.genome += seqs[start];
        res.positions[start] = chain_start;
        used[start] = true;

        int max_chain = n; // safety
        while (max_chain-- > 0) {
            int genome_len = (int)res.genome.size();
            if (genome_len < k) break;
            std::string suffix = res.genome.substr(genome_len - k, k);
            auto it = prefix_hash.find(suffix);
            if (it == prefix_hash.end()) break;

            bool found = false;
            for (int next : it->second) {
                if (!used[next] && (int)seqs[next].size() > k) {
                    // Record position: this read starts at (genome_len - k)
                    res.positions[next] = (pos_t)(genome_len - k);
                    res.genome += seqs[next].substr(k);
                    used[next] = true;
                    found = true;
                    break;
                }
            }
            if (!found) break;
        }
    }
    return res;
}

std::string build_greedy_scs(const std::vector<std::string>& seqs, int k) {
    return build_greedy_scs_with_positions(seqs, k).genome;
}

// ── Full pipeline ──────────────────────────────────────────────────────────────
PseudogenomeResult build_pseudogenome(const std::vector<Read>& reads,
                                      const ARCSParams& params) {
    ARCS_CHECK(!reads.empty(), "No reads provided for pseudogenome construction");

    size_t read_len = reads[0].seq.size();

    // 1. Select k — use smaller k for short reads and low coverage
    int k = (params.k > 0) ? params.k
            : select_k(reads.size(), read_len, params.k_min, params.k_max);
    KmerEncoder enc(k);

    // 2. Count k-mers to estimate coverage
    KmerCounter counts = count_kmers(reads, enc, /*canonical=*/true);
    auto stats = counts.coverage_stats();
    double mean_cov = stats.mean_cov;

    // 3. Coverage-adaptive strategy selection
    // De Bruijn graph requires ≥5x coverage to build a useful assembly.
    // Below that, use greedy SCS which works at any coverage.
    static constexpr double DBG_COV_THRESHOLD = 5.0;
    std::string method;

    std::string genome;
    size_t n_kmers = 0;
    double est_genome = 0.0;

    if (mean_cov >= DBG_COV_THRESHOLD) {
        // High coverage: error-correct then build dBG + Eulerian path
        method = "debruijn";
        count_t thr = stats.threshold;
        std::vector<Read> corrected = reads;
        {
            ErrorCorrector ec(counts, enc, thr);
            ec.correct_batch(corrected);
        }
        // counts no longer needed after error correction
        KmerCounter counts2 = count_kmers(corrected, enc, true);
        DeBruijnGraph graph;
        graph.build(counts2, enc, thr);
        graph.remove_tips(2 * k);
        genome   = graph.eulerian_pseudogenome(enc);
        n_kmers  = graph.n_nodes();
        est_genome = (counts2.size() > 0 && mean_cov > 1.0)
            ? (double)counts2.size() * (1.0 - std::pow(1.0 - 1.0/mean_cov, (double)k))
            : (double)reads.size() * read_len;

    } else {
        // Low/medium coverage: greedy SCS. No second k-mer counting needed.
        method = "greedy_scs";
        int scs_k;
        if (mean_cov >= 2.0) {
            int avg_spacing = std::max(1, (int)(read_len / mean_cov));
            scs_k = std::max(12, (int)read_len - avg_spacing);
        } else {
            scs_k = std::max(12, (int)(read_len * 0.25));
        }
        scs_k = std::min(scs_k, (int)read_len - 5);

        std::vector<std::string> seqs;
        seqs.reserve(reads.size());
        for (const auto& r : reads) seqs.push_back(r.seq);

        auto scs = build_greedy_scs_with_positions(seqs, scs_k);
        n_kmers   = counts.size();
        est_genome = (double)reads.size() * read_len / std::max(mean_cov, 1.0);

        return { std::move(scs.genome), k, n_kmers, est_genome, method,
                 std::move(scs.positions) };
    }

    return { std::move(genome), k, n_kmers, est_genome, method, {} };
}
