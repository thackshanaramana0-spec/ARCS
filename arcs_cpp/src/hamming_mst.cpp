#include "hamming_mst.h"
#include <algorithm>
#include <random>
#include <unordered_map>
#include <queue>
#include <stack>
#include <cmath>
#include <numeric>

// ── Hamming+shift distance ─────────────────────────────────────────────────────
float hamming_shift_dist(const std::string& a, const std::string& b, int max_shift) {
    int la = (int)a.size(), lb = (int)b.size();
    if (la == 0 || lb == 0) return 1.0f;

    float best = 1.0f;

    for (int shift = -max_shift; shift <= max_shift; ++shift) {
        int start_a = std::max(0,  shift);
        int start_b = std::max(0, -shift);
        int len = std::min(la - start_a, lb - start_b);
        if (len <= 0) continue;

        int ham = 0;
        for (int i = 0; i < len; ++i) {
            if (a[start_a + i] != b[start_b + i]) ++ham;
        }
        float d = (float)ham / len;
        if (d < best) best = d;
    }
    return best;
}

// ── HammingLSH ────────────────────────────────────────────────────────────────
HammingLSH::HammingLSH(int read_len, int n_hashes, int bits_per_hash, uint64_t seed)
    : read_len_(read_len), n_hashes_(n_hashes), bits_(bits_per_hash) {

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> dist(0, read_len - 1);

    positions_.resize(n_hashes);
    for (int h = 0; h < n_hashes; ++h) {
        positions_[h].resize(bits_per_hash);
        for (int b = 0; b < bits_per_hash; ++b)
            positions_[h][b] = dist(rng);
    }
}

std::vector<uint64_t> HammingLSH::hash(const std::string& seq) const {
    std::vector<uint64_t> hashes(n_hashes_, 0);
    int slen = (int)seq.size();

    for (int h = 0; h < n_hashes_; ++h) {
        uint64_t hval = 0;
        for (int b = 0; b < bits_; ++b) {
            int pos = positions_[h][b];
            if (pos >= slen) continue;
            uint8_t base = encode_base(seq[pos]);
            // Use 2-bit encoding as hash bits (2 bits per position)
            hval = (hval << 2) | (base & 3);
        }
        hashes[h] = hval;
    }
    return hashes;
}

bool HammingLSH::is_candidate(const std::vector<uint64_t>& h1,
                               const std::vector<uint64_t>& h2) const {
    for (int h = 0; h < n_hashes_; ++h)
        if (h1[h] == h2[h]) return true;
    return false;
}

// ── build_knn_graph ────────────────────────────────────────────────────────────
std::vector<KNNEdge> build_knn_graph(const std::vector<std::string>& seqs,
                                      int k_neighbors,
                                      int n_lsh_hashes,
                                      int max_shift,
                                      uint64_t seed) {
    if (seqs.empty()) return {};
    int n = (int)seqs.size();
    int read_len = (int)seqs[0].size();

    HammingLSH lsh(read_len, n_lsh_hashes, 8, seed);

    // Compute all hash values
    std::vector<std::vector<uint64_t>> hashes(n);
    for (int i = 0; i < n; ++i)
        hashes[i] = lsh.hash(seqs[i]);

    // Group reads by hash bucket (per hash function)
    // For each pair sharing a bucket, compute exact distance
    std::unordered_map<uint64_t, std::vector<uint32_t>> buckets;
    std::vector<KNNEdge> candidates;

    for (int h = 0; h < n_lsh_hashes; ++h) {
        buckets.clear();
        for (int i = 0; i < n; ++i)
            buckets[hashes[i][h]].push_back((uint32_t)i);

        for (auto& [bucket_id, members] : buckets) {
            if (members.size() < 2 || members.size() > 200) continue; // skip huge buckets
            for (int a = 0; a < (int)members.size(); ++a) {
                for (int b = a+1; b < (int)members.size(); ++b) {
                    uint32_t u = members[a], v = members[b];
                    if (u > v) std::swap(u, v);
                    float d = hamming_shift_dist(seqs[u], seqs[v], max_shift);
                    candidates.push_back({u, v, d});
                }
            }
        }
    }

    // Deduplicate and keep k nearest per node
    std::sort(candidates.begin(), candidates.end(),
              [](const KNNEdge& a, const KNNEdge& b) {
                  return (a.u != b.u) ? a.u < b.u : a.dist < b.dist;
              });

    // Keep k_neighbors best edges per node
    std::vector<int> node_count(n, 0);
    std::vector<KNNEdge> result;
    result.reserve(n * k_neighbors);

    for (const auto& e : candidates) {
        if (node_count[e.u] < k_neighbors && node_count[e.v] < k_neighbors) {
            result.push_back(e);
            ++node_count[e.u];
            ++node_count[e.v];
        }
    }

    return result;
}

// ── UnionFind ─────────────────────────────────────────────────────────────────
UnionFind::UnionFind(size_t n) : parent_(n), rank_(n, 0) {
    std::iota(parent_.begin(), parent_.end(), 0);
}

int UnionFind::find(int x) {
    while (parent_[x] != x) {
        parent_[x] = parent_[parent_[x]]; // path compression
        x = parent_[x];
    }
    return x;
}

bool UnionFind::unite(int x, int y) {
    int rx = find(x), ry = find(y);
    if (rx == ry) return false;
    if (rank_[rx] < rank_[ry]) std::swap(rx, ry);
    parent_[ry] = rx;
    if (rank_[rx] == rank_[ry]) ++rank_[rx];
    return true;
}

// ── Kruskal's MST ─────────────────────────────────────────────────────────────
std::vector<KNNEdge> kruskal_mst(std::vector<KNNEdge> edges, size_t n_nodes) {
    std::sort(edges.begin(), edges.end(),
              [](const KNNEdge& a, const KNNEdge& b) { return a.dist < b.dist; });

    UnionFind uf(n_nodes);
    std::vector<KNNEdge> mst;
    mst.reserve(n_nodes - 1);

    for (const auto& e : edges) {
        if (uf.unite((int)e.u, (int)e.v)) {
            mst.push_back(e);
            if (mst.size() == n_nodes - 1) break;
        }
    }
    return mst;
}

// ── DFS order from MST ────────────────────────────────────────────────────────
std::vector<uint32_t> mst_dfs_order(const std::vector<KNNEdge>& mst, size_t n_nodes) {
    // Build adjacency list of MST
    std::vector<std::vector<uint32_t>> adj(n_nodes);
    for (const auto& e : mst) {
        adj[e.u].push_back(e.v);
        adj[e.v].push_back(e.u);
    }

    // DFS from node 0
    std::vector<bool>     visited(n_nodes, false);
    std::vector<uint32_t> order;
    order.reserve(n_nodes);

    std::stack<uint32_t> stack;
    stack.push(0);

    while (!stack.empty()) {
        uint32_t v = stack.top(); stack.pop();
        if (visited[v]) continue;
        visited[v] = true;
        order.push_back(v);
        // Push neighbors in reverse order for consistent traversal
        for (int i = (int)adj[v].size()-1; i >= 0; --i)
            if (!visited[adj[v][i]])
                stack.push(adj[v][i]);
    }

    // Append any unvisited nodes (disconnected components)
    for (uint32_t i = 0; i < (uint32_t)n_nodes; ++i)
        if (!visited[i]) order.push_back(i);

    return order;
}

// ── Full pipeline ─────────────────────────────────────────────────────────────
std::vector<uint32_t> order_unmapped_reads(const std::vector<std::string>& seqs,
                                            const ARCSParams& params) {
    int n = (int)seqs.size();
    if (n <= 1) {
        std::vector<uint32_t> trivial(n);
        std::iota(trivial.begin(), trivial.end(), 0);
        return trivial;
    }

    auto edges = build_knn_graph(seqs, params.knn, params.n_lsh);
    auto mst   = kruskal_mst(std::move(edges), n);
    return mst_dfs_order(mst, n);
}
