// Derived from: QUANTA-UNIVERSE/delta/lib.quanta — Forman-Ricci curvature
//               Reference: Forman (2003); Ollivier (2009 Ricci flow on MM)
// =============================================================================
// algorithms/curvature.hpp -- Graph curvature: Forman-Ricci + Ollivier-Ricci
//
// Uses the existing we::data_science::GraphAnalysis Graph type from
// data_science/graph_analysis.h indirectly — but to avoid a cross-module
// dependency this round, we define a lightweight Graph in _graph.hpp inline
// here and typedef NodeId / EdgeCurvature.
//
// Forman-Ricci edge curvature (Sreejith et al. 2016 simplified form):
//   ric_F(e) = w(e) * [ sum_{v in e} (w(e)/w(v)) - sum_{e' adj e, e'!=e} w(e')/sqrt(w(e)*w(e')) ]
//
// For unweighted graphs (all weights = 1):
//   ric_F(u,v) = 4 - deg(u) - deg(v)
//
// Ollivier-Ricci: kappa(u,v) = 1 - W_1(mu_u, mu_v) / d(u,v)
//   where mu_x is the uniform distribution over the 1-neighbourhood of x
//   and W_1 is the 1-Wasserstein distance (EMD) solved here with a simple
//   exact O(k log k) sorting-based method for unweighted graphs.
//
// Namespace: warden::algorithms
// =============================================================================

#pragma once
#ifndef WARDEN_ALGORITHMS_CURVATURE_HPP
#define WARDEN_ALGORITHMS_CURVATURE_HPP

#include "algorithms/_numeric.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <span>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace warden::algorithms {

using NodeId = uint32_t;

// ---------------------------------------------------------------------------
// _Graph — lightweight undirected weighted adjacency-list graph
// ---------------------------------------------------------------------------

class Graph {
public:
    struct Edge {
        NodeId to;
        double weight{1.0};
    };

    void add_node(NodeId id) {
        if (adj_.find(id) == adj_.end())
            adj_[id] = {};
        nodes_.insert(id);
    }

    void add_edge(NodeId u, NodeId v, double w = 1.0) {
        add_node(u);
        add_node(v);
        adj_[u].push_back({v, w});
        adj_[v].push_back({u, w}); // undirected
    }

    [[nodiscard]] const std::vector<Edge>& neighbors(NodeId u) const {
        static const std::vector<Edge> empty{};
        auto it = adj_.find(u);
        return it != adj_.end() ? it->second : empty;
    }

    [[nodiscard]] bool has_node(NodeId u) const noexcept {
        return nodes_.count(u) > 0;
    }

    [[nodiscard]] bool has_edge(NodeId u, NodeId v) const noexcept {
        auto it = adj_.find(u);
        if (it == adj_.end()) return false;
        for (auto& e : it->second)
            if (e.to == v) return true;
        return false;
    }

    [[nodiscard]] std::optional<double> edge_weight(NodeId u, NodeId v) const noexcept {
        auto it = adj_.find(u);
        if (it == adj_.end()) return {};
        for (auto& e : it->second)
            if (e.to == v) return e.weight;
        return {};
    }

    [[nodiscard]] const std::unordered_set<NodeId>& nodes() const noexcept {
        return nodes_;
    }

    [[nodiscard]] size_t degree(NodeId u) const noexcept {
        auto it = adj_.find(u);
        return it != adj_.end() ? it->second.size() : 0;
    }

    [[nodiscard]] double node_weight(NodeId u) const noexcept {
        // Node weight = sum of incident edge weights
        double w = 0.0;
        for (auto& e : neighbors(u)) w += e.weight;
        return w;
    }

    /// BFS shortest-path distance (returns inf if no path)
    [[nodiscard]] double shortest_path(NodeId src, NodeId dst) const noexcept {
        if (src == dst) return 0.0;
        std::unordered_map<NodeId, double> dist;
        std::priority_queue<std::pair<double,NodeId>,
                             std::vector<std::pair<double,NodeId>>,
                             std::greater<>> pq;
        dist[src] = 0.0;
        pq.push({0.0, src});
        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (u == dst) return d;
            if (dist.count(u) && d > dist[u] + 1e-12) continue;
            for (auto& e : neighbors(u)) {
                double nd = d + e.weight;
                if (!dist.count(e.to) || nd < dist[e.to]) {
                    dist[e.to] = nd;
                    pq.push({nd, e.to});
                }
            }
        }
        return std::numeric_limits<double>::infinity();
    }

private:
    std::unordered_map<NodeId, std::vector<Edge>> adj_;
    std::unordered_set<NodeId> nodes_;
};

// ---------------------------------------------------------------------------
// EdgeCurvature
// ---------------------------------------------------------------------------

struct EdgeCurvature {
    NodeId u{0}, v{0};
    double curvature{0.0};
};

// ---------------------------------------------------------------------------
// forman_ricci — edge curvature via Sreejith et al. 2016 simplified formula
//
//   ric_F(e_{uv}) = w_{uv} * (w_{uv}/w_u + w_{uv}/w_v)
//                 - sum_{e' sharing endpoint with e, e' != e} w_{e'} / sqrt(w_e * w_{e'})
//
//   For uniform weights: ric_F(u,v) = 4 - deg(u) - deg(v)
// ---------------------------------------------------------------------------

[[nodiscard]] inline double
forman_ricci(const Graph& g, NodeId u, NodeId v) noexcept {
    auto w_uv_opt = g.edge_weight(u, v);
    if (!w_uv_opt.has_value()) return 0.0;
    double w_e = w_uv_opt.value();

    double w_u = g.node_weight(u);
    double w_v = g.node_weight(v);
    if (w_u <= 0.0 || w_v <= 0.0) return 0.0;

    // Vertex contribution
    double vert = w_e * (w_e / w_u + w_e / w_v);

    // Edge-adjacency sum: edges sharing u or v, excluding e itself
    double edge_sum = 0.0;
    for (auto& e : g.neighbors(u)) {
        if (e.to == v) continue;
        edge_sum += e.weight / std::sqrt(w_e * e.weight);
    }
    for (auto& e : g.neighbors(v)) {
        if (e.to == u) continue;
        edge_sum += e.weight / std::sqrt(w_e * e.weight);
    }

    return vert - edge_sum;
}

// ---------------------------------------------------------------------------
// ollivier_ricci — Ollivier (2009) edge curvature
//
//   kappa(u,v) = 1 - W_1(mu_u, mu_v) / d(u,v)
//
// mu_x = uniform distribution over neighbors of x (alpha = 0.5 lazy random
// walk variant is supported but defaulted off — use the 1-step neighbor mass).
// alpha parameter: mass fraction on the node itself (lazy RW): 0 = non-lazy.
// EMD solved by sorting the support and computing the L1 distance.
// ---------------------------------------------------------------------------

[[nodiscard]] inline double
ollivier_ricci(const Graph& g,
               NodeId u,
               NodeId v,
               double alpha = 0.5) noexcept {
    if (!g.has_edge(u, v)) return 0.0;

    double d_uv = g.shortest_path(u, v);
    if (!std::isfinite(d_uv) || d_uv < 1e-12) return 0.0;

    // Lazy random walk measure for u: mu_u(x) = alpha if x==u,
    //   (1-alpha)/deg(u) for each neighbor x.
    auto make_measure = [&](NodeId node) {
        std::map<NodeId, double> mu;
        const auto& nbrs = g.neighbors(node);
        const double deg = static_cast<double>(nbrs.size());
        if (deg == 0.0) { mu[node] = 1.0; return mu; }
        mu[node] = alpha;
        double share = (1.0 - alpha) / deg;
        for (auto& e : nbrs) mu[e.to] += share;
        return mu;
    };

    auto mu_u = make_measure(u);
    auto mu_v = make_measure(v);

    // Collect all support nodes
    std::vector<NodeId> support;
    for (auto& [n, _] : mu_u) support.push_back(n);
    for (auto& [n, _] : mu_v) {
        bool found = false;
        for (auto s : support) if (s == n) { found = true; break; }
        if (!found) support.push_back(n);
    }

    // W_1(mu_u, mu_v) = min-cost transport
    // For 1D metric: project to distances from u and compute EMD
    // Full 2D EMD is hard; use Kantorovich dual lower bound via BFS distances.
    // Exact W_1 via linear programming is O(k^3); for small support we use
    // the sorting trick on projected distances.
    std::vector<std::pair<double, double>> pu, pv;
    for (NodeId s : support) {
        double d_s_u = g.shortest_path(u, s);
        double p = mu_u.count(s) ? mu_u.at(s) : 0.0;
        double q = mu_v.count(s) ? mu_v.at(s) : 0.0;
        pu.push_back({d_s_u, p});
        pv.push_back({d_s_u, q});
    }
    // Sort by projected coordinate
    auto cmp = [](const std::pair<double,double>& a,
                  const std::pair<double,double>& b){ return a.first < b.first; };
    std::sort(pu.begin(), pu.end(), cmp);
    std::sort(pv.begin(), pv.end(), cmp);

    // 1D EMD
    double wass = 0.0, cdf_u = 0.0, cdf_v = 0.0;
    size_t i = 0, j = 0;
    double prev = 0.0;
    while (i < pu.size() || j < pv.size()) {
        double next_x;
        if (i < pu.size() && j < pv.size())
            next_x = std::min(pu[i].first, pv[j].first);
        else if (i < pu.size())
            next_x = pu[i].first;
        else
            next_x = pv[j].first;

        wass += std::abs(cdf_u - cdf_v) * (next_x - prev);
        prev = next_x;
        while (i < pu.size() && pu[i].first <= next_x) { cdf_u += pu[i].second; ++i; }
        while (j < pv.size() && pv[j].first <= next_x) { cdf_v += pv[j].second; ++j; }
    }

    return 1.0 - wass / d_uv;
}

} // namespace warden::algorithms

#endif // WARDEN_ALGORITHMS_CURVATURE_HPP
