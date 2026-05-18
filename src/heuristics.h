#pragma once

#include "graph.h"

#include <cstdint>
#include <vector>

namespace hs {

class SmallExactMWIS {
public:
    SmallExactMWIS(const std::vector<ll>& weight_value, const std::vector<std::vector<int>>& adj);

    SolveResult solve();

private:
    bool is_comp_edge(int u, int v) const;
    void color_sort(const std::vector<int>& candidates, std::vector<int>& order, std::vector<ll>& bound) const;
    void expand(std::vector<int>& candidates, ll current_weight);

    int n;
    int words;
    std::vector<ll> weight;
    std::vector<std::vector<ull>> comp_adj;
    std::vector<int> comp_degree;
    ll best_weight = 0;
    std::vector<int> current;
    std::vector<int> best;
};

class HeuristicSolver {
public:
    HeuristicSolver(const Component& component, uint64_t seed);

    SolveResult solve(const Timer& timer, long long budget_ms);

private:
    struct State {
        std::vector<char> in_solution;
        std::vector<int> blocked_count;
        std::vector<ll> blocked_weight;
        ll total_weight = 0;
    };

    static double key_value(ll weight, int alive_degree, double alpha, double noise_scale, double noise_unit);
    std::vector<char> construct_dynamic(double alpha, double noise_scale);
    std::vector<char> construct_forest_dp();
    void add_vertex(State& st, int v) const;
    void remove_vertex(State& st, int v) const;
    State build_state(const std::vector<char>& chosen) const;
    void greedy_augment(State& st) const;
    bool improve_heavy_insert(State& st, const Timer& timer, long long safety_ms) const;
    SolveResult local_subset_greedy(const std::vector<ll>& weights, const std::vector<std::vector<int>>& adj_local) const;
    bool improve_one_to_many(State& st, const Timer& timer, long long safety_ms) const;
    SolveResult improve_from_initial(std::vector<char> initial, const Timer& timer) const;

    const Component& comp;
    int n;
    std::vector<int> degree;
    std::vector<int> augment_order;
    SplitMix64 rng;
};

}  // namespace hs
