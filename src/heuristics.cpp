#include "heuristics.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <queue>
#include <utility>

namespace hs {

SmallExactMWIS::SmallExactMWIS(const std::vector<ll>& weight_value, const std::vector<std::vector<int>>& adj)
    : n(static_cast<int>(weight_value.size())),
      words((n + 63) / 64),
      weight(weight_value),
      comp_adj(n, std::vector<ull>(words, 0)),
      comp_degree(n, 0) {
    for (int u = 0; u < n; ++u) {
        for (int word = 0; word < words; ++word) {
            comp_adj[u][word] = ~0ULL;
        }
    }

    if (words > 0) {
        int last_bits = n % 64;
        if (last_bits != 0) {
            ull mask = (1ULL << last_bits) - 1ULL;
            for (int u = 0; u < n; ++u) {
                comp_adj[u][words - 1] &= mask;
            }
        }
    }

    for (int u = 0; u < n; ++u) {
        comp_adj[u][u >> 6] &= ~(1ULL << (u & 63));
    }

    for (int u = 0; u < n; ++u) {
        for (int v : adj[u]) {
            if (u < v) {
                comp_adj[u][v >> 6] &= ~(1ULL << (v & 63));
                comp_adj[v][u >> 6] &= ~(1ULL << (u & 63));
            }
        }
    }

    for (int u = 0; u < n; ++u) {
        int degree_count = 0;
        for (int word = 0; word < words; ++word) {
            degree_count += __builtin_popcountll(comp_adj[u][word]);
        }
        comp_degree[u] = degree_count;
    }
}

bool SmallExactMWIS::is_comp_edge(int u, int v) const {
    return ((comp_adj[u][v >> 6] >> (v & 63)) & 1ULL) != 0ULL;
}

void SmallExactMWIS::color_sort(const std::vector<int>& candidates,
                                std::vector<int>& order,
                                std::vector<ll>& bound) const {
    std::vector<std::vector<int>> color_vertices;
    std::vector<std::vector<ull>> color_bits;
    std::vector<ll> color_max;

    order.clear();
    bound.clear();

    for (int v : candidates) {
        int color = 0;
        for (; color < static_cast<int>(color_vertices.size()); ++color) {
            bool conflict = false;
            for (int word = 0; word < words; ++word) {
                if ((comp_adj[v][word] & color_bits[color][word]) != 0ULL) {
                    conflict = true;
                    break;
                }
            }
            if (!conflict) {
                break;
            }
        }

        if (color == static_cast<int>(color_vertices.size())) {
            color_vertices.push_back({});
            color_bits.push_back(std::vector<ull>(words, 0ULL));
            color_max.push_back(0);
        }

        color_vertices[color].push_back(v);
        color_bits[color][v >> 6] |= (1ULL << (v & 63));
        color_max[color] = std::max(color_max[color], weight[v]);
    }

    ll prefix = 0;
    for (int color = 0; color < static_cast<int>(color_vertices.size()); ++color) {
        auto& verts = color_vertices[color];
        std::sort(verts.begin(), verts.end(), [&](int a, int b) {
            if (weight[a] != weight[b]) {
                return weight[a] < weight[b];
            }
            return comp_degree[a] < comp_degree[b];
        });
        prefix += color_max[color];
        for (int v : verts) {
            order.push_back(v);
            bound.push_back(prefix);
        }
    }
}

void SmallExactMWIS::expand(std::vector<int>& candidates, ll current_weight) {
    if (candidates.empty()) {
        if (current_weight > best_weight) {
            best_weight = current_weight;
            best = current;
        }
        return;
    }

    std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        if (weight[a] != weight[b]) {
            return weight[a] > weight[b];
        }
        return comp_degree[a] > comp_degree[b];
    });

    std::vector<int> ordered;
    std::vector<ll> bound;
    color_sort(candidates, ordered, bound);

    for (int idx = static_cast<int>(ordered.size()) - 1; idx >= 0; --idx) {
        if (current_weight + bound[idx] <= best_weight) {
            return;
        }

        int v = ordered[idx];
        current.push_back(v);

        std::vector<int> next;
        next.reserve(idx);
        for (int j = 0; j < idx; ++j) {
            int u = ordered[j];
            if (is_comp_edge(v, u)) {
                next.push_back(u);
            }
        }

        ll new_weight = current_weight + weight[v];
        if (next.empty()) {
            if (new_weight > best_weight) {
                best_weight = new_weight;
                best = current;
            }
        } else {
            expand(next, new_weight);
        }

        current.pop_back();
    }
}

SolveResult SmallExactMWIS::solve() {
    std::vector<int> all(n);
    std::iota(all.begin(), all.end(), 0);
    expand(all, 0);

    SolveResult result;
    result.weight = best_weight;
    result.chosen.assign(n, 0);
    for (int v : best) {
        result.chosen[v] = 1;
    }
    return result;
}

HeuristicSolver::HeuristicSolver(const Component& component, uint64_t seed)
    : comp(component),
      n(static_cast<int>(component.weight.size())),
      degree(n),
      augment_order(n),
      rng(seed) {
    for (int i = 0; i < n; ++i) {
        degree[i] = static_cast<int>(comp.adj[i].size());
        augment_order[i] = i;
    }

    std::sort(augment_order.begin(), augment_order.end(), [&](int a, int b) {
        long double score_a =
            static_cast<long double>(comp.weight[a]) / std::pow(static_cast<long double>(degree[a] + 1), 0.7L);
        long double score_b =
            static_cast<long double>(comp.weight[b]) / std::pow(static_cast<long double>(degree[b] + 1), 0.7L);
        if (score_a != score_b) {
            return score_a > score_b;
        }
        return comp.weight[a] > comp.weight[b];
    });
}

double HeuristicSolver::key_value(ll weight,
                                  int alive_degree,
                                  double alpha,
                                  double noise_scale,
                                  double noise_unit) {
    double perturbed = static_cast<double>(weight) * (1.0 + noise_scale * (noise_unit - 0.5));
    return perturbed / std::pow(static_cast<double>(alive_degree + 1), alpha);
}

std::vector<char> HeuristicSolver::construct_dynamic(double alpha, double noise_scale) {
    struct Node {
        double key;
        int v;
        int snapshot_degree;

        bool operator<(const Node& other) const { return key < other.key; }
    };

    std::vector<int> alive_degree(n);
    std::vector<char> status(n, 0);
    std::vector<double> noise(n, 0.5);

    for (int v = 0; v < n; ++v) {
        alive_degree[v] = degree[v];
        noise[v] = rng.next_double();
    }

    auto push_vertex = [&](std::priority_queue<Node>& pq, int v) {
        pq.push(Node{key_value(comp.weight[v], alive_degree[v], alpha, noise_scale, noise[v]), v, alive_degree[v]});
    };

    std::priority_queue<Node> pq;
    for (int v = 0; v < n; ++v) {
        push_vertex(pq, v);
    }

    auto erase_from_alive = [&](int v) {
        for (int to : comp.adj[v]) {
            if (status[to] == 0) {
                --alive_degree[to];
                push_vertex(pq, to);
            }
        }
    };

    while (!pq.empty()) {
        Node node = pq.top();
        pq.pop();
        int v = node.v;
        if (status[v] != 0 || node.snapshot_degree != alive_degree[v]) {
            continue;
        }

        status[v] = 1;
        erase_from_alive(v);

        for (int to : comp.adj[v]) {
            if (status[to] != 0) {
                continue;
            }
            status[to] = 2;
            erase_from_alive(to);
        }
    }

    std::vector<char> chosen(n, 0);
    for (int v = 0; v < n; ++v) {
        if (status[v] == 1) {
            chosen[v] = 1;
        }
    }
    return chosen;
}

std::vector<char> HeuristicSolver::construct_forest_dp() {
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    for (int i = n - 1; i > 0; --i) {
        int j = static_cast<int>(rng.next_u64() % static_cast<uint64_t>(i + 1));
        std::swap(perm[i], perm[j]);
    }

    std::vector<char> seen(n, 0);
    std::vector<char> in_forest(n, 0);
    std::vector<int> forest_vertices;
    forest_vertices.reserve(n);

    for (int v : perm) {
        int seen_neighbors = 0;
        for (int to : comp.adj[v]) {
            if (seen[to]) {
                ++seen_neighbors;
            }
        }
        if (seen_neighbors <= 1) {
            in_forest[v] = 1;
            forest_vertices.push_back(v);
        }
        seen[v] = 1;
    }

    std::vector<int> parent(n, -1);
    std::vector<int> order;
    order.reserve(forest_vertices.size());
    std::vector<int> roots;
    roots.reserve(forest_vertices.size());

    for (int root : forest_vertices) {
        if (parent[root] != -1) {
            continue;
        }
        parent[root] = root;
        roots.push_back(root);
        std::vector<int> stack;
        stack.push_back(root);
        while (!stack.empty()) {
            int v = stack.back();
            stack.pop_back();
            order.push_back(v);
            for (int to : comp.adj[v]) {
                if (!in_forest[to] || parent[to] != -1) {
                    continue;
                }
                parent[to] = v;
                stack.push_back(to);
            }
        }
    }

    std::vector<ll> take(n, 0);
    std::vector<ll> skip(n, 0);
    for (int idx = static_cast<int>(order.size()) - 1; idx >= 0; --idx) {
        int v = order[idx];
        take[v] = comp.weight[v];
        skip[v] = 0;
        for (int to : comp.adj[v]) {
            if (in_forest[to] && parent[to] == v) {
                take[v] += skip[to];
                skip[v] += std::max(take[to], skip[to]);
            }
        }
    }

    std::vector<char> chosen(n, 0);
    std::vector<std::pair<int, bool>> stack;
    stack.reserve(forest_vertices.size());
    for (int root : roots) {
        stack.push_back({root, false});
    }

    while (!stack.empty()) {
        std::pair<int, bool> state = stack.back();
        stack.pop_back();
        int v = state.first;
        bool parent_taken = state.second;
        bool take_v = false;
        if (!parent_taken && take[v] >= skip[v]) {
            take_v = true;
            chosen[v] = 1;
        }
        for (int to : comp.adj[v]) {
            if (in_forest[to] && parent[to] == v) {
                stack.push_back({to, take_v});
            }
        }
    }

    return chosen;
}

void HeuristicSolver::add_vertex(State& st, int v) const {
    st.in_solution[v] = 1;
    st.total_weight += comp.weight[v];
    for (int to : comp.adj[v]) {
        ++st.blocked_count[to];
        st.blocked_weight[to] += comp.weight[v];
    }
}

void HeuristicSolver::remove_vertex(State& st, int v) const {
    st.in_solution[v] = 0;
    st.total_weight -= comp.weight[v];
    for (int to : comp.adj[v]) {
        --st.blocked_count[to];
        st.blocked_weight[to] -= comp.weight[v];
    }
}

HeuristicSolver::State HeuristicSolver::build_state(const std::vector<char>& chosen) const {
    State st;
    st.in_solution.assign(n, 0);
    st.blocked_count.assign(n, 0);
    st.blocked_weight.assign(n, 0);
    st.total_weight = 0;

    for (int v = 0; v < n; ++v) {
        if (chosen[v]) {
            add_vertex(st, v);
        }
    }
    return st;
}

void HeuristicSolver::greedy_augment(State& st) const {
    for (int v : augment_order) {
        if (!st.in_solution[v] && st.blocked_count[v] == 0) {
            add_vertex(st, v);
        }
    }
}

bool HeuristicSolver::improve_heavy_insert(State& st, const Timer& timer, long long safety_ms) const {
    std::vector<std::pair<ll, int>> candidates;
    candidates.reserve(n);
    for (int v = 0; v < n; ++v) {
        if (!st.in_solution[v] && st.blocked_count[v] > 0 && st.blocked_weight[v] < comp.weight[v]) {
            candidates.push_back({comp.weight[v] - st.blocked_weight[v], v});
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) {
            return a.first > b.first;
        }
        return a.second < b.second;
    });

    bool changed = false;
    std::vector<int> to_remove;

    for (const auto& entry : candidates) {
        if (timer.expired(safety_ms)) {
            break;
        }

        int v = entry.second;
        if (st.in_solution[v] || st.blocked_count[v] == 0 || st.blocked_weight[v] >= comp.weight[v]) {
            continue;
        }

        to_remove.clear();
        for (int to : comp.adj[v]) {
            if (st.in_solution[to]) {
                to_remove.push_back(to);
            }
        }

        ll removed_weight = 0;
        for (int u : to_remove) {
            removed_weight += comp.weight[u];
        }
        if (removed_weight >= comp.weight[v]) {
            continue;
        }

        for (int u : to_remove) {
            remove_vertex(st, u);
        }
        if (st.blocked_count[v] == 0) {
            add_vertex(st, v);
            changed = true;
        }
    }

    return changed;
}

SolveResult HeuristicSolver::local_subset_greedy(const std::vector<ll>& weights,
                                                 const std::vector<std::vector<int>>& adj_local) const {
    int k = static_cast<int>(weights.size());
    std::vector<int> deg_local(k);
    std::vector<int> order(k);
    std::iota(order.begin(), order.end(), 0);

    for (int i = 0; i < k; ++i) {
        deg_local[i] = static_cast<int>(adj_local[i].size());
    }

    auto solve_with_alpha = [&](double alpha) {
        std::vector<int> perm = order;
        std::sort(perm.begin(), perm.end(), [&](int a, int b) {
            long double score_a =
                static_cast<long double>(weights[a]) / std::pow(static_cast<long double>(deg_local[a] + 1), alpha);
            long double score_b =
                static_cast<long double>(weights[b]) / std::pow(static_cast<long double>(deg_local[b] + 1), alpha);
            if (score_a != score_b) {
                return score_a > score_b;
            }
            return weights[a] > weights[b];
        });

        std::vector<char> chosen(k, 0);
        std::vector<char> banned(k, 0);
        ll total = 0;
        for (int v : perm) {
            if (banned[v]) {
                continue;
            }
            chosen[v] = 1;
            total += weights[v];
            banned[v] = 1;
            for (int to : adj_local[v]) {
                banned[to] = 1;
            }
        }
        return SolveResult{total, chosen};
    };

    SolveResult best = solve_with_alpha(0.0);
    for (double alpha : {0.4, 0.8, 1.2}) {
        SolveResult cand = solve_with_alpha(alpha);
        if (cand.weight > best.weight) {
            best = std::move(cand);
        }
    }
    return best;
}

bool HeuristicSolver::improve_one_to_many(State& st, const Timer& timer, long long safety_ms) const {
    std::vector<int> selected_vertices;
    selected_vertices.reserve(n);
    for (int v = 0; v < n; ++v) {
        if (st.in_solution[v]) {
            selected_vertices.push_back(v);
        }
    }

    std::sort(selected_vertices.begin(), selected_vertices.end(), [&](int a, int b) {
        if (comp.weight[a] != comp.weight[b]) {
            return comp.weight[a] < comp.weight[b];
        }
        return degree[a] < degree[b];
    });

    bool changed = false;

    for (int s : selected_vertices) {
        if (timer.expired(safety_ms)) {
            break;
        }
        if (!st.in_solution[s]) {
            continue;
        }

        std::vector<int> candidates;
        ll candidate_weight_sum = 0;
        for (int to : comp.adj[s]) {
            if (!st.in_solution[to] && st.blocked_count[to] == 1) {
                candidates.push_back(to);
                candidate_weight_sum += comp.weight[to];
            }
        }

        if (candidates.empty() || candidate_weight_sum <= comp.weight[s]) {
            continue;
        }

        int k = static_cast<int>(candidates.size());
        std::vector<int> pos(n, -1);
        std::vector<ll> local_weights(k);
        std::vector<std::vector<int>> local_adj(k);
        ll edge_count = 0;

        for (int i = 0; i < k; ++i) {
            pos[candidates[i]] = i;
            local_weights[i] = comp.weight[candidates[i]];
        }

        for (int i = 0; i < k; ++i) {
            int u = candidates[i];
            for (int to : comp.adj[u]) {
                int j = (to >= 0 && to < n) ? pos[to] : -1;
                if (j != -1 && i < j) {
                    local_adj[i].push_back(j);
                    local_adj[j].push_back(i);
                    ++edge_count;
                }
            }
        }

        for (int u : candidates) {
            pos[u] = -1;
        }

        SolveResult local_best;
        if (edge_count == 0) {
            local_best.weight = candidate_weight_sum;
            local_best.chosen.assign(k, 1);
        } else if (k <= 44) {
            SmallExactMWIS exact(local_weights, local_adj);
            local_best = exact.solve();
        } else {
            local_best = local_subset_greedy(local_weights, local_adj);
        }

        if (local_best.weight <= comp.weight[s]) {
            continue;
        }

        remove_vertex(st, s);
        for (int i = 0; i < k; ++i) {
            if (local_best.chosen[i] && st.blocked_count[candidates[i]] == 0) {
                add_vertex(st, candidates[i]);
            }
        }
        changed = true;
    }

    return changed;
}

SolveResult HeuristicSolver::improve_from_initial(std::vector<char> initial, const Timer& timer) const {
    State st = build_state(initial);
    greedy_augment(st);

    for (int iter = 0; iter < 8 && !timer.expired(150); ++iter) {
        bool changed = false;
        changed |= improve_heavy_insert(st, timer, 120);
        greedy_augment(st);
        changed |= improve_one_to_many(st, timer, 120);
        greedy_augment(st);
        changed |= improve_heavy_insert(st, timer, 120);
        greedy_augment(st);
        if (!changed) {
            break;
        }
    }

    return SolveResult{st.total_weight, st.in_solution};
}

SolveResult HeuristicSolver::solve(const Timer& timer, long long budget_ms) {
    ll best_weight = -1;
    std::vector<char> best_chosen;
    int stale_runs = 0;
    int runs = 0;
    const int min_runs = 6;
    const int stale_limit = 20;
    long long start_ms = timer.elapsed_ms();

    const std::array<double, 6> alphas{0.0, 0.35, 0.6, 0.85, 1.1, 1.35};
    const std::array<double, 4> noises{0.0, 0.04, 0.08, 0.12};

    if (budget_ms >= 400 && n <= 30000 && !timer.expired(150)) {
        SolveResult forest_seed = improve_from_initial(construct_forest_dp(), timer);
        best_weight = forest_seed.weight;
        best_chosen = std::move(forest_seed.chosen);
    }

    while (!timer.expired(150)) {
        if (timer.elapsed_ms() - start_ms >= budget_ms) {
            break;
        }

        double alpha = alphas[runs % static_cast<int>(alphas.size())];
        double noise = noises[(runs / static_cast<int>(alphas.size())) % static_cast<int>(noises.size())];

        SolveResult candidate = improve_from_initial(construct_dynamic(alpha, noise), timer);

        if (candidate.weight > best_weight) {
            best_weight = candidate.weight;
            best_chosen = std::move(candidate.chosen);
            stale_runs = 0;
        } else {
            ++stale_runs;
        }

        ++runs;
        if (runs >= min_runs && stale_runs >= stale_limit) {
            break;
        }
    }

    if (best_weight < 0) {
        best_weight = 0;
        best_chosen.assign(n, 0);
    }

    return {best_weight, best_chosen};
}

}  // namespace hs
