#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace hs {

using ll = long long;
using ull = unsigned long long;

constexpr ll NEG_INF = std::numeric_limits<ll>::min() / 4;

struct SolveResult {
    ll weight = 0;
    std::vector<char> chosen;
};

struct Timer {
    std::chrono::steady_clock::time_point start;
    long long limit_ms;

    explicit Timer(long long limit_ms_value);

    long long elapsed_ms() const;
    long long remaining_ms() const;
    bool expired(long long safety_ms = 0) const;
};

struct SplitMix64 {
    uint64_t state;

    explicit SplitMix64(uint64_t seed);

    uint64_t next_u64();
    double next_double();
};

long long read_time_limit_ms(int argc, char** argv);

struct Component {
    std::vector<int> global;
    std::vector<ll> weight;
    std::vector<std::vector<int>> adj;
    int edges = 0;
};

struct GraphData {
    std::vector<int> degree;
    std::vector<std::vector<int>> graph;
    std::vector<std::vector<int>> component_vertices;
};

GraphData build_graph_data(int n, const std::vector<std::pair<int, int>>& edges);

Component build_component(const std::vector<int>& vertices,
                          const std::vector<ll>& weight_global,
                          const std::vector<std::vector<int>>& graph,
                          std::vector<int>& local_index);

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

struct ProblemInstance {
    int n = 0;
    int m = 0;
    std::vector<ll> weight;
    std::vector<std::pair<int, int>> edges;
};

struct FinalSolution {
    ll total_weight = 0;
    std::vector<int> picked;
};

FinalSolution solve_problem(const ProblemInstance& instance, long long time_limit_ms);

Timer::Timer(long long limit_ms_value)
    : start(std::chrono::steady_clock::now()), limit_ms(limit_ms_value) {}

long long Timer::elapsed_ms() const {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
}

long long Timer::remaining_ms() const {
    long long remaining = limit_ms - elapsed_ms();
    return remaining > 0 ? remaining : 0;
}

bool Timer::expired(long long safety_ms) const { return elapsed_ms() + safety_ms >= limit_ms; }

SplitMix64::SplitMix64(uint64_t seed) : state(seed) {}

uint64_t SplitMix64::next_u64() {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31U);
}

double SplitMix64::next_double() {
    return static_cast<double>(next_u64() >> 11U) * (1.0 / static_cast<double>(1ULL << 53U));
}

long long read_time_limit_ms(int argc, char** argv) {
    const char* env = std::getenv("HS_TIME_LIMIT_MS");
    if (env != nullptr) {
        return std::max(1000LL, std::atoll(env));
    }
    if (argc >= 2) {
        return std::max(1000LL, std::atoll(argv[1]));
    }
    return 290000LL;
}

GraphData build_graph_data(int n, const std::vector<std::pair<int, int>>& edges) {
    GraphData data;
    data.degree.assign(n, 0);

    for (const auto& edge : edges) {
        ++data.degree[edge.first];
        ++data.degree[edge.second];
    }

    data.graph.assign(n, {});
    for (int i = 0; i < n; ++i) {
        data.graph[i].reserve(data.degree[i]);
    }
    for (const auto& edge : edges) {
        int u = edge.first;
        int v = edge.second;
        data.graph[u].push_back(v);
        data.graph[v].push_back(u);
    }

    data.component_vertices.reserve(n);
    std::vector<char> seen(n, 0);

    for (int start = 0; start < n; ++start) {
        if (seen[start]) {
            continue;
        }

        std::vector<int> vertices;
        std::queue<int> q;
        q.push(start);
        seen[start] = 1;

        while (!q.empty()) {
            int v = q.front();
            q.pop();
            vertices.push_back(v);
            for (int to : data.graph[v]) {
                if (!seen[to]) {
                    seen[to] = 1;
                    q.push(to);
                }
            }
        }

        data.component_vertices.push_back(std::move(vertices));
    }

    return data;
}

Component build_component(const std::vector<int>& vertices,
                          const std::vector<ll>& weight_global,
                          const std::vector<std::vector<int>>& graph,
                          std::vector<int>& local_index) {
    Component comp;
    comp.global = vertices;

    int n = static_cast<int>(vertices.size());
    comp.weight.resize(n);
    comp.adj.assign(n, {});

    for (int i = 0; i < n; ++i) {
        local_index[vertices[i]] = i;
        comp.weight[i] = weight_global[vertices[i]];
    }

    int edges = 0;
    for (int i = 0; i < n; ++i) {
        int gv = vertices[i];
        for (int to_global : graph[gv]) {
            int j = local_index[to_global];
            if (j != -1) {
                comp.adj[i].push_back(j);
                if (i < j) {
                    ++edges;
                }
            }
        }
    }

    for (int v : vertices) {
        local_index[v] = -1;
    }

    comp.edges = edges;
    return comp;
}

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

namespace {

struct Dinic {
    struct Edge {
        int to;
        int rev;
        ll cap;
    };

    int n;
    std::vector<std::vector<Edge>> graph;
    std::vector<int> level;
    std::vector<int> ptr;

    explicit Dinic(int n_value) : n(n_value), graph(n_value), level(n_value), ptr(n_value) {}

    void add_edge(int u, int v, ll cap) {
        Edge a{v, static_cast<int>(graph[v].size()), cap};
        Edge b{u, static_cast<int>(graph[u].size()), 0};
        graph[u].push_back(a);
        graph[v].push_back(b);
    }

    bool bfs(int source, int sink) {
        std::fill(level.begin(), level.end(), -1);
        std::queue<int> q;
        level[source] = 0;
        q.push(source);
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (const Edge& e : graph[v]) {
                if (e.cap > 0 && level[e.to] == -1) {
                    level[e.to] = level[v] + 1;
                    q.push(e.to);
                }
            }
        }
        return level[sink] != -1;
    }

    ll dfs(int v, int sink, ll pushed) {
        if (pushed == 0) {
            return 0;
        }
        if (v == sink) {
            return pushed;
        }
        for (int& cid = ptr[v]; cid < static_cast<int>(graph[v].size()); ++cid) {
            Edge& e = graph[v][cid];
            if (e.cap <= 0 || level[e.to] != level[v] + 1) {
                continue;
            }
            ll tr = dfs(e.to, sink, std::min(pushed, e.cap));
            if (tr == 0) {
                continue;
            }
            e.cap -= tr;
            graph[e.to][e.rev].cap += tr;
            return tr;
        }
        return 0;
    }

    ll max_flow(int source, int sink) {
        ll flow = 0;
        while (bfs(source, sink)) {
            std::fill(ptr.begin(), ptr.end(), 0);
            while (true) {
                ll pushed = dfs(source, sink, std::numeric_limits<ll>::max() / 4);
                if (pushed == 0) {
                    break;
                }
                flow += pushed;
            }
        }
        return flow;
    }

    std::vector<char> reachable_from_source(int source) const {
        std::vector<char> seen(n, 0);
        std::queue<int> q;
        q.push(source);
        seen[source] = 1;
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (const Edge& e : graph[v]) {
                if (e.cap > 0 && !seen[e.to]) {
                    seen[e.to] = 1;
                    q.push(e.to);
                }
            }
        }
        return seen;
    }
};

struct LowFeedbackForest {
    int n = 0;
    std::vector<char> is_cut;
    std::vector<int> cut_vertices;
    std::vector<int> cut_index;
    std::vector<ull> cut_adj_mask;
    std::vector<int> remaining_vertices;
    std::vector<int> remaining_index;
    std::vector<std::vector<int>> children;
    std::vector<int> roots;
    std::vector<int> order;
};

struct HardComponent {
    Component comp;
    int id = 0;
    long long measure = 0;
};

ll chosen_weight(const std::vector<ll>& weight, const std::vector<char>& chosen) {
    ll total = 0;
    for (int i = 0; i < static_cast<int>(chosen.size()); ++i) {
        if (chosen[i]) {
            total += weight[i];
        }
    }
    return total;
}

bool validate_component_solution(const Component& comp, const std::vector<char>& chosen) {
    if (static_cast<int>(chosen.size()) != static_cast<int>(comp.weight.size())) {
        return false;
    }
    for (int u = 0; u < static_cast<int>(comp.weight.size()); ++u) {
        if (!chosen[u]) {
            continue;
        }
        for (int v : comp.adj[u]) {
            if (u < v && chosen[v]) {
                return false;
            }
        }
    }
    return true;
}

SolveResult solve_tree_component(const Component& comp) {
    int n = static_cast<int>(comp.weight.size());
    std::vector<int> parent(n, -1);
    std::vector<int> order;
    order.reserve(n);

    std::queue<int> q;
    q.push(0);
    parent[0] = 0;

    while (!q.empty()) {
        int v = q.front();
        q.pop();
        order.push_back(v);
        for (int to : comp.adj[v]) {
            if (parent[to] == -1) {
                parent[to] = v;
                q.push(to);
            }
        }
    }

    std::vector<ll> take(n, 0);
    std::vector<ll> skip(n, 0);
    for (int idx = n - 1; idx >= 0; --idx) {
        int v = order[idx];
        take[v] = comp.weight[v];
        skip[v] = 0;
        for (int to : comp.adj[v]) {
            if (parent[to] == v) {
                take[v] += skip[to];
                skip[v] += std::max(take[to], skip[to]);
            }
        }
    }

    std::vector<char> chosen(n, 0);
    std::vector<std::pair<int, bool>> st;
    st.push_back({0, false});
    while (!st.empty()) {
        std::pair<int, bool> state = st.back();
        st.pop_back();
        int v = state.first;
        bool parent_taken = state.second;
        bool take_v = false;
        if (!parent_taken && take[v] >= skip[v]) {
            take_v = true;
            chosen[v] = 1;
        }
        for (int to : comp.adj[v]) {
            if (parent[to] == v) {
                st.push_back({to, take_v});
            }
        }
    }

    return {chosen_weight(comp.weight, chosen), chosen};
}

std::pair<bool, std::vector<int>> bipartite_coloring(const Component& comp) {
    int n = static_cast<int>(comp.weight.size());
    std::vector<int> color(n, -1);
    std::queue<int> q;

    for (int start = 0; start < n; ++start) {
        if (color[start] != -1) {
            continue;
        }
        color[start] = 0;
        q.push(start);
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (int to : comp.adj[v]) {
                if (color[to] == -1) {
                    color[to] = color[v] ^ 1;
                    q.push(to);
                } else if (color[to] == color[v]) {
                    return {false, {}};
                }
            }
        }
    }

    return {true, color};
}

SolveResult solve_bipartite_component(const Component& comp, const std::vector<int>& color) {
    int n = static_cast<int>(comp.weight.size());
    int source = n;
    int sink = n + 1;
    Dinic dinic(n + 2);
    ll total_weight = 0;
    const ll inf = std::numeric_limits<ll>::max() / 4;

    for (int v = 0; v < n; ++v) {
        total_weight += comp.weight[v];
        if (color[v] == 0) {
            dinic.add_edge(source, v, comp.weight[v]);
        } else {
            dinic.add_edge(v, sink, comp.weight[v]);
        }
    }

    for (int u = 0; u < n; ++u) {
        for (int v : comp.adj[u]) {
            if (u < v) {
                if (color[u] == 0) {
                    dinic.add_edge(u, v, inf);
                } else {
                    dinic.add_edge(v, u, inf);
                }
            }
        }
    }

    ll cover_weight = dinic.max_flow(source, sink);
    std::vector<char> reach = dinic.reachable_from_source(source);

    std::vector<char> chosen(n, 0);
    for (int v = 0; v < n; ++v) {
        if (color[v] == 0) {
            if (reach[v]) {
                chosen[v] = 1;
            }
        } else if (!reach[v]) {
            chosen[v] = 1;
        }
    }

    return {total_weight - cover_weight, chosen};
}

bool build_low_feedback_forest(const Component& comp,
                               int max_cut_size,
                               long long max_work,
                               LowFeedbackForest& out) {
    int n = static_cast<int>(comp.weight.size());
    std::vector<int> parent(n, -1);
    std::vector<int> depth(n, 0);
    std::vector<int> visit_order;
    visit_order.reserve(n);

    std::queue<int> q;
    q.push(0);
    parent[0] = 0;
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        visit_order.push_back(v);
        for (int to : comp.adj[v]) {
            if (parent[to] == -1) {
                parent[to] = v;
                depth[to] = depth[v] + 1;
                q.push(to);
            }
        }
    }

    std::vector<std::pair<int, int>> extra_edges;
    for (int u = 0; u < n; ++u) {
        for (int v : comp.adj[u]) {
            if (u < v && parent[u] != v && parent[v] != u) {
                extra_edges.push_back({u, v});
            }
        }
    }

    if (extra_edges.empty()) {
        return false;
    }

    std::vector<char> covered(extra_edges.size(), 0);
    std::vector<char> in_cut(n, 0);
    int uncovered = static_cast<int>(extra_edges.size());

    while (uncovered > 0) {
        std::vector<int> freq(n, 0);
        for (int i = 0; i < static_cast<int>(extra_edges.size()); ++i) {
            if (!covered[i]) {
                ++freq[extra_edges[i].first];
                ++freq[extra_edges[i].second];
            }
        }

        int best_v = -1;
        for (int v = 0; v < n; ++v) {
            if (in_cut[v]) {
                continue;
            }
            if (best_v == -1 || freq[v] > freq[best_v]) {
                best_v = v;
            }
        }

        if (best_v == -1 || freq[best_v] == 0) {
            break;
        }

        in_cut[best_v] = 1;
        for (int i = 0; i < static_cast<int>(extra_edges.size()); ++i) {
            if (!covered[i] &&
                (extra_edges[i].first == best_v || extra_edges[i].second == best_v)) {
                covered[i] = 1;
                --uncovered;
            }
        }
    }

    std::vector<int> cut_vertices;
    for (int v = 0; v < n; ++v) {
        if (in_cut[v]) {
            cut_vertices.push_back(v);
        }
    }

    if (uncovered > 0) {
        return false;
    }
    if (cut_vertices.empty() || static_cast<int>(cut_vertices.size()) > max_cut_size) {
        return false;
    }

    long long work_estimate = static_cast<long long>(n) * (1LL << cut_vertices.size());
    if (work_estimate > max_work) {
        return false;
    }

    out = LowFeedbackForest{};
    out.n = n;
    out.is_cut = std::move(in_cut);
    out.cut_vertices = cut_vertices;
    out.cut_index.assign(n, -1);
    out.cut_adj_mask.assign(n, 0ULL);

    for (int i = 0; i < static_cast<int>(out.cut_vertices.size()); ++i) {
        out.cut_index[out.cut_vertices[i]] = i;
    }

    for (int v = 0; v < n; ++v) {
        if (out.is_cut[v]) {
            continue;
        }
        ull mask = 0ULL;
        for (int to : comp.adj[v]) {
            int id = out.cut_index[to];
            if (id != -1) {
                mask |= (1ULL << id);
            }
        }
        out.cut_adj_mask[v] = mask;
    }

    out.remaining_index.assign(n, -1);
    for (int v = 0; v < n; ++v) {
        if (!out.is_cut[v]) {
            out.remaining_index[v] = static_cast<int>(out.remaining_vertices.size());
            out.remaining_vertices.push_back(v);
        }
    }

    int rem = static_cast<int>(out.remaining_vertices.size());
    out.children.assign(rem, {});
    std::vector<int> parent_rem(rem, -1);
    std::vector<char> seen(rem, 0);

    for (int root_pos = 0; root_pos < rem; ++root_pos) {
        if (seen[root_pos]) {
            continue;
        }
        seen[root_pos] = 1;
        out.roots.push_back(root_pos);
        std::queue<int> bfs;
        bfs.push(root_pos);

        while (!bfs.empty()) {
            int pos = bfs.front();
            bfs.pop();
            int v = out.remaining_vertices[pos];
            out.order.push_back(pos);
            for (int to : comp.adj[v]) {
                int to_pos = out.remaining_index[to];
                if (to_pos == -1 || seen[to_pos]) {
                    continue;
                }
                seen[to_pos] = 1;
                parent_rem[to_pos] = pos;
                out.children[pos].push_back(to_pos);
                bfs.push(to_pos);
            }
        }
    }

    return true;
}

SolveResult solve_low_feedback_component(const Component& comp, const LowFeedbackForest& forest) {
    int n = forest.n;
    int rem = static_cast<int>(forest.remaining_vertices.size());
    int cut_n = static_cast<int>(forest.cut_vertices.size());

    std::vector<ull> cut_graph(cut_n, 0ULL);
    for (int i = 0; i < cut_n; ++i) {
        int u = forest.cut_vertices[i];
        for (int v : comp.adj[u]) {
            int j = forest.cut_index[v];
            if (j != -1) {
                cut_graph[i] |= (1ULL << j);
            }
        }
    }

    std::vector<ll> take(rem, 0);
    std::vector<ll> skip(rem, 0);
    ull best_mask = 0;
    ll best_weight = -1;

    std::function<void(int, ull, ll)> dfs = [&](int idx, ull selected_mask, ll cut_weight) {
        if (idx == cut_n) {
            for (int pos_idx = rem - 1; pos_idx >= 0; --pos_idx) {
                int pos = forest.order[pos_idx];
                int v = forest.remaining_vertices[pos];
                bool blocked = (forest.cut_adj_mask[v] & selected_mask) != 0ULL;
                ll cur_take = blocked ? NEG_INF : comp.weight[v];
                ll cur_skip = 0;
                for (int child : forest.children[pos]) {
                    cur_skip += std::max(take[child], skip[child]);
                    if (!blocked) {
                        cur_take += skip[child];
                    }
                }
                take[pos] = cur_take;
                skip[pos] = cur_skip;
            }

            ll total = cut_weight;
            for (int root : forest.roots) {
                total += std::max(take[root], skip[root]);
            }

            if (total > best_weight) {
                best_weight = total;
                best_mask = selected_mask;
            }
            return;
        }

        dfs(idx + 1, selected_mask, cut_weight);
        if ((selected_mask & cut_graph[idx]) == 0ULL) {
            dfs(idx + 1, selected_mask | (1ULL << idx), cut_weight + comp.weight[forest.cut_vertices[idx]]);
        }
    };

    dfs(0, 0ULL, 0LL);

    std::vector<char> chosen(n, 0);
    for (int i = 0; i < cut_n; ++i) {
        if ((best_mask >> i) & 1ULL) {
            chosen[forest.cut_vertices[i]] = 1;
        }
    }

    for (int pos_idx = rem - 1; pos_idx >= 0; --pos_idx) {
        int pos = forest.order[pos_idx];
        int v = forest.remaining_vertices[pos];
        bool blocked = (forest.cut_adj_mask[v] & best_mask) != 0ULL;
        ll cur_take = blocked ? NEG_INF : comp.weight[v];
        ll cur_skip = 0;
        for (int child : forest.children[pos]) {
            cur_skip += std::max(take[child], skip[child]);
            if (!blocked) {
                cur_take += skip[child];
            }
        }
        take[pos] = cur_take;
        skip[pos] = cur_skip;
    }

    std::vector<std::pair<int, bool>> stack;
    for (int root : forest.roots) {
        stack.push_back({root, false});
    }

    while (!stack.empty()) {
        std::pair<int, bool> state = stack.back();
        stack.pop_back();
        int pos = state.first;
        bool parent_taken = state.second;
        int v = forest.remaining_vertices[pos];
        bool blocked = (forest.cut_adj_mask[v] & best_mask) != 0ULL;
        bool take_v = false;
        if (!parent_taken && !blocked && take[pos] >= skip[pos]) {
            take_v = true;
            chosen[v] = 1;
        }
        for (int child : forest.children[pos]) {
            stack.push_back({child, take_v});
        }
    }

    return {best_weight, chosen};
}

void fail_with_message(const std::string& message) {
    std::cerr << message << '\n';
    std::exit(1);
}

void apply_component_solution(const Component& comp,
                              const SolveResult& result,
                              std::vector<char>& global_chosen) {
    for (int i = 0; i < static_cast<int>(comp.weight.size()); ++i) {
        if (result.chosen[i]) {
            global_chosen[comp.global[i]] = 1;
        }
    }
}

bool solve_exact_component(const Component& comp, int component_id, SolveResult& result) {
    constexpr int SMALL_EXACT_LIMIT = 90;
    constexpr int MAX_LOW_FEEDBACK_CUT = 22;
    constexpr long long MAX_LOW_FEEDBACK_WORK = 80000000LL;

    int comp_n = static_cast<int>(comp.weight.size());
    bool solved_exact = false;

    if (comp.edges == 0) {
        result.chosen.assign(comp_n, 1);
        result.weight = chosen_weight(comp.weight, result.chosen);
        solved_exact = true;
    } else if (comp.edges == comp_n - 1) {
        result = solve_tree_component(comp);
        solved_exact = true;
    } else if (comp_n <= SMALL_EXACT_LIMIT) {
        SmallExactMWIS exact(comp.weight, comp.adj);
        result = exact.solve();
        solved_exact = true;
    } else {
        std::pair<bool, std::vector<int>> bip = bipartite_coloring(comp);
        bool is_bipartite = bip.first;
        const std::vector<int>& color = bip.second;
        if (is_bipartite) {
            result = solve_bipartite_component(comp, color);
            solved_exact = true;
        } else {
            LowFeedbackForest forest;
            if (build_low_feedback_forest(comp, MAX_LOW_FEEDBACK_CUT, MAX_LOW_FEEDBACK_WORK, forest)) {
                result = solve_low_feedback_component(comp, forest);
                solved_exact = true;
            }
        }
    }

    if (solved_exact && !validate_component_solution(comp, result.chosen)) {
        fail_with_message("Invalid exact solution on component " + std::to_string(component_id));
    }

    return solved_exact;
}

std::vector<HardComponent> collect_hard_components(const ProblemInstance& instance,
                                                   const GraphData& graph_data,
                                                   std::vector<char>& global_chosen) {
    std::vector<HardComponent> hard_components;
    hard_components.reserve(graph_data.component_vertices.size());

    std::vector<int> reusable_local_index(instance.n, -1);

    for (int cid = 0; cid < static_cast<int>(graph_data.component_vertices.size()); ++cid) {
        Component comp = build_component(
            graph_data.component_vertices[cid], instance.weight, graph_data.graph, reusable_local_index);

        SolveResult result;
        if (solve_exact_component(comp, cid, result)) {
            apply_component_solution(comp, result, global_chosen);
        } else {
            int comp_n = static_cast<int>(comp.weight.size());
            hard_components.push_back(
                HardComponent{std::move(comp), cid, static_cast<long long>(comp_n) + 2LL * comp.edges});
        }
    }

    return hard_components;
}

void solve_hard_components(std::vector<HardComponent>& hard_components,
                           const Timer& timer,
                           std::vector<char>& global_chosen) {
    for (int i = 0; i < static_cast<int>(hard_components.size()); ++i) {
        if (timer.expired(150)) {
            break;
        }

        long long remaining_measure = 0;
        for (int j = i; j < static_cast<int>(hard_components.size()); ++j) {
            remaining_measure += hard_components[j].measure;
        }

        long long remaining_time = timer.remaining_ms();
        long long budget_ms = remaining_measure == 0
                                  ? 0
                                  : std::max(200LL,
                                             (remaining_time - 120LL) * hard_components[i].measure /
                                                 std::max(1LL, remaining_measure));

        HeuristicSolver heuristic(hard_components[i].comp, 0xC0D3D00DULL + static_cast<uint64_t>(hard_components[i].id));
        SolveResult result = heuristic.solve(timer, budget_ms);
        if (!validate_component_solution(hard_components[i].comp, result.chosen)) {
            fail_with_message("Invalid heuristic solution on component " +
                              std::to_string(hard_components[i].id));
        }
        apply_component_solution(hard_components[i].comp, result, global_chosen);
    }
}

FinalSolution build_final_solution(const ProblemInstance& instance, const std::vector<char>& global_chosen) {
    ll total_weight = 0;
    std::vector<int> picked;
    picked.reserve(instance.n);

    for (int v = 0; v < instance.n; ++v) {
        if (global_chosen[v]) {
            total_weight += instance.weight[v];
            picked.push_back(v + 1);
        }
    }

    for (const auto& edge : instance.edges) {
        int u = edge.first;
        int v = edge.second;
        if (global_chosen[u] && global_chosen[v]) {
            fail_with_message("Final solution is invalid");
        }
    }

    return FinalSolution{total_weight, picked};
}

}  // namespace

FinalSolution solve_problem(const ProblemInstance& instance, long long time_limit_ms) {
    GraphData graph_data = build_graph_data(instance.n, instance.edges);
    Timer timer(time_limit_ms);

    std::vector<char> global_chosen(instance.n, 0);
    std::vector<HardComponent> hard_components = collect_hard_components(instance, graph_data, global_chosen);
    solve_hard_components(hard_components, timer, global_chosen);
    return build_final_solution(instance, global_chosen);
}

}  // namespace hs

using std::cin;
using std::cout;

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    cin.tie(nullptr);

    hs::ProblemInstance instance;
    cin >> instance.n >> instance.m;

    instance.weight.resize(instance.n);
    for (int i = 0; i < instance.n; ++i) {
        cin >> instance.weight[i];
    }

    instance.edges.reserve(instance.m);
    for (int i = 0; i < instance.m; ++i) {
        int u = 0;
        int v = 0;
        cin >> u >> v;
        instance.edges.push_back({u - 1, v - 1});
    }

    long long time_limit_ms = hs::read_time_limit_ms(argc, argv);
    hs::FinalSolution solution = hs::solve_problem(instance, time_limit_ms);

    cout << solution.total_weight << '\n';
    for (int i = 0; i < static_cast<int>(solution.picked.size()); ++i) {
        if (i > 0) {
            cout << ' ';
        }
        cout << solution.picked[i];
    }
    cout << '\n';

    return 0;
}
