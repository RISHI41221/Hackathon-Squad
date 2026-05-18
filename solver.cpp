#include "solver.h"

#include "graph.h"
#include "heuristics.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace hs {
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
