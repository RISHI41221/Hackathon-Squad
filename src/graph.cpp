#include "graph.h"

#include <queue>

namespace hs {

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

}  // namespace hs
