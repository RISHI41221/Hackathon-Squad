#pragma once

#include "utils.h"

#include <utility>
#include <vector>

namespace hs {

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

}  // namespace hs
