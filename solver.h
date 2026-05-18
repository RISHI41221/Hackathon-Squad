#pragma once

#include "utils.h"

#include <utility>
#include <vector>

namespace hs {

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

}  // namespace hs
