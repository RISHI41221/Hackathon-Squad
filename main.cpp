#include "solver.h"
#include "utils.h"

#include <iostream>
#include <utility>
#include <vector>

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
