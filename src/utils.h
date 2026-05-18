#pragma once

#include <chrono>
#include <cstdint>
#include <limits>
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

}  // namespace hs
