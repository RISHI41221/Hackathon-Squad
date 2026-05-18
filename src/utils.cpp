#include "utils.h"

#include <algorithm>
#include <cstdlib>

namespace hs {

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

}  // namespace hs
