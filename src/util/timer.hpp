#pragma once
#include <chrono>
#include <cstdio>
#include "util/context.hpp"

class ScopeTimer
{
public:
    explicit ScopeTimer(const char* label) noexcept
        : label_{label},
          start_{std::chrono::steady_clock::now()}
    {}

    ~ScopeTimer() noexcept
    {
        auto &C = ctx();

        if (!C.timingEnabled) return;

        using namespace std::chrono;
        const auto ms =
            duration_cast<milliseconds>(steady_clock::now() - start_).count();

        std::fprintf(stderr, "[timer] %s : %ld ms\n", label_, ms);
    }

    ScopeTimer(const ScopeTimer&)            = delete;
    ScopeTimer& operator=(const ScopeTimer&) = delete;

private:
    const char*                                         label_;
    const std::chrono::steady_clock::time_point         start_;
};

#define TIME_BLOCK(label)  ::ScopeTimer timer_##__COUNTER__{label}

#define TIME_BLOCK_SCOPE(label, code)                     \
    {                                                     \
        TIME_BLOCK(label);                                \
        code                                              \
    }
