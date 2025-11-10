#pragma once
#include <atomic>
#include <chrono>
#include "util/logger.hpp"

namespace phase {

inline std::atomic<long long>& slot_io()    { static std::atomic<long long> v{0}; return v; }
inline std::atomic<long long>& slot_build() { static std::atomic<long long> v{0}; return v; }
inline std::atomic<long long>& slot_logic() { static std::atomic<long long> v{0}; return v; }

struct Accum {
    std::atomic<long long>& slot;
    std::chrono::steady_clock::time_point t0;
    explicit Accum(std::atomic<long long>& s)
        : slot(s), t0(std::chrono::steady_clock::now()) {}
    ~Accum() {
        const auto t1 = std::chrono::steady_clock::now();
        const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
        slot.fetch_add(ns, std::memory_order_relaxed);
    }
};

inline void reset() {
    slot_io   ().store(0, std::memory_order_relaxed);
    slot_build().store(0, std::memory_order_relaxed);
    slot_logic().store(0, std::memory_order_relaxed);
}

inline void report() {
    const double io_s    = slot_io   ().load(std::memory_order_relaxed) / 1e9;
    const double build_s = slot_build().load(std::memory_order_relaxed) / 1e9;
    const double logic_s = slot_logic().load(std::memory_order_relaxed) / 1e9;
    logger::info("Phase times (global): I/O {:.3f}s | BUILD (BC+SPQR) {:.3f}s | LOGIC {:.3f}s",
                 io_s, build_s, logic_s);
}

} 

#define ACCUM_IO()    phase::Accum __acc_io_##__LINE__    ( phase::slot_io() )
#define ACCUM_BUILD() phase::Accum __acc_build_##__LINE__ ( phase::slot_build() )
#define ACCUM_LOGIC() phase::Accum __acc_logic_##__LINE__ ( phase::slot_logic() )

#define ACCUM_BC()    ACCUM_BUILD()
#define ACCUM_SPQR()  ACCUM_BUILD()
