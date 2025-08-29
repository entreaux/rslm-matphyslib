#pragma once
/**
 * RSLM Maths — rng.hpp
 * --------------------
 * PCG32-based RNG with Gaussian (Box–Muller) and U[0,1) helpers.
 * Deterministic, cheap. Per-sample logging is THROTTLED by default.
 * Enable heavy tracing with: export RSLM_TRACE_HEAVY=1
 */

#include "config.hpp"
#include "trace.hpp"
#include <cstdint>
#include <cmath>
#include <cstdlib>

namespace rslm::rng {

using rslm::cfg::real;

// --- Heavy trace gate (env controlled) --------------------------------------
inline bool heavy_trace_enabled() {
    static int v = []{
        if (const char* e = std::getenv("RSLM_TRACE_HEAVY")) {
            return (std::string(e) != "0") ? 1 : 0;
        }
        return 0; // default OFF
    }();
    return v != 0;
}

// Throttled emitter: logs first K, then every 'stride'
struct Throttle {
    std::uint64_t n = 0;
    std::uint64_t first = 16;
    std::uint64_t stride = 10000;

    bool tick() {
        ++n;
        if (n <= first) return true;
        return (stride > 0) && (n % stride == 0);
    }
};

// Minimal PCG32 (state+inc)
struct PCG32 {
    std::uint64_t state = 0x853c49e6748fea9bULL;
    std::uint64_t inc   = 0xda3e39cb94b95bdbULL; // must be odd

    // per-instance throttles so multiple RNGs don’t share counters
    Throttle th_u;
    Throttle th_n;

    explicit PCG32(std::uint64_t seed = 0u, std::uint64_t seq = 1u) { reseed(seed, seq); }

    void reseed(std::uint64_t seed, std::uint64_t seq=1u) {
        state = 0u;
        inc   = (seq << 1u) | 1u;
        next_u32();
        state += seed;
        next_u32();
        TRACE_INFO("pcg32_seed", "seed=" + std::to_string(seed) + " seq=" + std::to_string(seq));
        th_u = Throttle{}; th_n = Throttle{}; // reset throttles
    }

    inline std::uint32_t next_u32() {
        std::uint64_t oldstate = state;
        state = oldstate * 6364136223846793005ULL + inc;
        std::uint32_t xorshifted = static_cast<std::uint32_t>(((oldstate >> 18u) ^ oldstate) >> 27u);
        std::uint32_t rot = static_cast<std::uint32_t>(oldstate >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

    // U[0,1)
    inline real uniform01() {
        // 24 or 53 bits depending on real
        real u = (next_u32() + real(0.5)) / real(4294967296.0);
        if (heavy_trace_enabled() && th_u.tick()) {
            TRACE_DEBUG("uniform01", u);
        }
        return u;
    }

    // Standard normal via Box–Muller (basic)
    inline real normal01() {
        const real u1 = uniform01();
        const real u2 = uniform01();
        const real r  = std::sqrt(real(-2.0) * std::log(std::max(u1, real(1e-12))));
        const real th = real(6.2831853071795864769) * u2; // 2π
        real z = r * std::cos(th);
        if (heavy_trace_enabled() && th_n.tick()) {
            TRACE_DEBUG("normal01", z);
        }
        return z;
    }
};

} // namespace rslm::rng
