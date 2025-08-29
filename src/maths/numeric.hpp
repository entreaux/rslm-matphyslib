#pragma once
/**
 * RSLM Maths — numeric.hpp/.cpp
 * -----------------------------
 * Small, branch-safe numeric helpers with robust edge-case handling.
 * All functions emit telemetry for inputs/outputs (summarized when heavy).
 */

#include "config.hpp"
#include "trace.hpp"
#include <cmath>
#include <limits>
#include <type_traits>

namespace rslm::num {

using rslm::cfg::real;

// Generic clamp
template <typename T>
inline T clamp(T x, T lo, T hi) {
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    TRACE_DEBUG("clamp_in", std::to_string(double(x)) + " lo=" + std::to_string(double(lo)) + " hi=" + std::to_string(double(hi)));
    return x;
}

// Linear interpolation (stable at t=0/1)
inline real lerp(real a, real b, real t) {
    // (1-t)*a + t*b loses precision when a,b large and t small; use FMA-ish split
    real r = a + (b - a) * t;
    TRACE_DEBUG("lerp", "a=" + std::to_string(double(a)) + " b=" + std::to_string(double(b)) + " t=" + std::to_string(double(t)) +
                        " -> " + std::to_string(double(r)));
    return r;
}

// Safe sqrt: clamps tiny negatives from roundoff to zero
inline real safe_sqrt(real x) {
    real y = x;
    if (RSLM_UNLIKELY(y < real(0))) {
        if (y > -rslm::cfg::real(1) * std::numeric_limits<real>::epsilon() * rslm::cfg::real(10)) y = real(0);
    }
    real r = std::sqrt(y);
    TRACE_DEBUG("safe_sqrt", "x=" + std::to_string(double(x)) + " y=" + std::to_string(double(y)) + " r=" + std::to_string(double(r)));
    return r;
}

// Stable hypot for (a,b) to avoid overflow/underflow
inline real stable_hypot(real a, real b) {
    real r = std::hypot(a, b);
    TRACE_DEBUG("stable_hypot", "a=" + std::to_string(double(a)) + " b=" + std::to_string(double(b)) + " r=" + std::to_string(double(r)));
    return r;
}

// Saturating tanh that never returns exactly ±1
inline real saturating_tanh(real x) {
    real t = std::tanh(x);
    // Pull back from ±1 by epsilon to avoid later divisions by (1-t^2)=0
    const real eps = rslm::cfg::real(4) * std::numeric_limits<real>::epsilon();
    if (t >= real(1)) t = real(1) - eps;
    if (t <= real(-1)) t = real(-1) + eps;
    TRACE_DEBUG("saturating_tanh", "x=" + std::to_string(double(x)) + " t=" + std::to_string(double(t)));
    return t;
}

// Almost-equal with mixed abs/rel tolerances
inline bool almost_equal(real a, real b, real atol, real rtol) {
    real diff = std::fabs(a - b);
    real thr  = atol + rtol * std::max(std::fabs(a), std::fabs(b));
    bool ok = diff <= thr;
    TRACE_DEBUG("almost_equal", "a=" + std::to_string(double(a)) + " b=" + std::to_string(double(b)) +
                               " diff=" + std::to_string(double(diff)) + " thr=" + std::to_string(double(thr)) +
                               " ok=" + std::to_string(ok));
    return ok;
}

// Finite checker with logging
inline bool is_finite(real x) {
    bool ok = std::isfinite(x);
    TRACE_DEBUG("is_finite", "x=" + std::to_string(double(x)) + " ok=" + std::to_string(ok));
    return ok;
}

} // namespace rslm::num
