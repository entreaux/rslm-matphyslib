#pragma once
/**
 * RSLM Maths — units.hpp
 * ----------------------
 * Physical constants, numeric epsilons, and step sizes used across geometry.
 * We keep c explicit as a *hyperparameter*, defaulting to 1 in natural units.
 */

#include "config.hpp"
#include "trace.hpp"
#include <cmath>
#include <cstdlib>

namespace rslm::units {

using rslm::cfg::real;

// --------- Physical constants (runtime-modifiable where it helps) ------------
struct Constants {
    real c = static_cast<real>(1.0);     // speed of light (hyperparam, default 1)
    real g_eps = static_cast<real>(1e-12); // generic epsilon for metric ops
    real fd_h  = static_cast<real>(1e-4);  // finite-diff step for ∂g
    real dtau  = static_cast<real>(0.5);   // integrator step (blueprint default)
};

// Global constants instance (header-only, ODR-safe since inline)
inline Constants& C() {
    static Constants k;
    return k;
}

// --------- Numeric epsilons tuned to chosen 'real' ---------------------------
inline constexpr real EPS()       { return rslm::cfg::kRealName[0]=='f' ? real(1e-6)  : real(1e-12); }
inline constexpr real SQRT_EPS()  { return rslm::cfg::kRealName[0]=='f' ? real(1e-3)  : real(1e-6);  }
inline constexpr real TINY()      { return rslm::cfg::kRealName[0]=='f' ? real(1e-12) : real(1e-18); }

// Snapshot current constants to logs (for reproducibility)
inline void log_units_snapshot() {
    TRACE_SCOPE("units_snapshot");
    TRACE_INFO("real_type", rslm::cfg::kRealName);
    TRACE_INFO("c", C().c);
    TRACE_INFO("g_eps", C().g_eps);
    TRACE_INFO("fd_h", C().fd_h);
    TRACE_INFO("dtau", C().dtau);
}

} // namespace rslm::units
