#pragma once
/**
 * RSLM Maths — config.hpp
 * -----------------------
 * Central compile-time knobs and tiny runtime flags used across maths.
 * Keep this header lean and dependency-free.
 */

#include <cstdint>
#include <cstddef>

namespace rslm::cfg {

/* ------------------- CONFIGURATION ------------------*/
//-----------------------------------------------------//

/* ------------------- Physics constants --------------*/

constexpr float                   c = 1.0;    // Lightspeed constant
constexpr float                  c2 = c*c;    // Squared c constant
constexpr float               kappa = 1e-1;   // Kappa Einstein-fit coupling G_μν​≈κT_μν
constexpr float                 eta = 1e-2;   // Stabiliser in G_μν
constexpr float               sigma = 1e-0;   // Kernel width

/* ---------------- Geomeetry constants ---------------*/

constexpr float       gaussian_bump = 5e-2;   // Amplitude for the metric bump
constexpr float      gaussian_sigma = 5e-3;   // Spatial width of the bump

/* ------------ Differentiation constants -------------*/

constexpr float                fd_h = 1e-4;   // Finite-difference step
constexpr float              g_eps = 1e-12;   // G epsilon
constexpr float               dtau = 1e-2;    // Default proper-time


//---------------------------- Precision ---------------------------------------
// Default to double for geometry robustness; flip to float for speed/memory.
#if defined(RSLM_REAL_FLOAT)
using real = float;
#else
using real = double;
#endif

// Compile-time banner
inline constexpr const char* kRealName =
#if defined(RSLM_REAL_FLOAT)
    "float";
#else
    "double";
#endif

//---------------------------- Sanity ------------------------------------------
#if defined(__clang__)
#  define RSLM_ASSUME(x) __builtin_assume(x)
#else
#  define RSLM_ASSUME(x) ((void)0)
#endif

#if defined(__clang__)
#  define RSLM_LIKELY(x)   __builtin_expect(!!(x), 1)
#  define RSLM_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#  define RSLM_LIKELY(x)   (x)
#  define RSLM_UNLIKELY(x) (x)
#endif

// Turn on extra expensive argument logging in hot helpers (off by default).
// You can export RSLM_TRACE_HEAVY=1 at runtime; see numeric/rng for usage.
inline bool trace_heavy_default() noexcept { return false; }

} // namespace rslm::cfg
