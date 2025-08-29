#pragma once
/**
 * RSLM Maths — diagnostics/grid.hpp
 * ---------------------------------
 * Sample scalar diagnostics (e.g., curvature) over 2D spacetime slices.
 * - Grid2D stores an axis-aligned regular lattice and values.
 * - sample_xy() samples a scalar function f(t,x,y,z) on (x,y) at fixed (t0,z0).
 * - curv_scalar() and curv_riemann_frob() compute curvature scalars at a point.
 *
 * No per-sample logging; only summary stats are emitted by the smoke test.
 */

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>

#include "config.hpp"
#include "linalg.hpp"
#include "curvature.hpp"
#include "connection.hpp"
#include "field.hpp"

namespace rslm::diag {

using rslm::cfg::real;
using rslm::linalg::vec4;
using rslm::linalg::mat4;
using rslm::field::IMetricField;

struct Grid2D {
    // origin at (x0,y0), spacing (dx,dy), size nx × ny (row-major: i=row/y, j=col/x)
    std::size_t nx{0}, ny{0};
    real x0{0}, y0{0}, dx{1}, dy{1};
    real t0{0}, z0{0};          // the fixed coordinates of the slice
    std::vector<real> val;      // size = nx*ny

    inline real& at(std::size_t i, std::size_t j)       { return val[i*ny + j]; }
    inline const real& at(std::size_t i, std::size_t j) const { return val[i*ny + j]; }
};

// ---- Curvature scalars ------------------------------------------------------

// Scalar curvature R(x) using field F
inline real curv_scalar(const IMetricField& F, const vec4& x) {
    auto P = rslm::conn::prepare_metric(F, x);
    auto Rm = rslm::curv::riemann_at(F, x);
    auto Rc = rslm::curv::ricci(Rm);
    return rslm::curv::scalar(P.g_inv, Rc);
}

// Frobenius norm ||R||_F at x
inline real curv_riemann_frob(const IMetricField& F, const vec4& x) {
    auto Rm = rslm::curv::riemann_at(F, x);
    return rslm::curv::frob_riemann(Rm);
}

// ---- Generic XY sampler -----------------------------------------------------

/**
 * Sample a scalar function s(x) on an XY slice at fixed (t0,z0).
 * f must be: real f(const IMetricField&, const vec4&).
 */
template <typename ScalarFn>
inline Grid2D sample_xy(const IMetricField& F, real t0, real z0,
                        real x0, real y0, real dx, real dy,
                        std::size_t nx, std::size_t ny,
                        ScalarFn f)
{
    Grid2D G; G.nx=nx; G.ny=ny; G.x0=x0; G.y0=y0; G.dx=dx; G.dy=dy; G.t0=t0; G.z0=z0;
    G.val.assign(nx*ny, real(0));
    vec4 x(t0, x0, y0, z0);
    for (std::size_t i=0;i<nx;++i) {
        x.v[2] = y0 + real(i)*dy;       // y row
        for (std::size_t j=0;j<ny;++j) {
            x.v[1] = x0 + real(j)*dx;   // x col
            G.at(i,j) = f(F, x);
        }
    }
    return G;
}

// ---- Simple stats -----------------------------------------------------------

struct Stats {
    real vmin{0}, vmax{0}, mean{0};
};

inline Stats stats(const Grid2D& G) {
    Stats s;
    if (G.val.empty()) return s;
    s.vmin = s.vmax = G.val[0];
    long double acc = 0;
    for (real v : G.val) { if (v<s.vmin) s.vmin=v; if (v>s.vmax) s.vmax=v; acc += v; }
    s.mean = static_cast<real>(acc / static_cast<long double>(G.val.size()));
    return s;
}

} // namespace rslm::diag
