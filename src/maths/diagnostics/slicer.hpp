#pragma once
/**
 * RSLM Maths — diagnostics/slicer.hpp
 * -----------------------------------
 * Generic plane sampling over (axis i, axis j) with the other two fixed.
 * Convenience wrappers: sample_xy, sample_xz, sample_ty.
 */

#include <cstddef>
#include <array>
#include <functional>
#include "config.hpp"
#include "linalg.hpp"
#include "field.hpp"
#include "grid.hpp"

namespace rslm::diag {

using rslm::cfg::real;
using rslm::linalg::vec4;
using rslm::field::IMetricField;

template <typename ScalarFn>
inline Grid2D sample_plane(const IMetricField& F,
                           int axi, int axj,
                           const std::array<real,4>& fixed, // all 4 coords; axi/axj overwritten
                           real u0, real v0, real du, real dv,
                           std::size_t nu, std::size_t nv,
                           ScalarFn f)
{
    Grid2D G;
    G.nx = nu; G.ny = nv;
    G.x0 = v0; G.y0 = u0; G.dx = dv; G.dy = du;
    G.t0 = fixed[0]; G.z0 = fixed[3];
    G.val.assign(nu*nv, real(0));

    vec4 x(fixed[0], fixed[1], fixed[2], fixed[3]);

    for (std::size_t i=0;i<nu;++i) {
        for (std::size_t j=0;j<nv;++j) {
            x.v[axi] = u0 + real(i)*du;
            x.v[axj] = v0 + real(j)*dv;
            G.at(i,j) = f(F, x);
        }
    }
    return G;
}

/** XZ slice at fixed (t0, y0). */
template <typename ScalarFn>
inline Grid2D sample_xz(const IMetricField& F,
                        real t0, real y0,
                        real x0, real z0, real dx, real dz,
                        std::size_t nx, std::size_t nz,
                        ScalarFn f)
{
    return sample_plane(F, /*axi=z*/3, /*axj=x*/1,
                        std::array<real,4>{t0, x0, y0, z0},
                        /*u0=*/z0, /*v0=*/x0, /*du=*/dz, /*dv=*/dx,
                        /*nu=*/nz, /*nv=*/nx, f);
}

/** TY slice at fixed (x0, z0). */
template <typename ScalarFn>
inline Grid2D sample_ty(const IMetricField& F,
                        real x0, real z0,
                        real t0, real y0, real dt, real dy,
                        std::size_t nt, std::size_t ny,
                        ScalarFn f)
{
    return sample_plane(F, /*axi=y*/2, /*axj=t*/0,
                        std::array<real,4>{t0, x0, y0, z0},
                        /*u0=*/y0, /*v0=*/t0, /*du=*/dy, /*dv=*/dt,
                        /*nu=*/ny, /*nv=*/nt, f);
}


} // namespace rslm::diag
