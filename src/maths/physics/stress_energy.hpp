#pragma once
/**
 * RSLM Maths — physics/stress_energy.hpp
 * --------------------------------------
 * Minimal semantic stress–energy analogue for diagnostics:
 *
 *   T_{μν}(x) = Σ_i φσ(d̃(x, x_i)) [ E_i u_{iμ} u_{iν} + η m_i c^2 g_{μν}(x) ]
 *
 * where:
 *  - φσ(r) = exp( - r^2 / (2 σ^2) )
 *  - d̃ is a *PD* local distance defined with \tilde g(x) = g(x)^2:
 *      d̃^2 = (x - x_i)^T \tilde g(x) (x - x_i)
 *    This avoids Lorentzian signature issues while still being geometry-aware.
 *
 * Notes:
 *  - u_{iμ} = g_{μα}(x) u_i^α (we “lower” with *local* g at the query point).
 *  - Units: c=cfg::c (hyperparameter), η small stabilizer.
 */

#include <vector>
#include <cmath>
#include <algorithm>

#include "units.hpp"
#include "config.hpp"
#include "linalg.hpp"
#include "metric.hpp"
#include "quadform.hpp"
#include "field.hpp"
#include "pd_proxy.hpp"

namespace rslm::phys {

using rslm::cfg::real;
using rslm::linalg::vec4;
using rslm::linalg::mat4;
using rslm::linalg::sym4;
using rslm::field::IMetricField;
using rslm::units::real;

struct Event {
    vec4 x;     // position (t,x,y,z)
    vec4 u;     // 4-velocity (approximately timelike, not required normalized)
    real E{1};  // semantic energy ≥ 0
    real m{1};  // “rest mass” > 0
};

struct TSParams {
    real sigma{1.0};   // kernel width
    real eta{1e-2};    // stabilizer for g-term
    real kappa{0.1};   // coupling constant
    real c2{1.0};
};

inline mat4 lower_u_outer(const sym4& g, const vec4& u) {
    // u_mu = g_{μa} u^a ; return u_mu u_nu outer product
    real ulim[4]{};
    for (int mu=0; mu<4; ++mu) {
        real s = 0;
        for (int a=0; a<4; ++a) s += g.m[mu][a] * u.v[a];
        ulim[mu] = s;
    }
    mat4 out{};
    for (int mu=0; mu<4; ++mu)
        for (int nu=0; nu<4; ++nu)
            out.m[mu][nu] = ulim[mu] * ulim[nu];
    return out;
}

inline mat4 scaled(const mat4& A, real s) {
    mat4 B = A;
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) B.m[i][j] *= s;
    return B;
}

inline mat4 add(const mat4& A, const mat4& B) {
    mat4 C{};
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) C.m[i][j] = A.m[i][j] + B.m[i][j];
    return C;
}

// φσ(r) with r^2 (already squared); expects σ>0
inline real kernel_exp(real r2, real sigma) {
    const real inv = real(0.5) / (sigma*sigma);
    return std::exp(- r2 * inv);
}

/** Compute T_{μν}(x) from events at query x using local metric g(x). */
inline mat4 stress_energy_at(const IMetricField& F,
                             const std::vector<Event>& evs,
                             const vec4& x,
                             const TSParams& P)
{
    sym4 g = F.g(x);
    // PD proxy for local distance
    mat4 gtilde = rslm::audits::pd_proxy_square(g);

    auto d2 = [&](const vec4& xi) {
        // (x - xi)^T gtilde (x - xi) ; gtilde is PD
        vec4 d{ x.v[0]-xi.v[0], x.v[1]-xi.v[1], x.v[2]-xi.v[2], x.v[3]-xi.v[3] };
        // y = gtilde * d
        real y[4]{};
        for (int i=0;i<4;++i) {
            real s=0; for (int j=0;j<4;++j) s += gtilde.m[i][j]*d.v[j];
            y[i]=s;
        }
        real s=0; for (int i=0;i<4;++i) s += d.v[i]*y[i];
        return s;
    };

    mat4 T{}; // accumulate
    for (const auto& e : evs) {
        real w = kernel_exp(d2(e.x), P.sigma);
        // lower with *local* g(x)
        mat4 umu_unu = lower_u_outer(g, e.u);
        mat4 term1   = scaled(umu_unu, e.E);
        // eta m c^2 g_{μν}  (promote g -> mat4)
        mat4 gmat{}; for (int i=0;i<4;++i) for (int j=0;j<4;++j) gmat.m[i][j]=g.m[i][j];
        mat4 term2 = scaled(gmat, P.eta * e.m * P.c2);
        T = add(T, scaled(add(term1, term2), w));
    }
    return T;
}

} // namespace rslm::physs
