#pragma once
/**
 * RSLM Maths — deriv.hpp
 * ----------------------
 * Central finite differences for metric fields and potentials.
 * We keep it explicit (no templates over tensor types) for clarity.
 */

#include "config.hpp"
#include "linalg.hpp"
#include "field.hpp"
#include "units.hpp"
#include "trace.hpp"

namespace rslm::deriv {

using rslm::cfg::real;
using rslm::linalg::vec4;
using rslm::linalg::mat4;
using rslm::linalg::sym4;
using rslm::field::IMetricField;
using rslm::field::IPotential;

// ∂_a g_{μν}(x) for a ∈ {0..3}. Returns a 4×4 matrix of partials (still symmetric).
inline mat4 dmetric(const IMetricField& F, const vec4& x, int a, real h = rslm::units::C().fd_h) {
    vec4 xp = x, xm = x;
    xp.v[a] += h; xm.v[a] -= h;
    sym4 gp = F.g(xp);
    sym4 gm = F.g(xm);
    mat4 d;
    const real s = real(0.5) / h;
    for (int r=0;r<4;++r) {
        for (int c=0;c<4;++c) {
            d.m[r][c] = (gp.m[r][c] - gm.m[r][c]) * s;
        }
    }
    return d;
}

// Full set of partials: dg[a] = ∂_a g
struct DMetric4 {
    mat4 dg[4];
};

inline DMetric4 dmetric4(const IMetricField& F, const vec4& x, real h = rslm::units::C().fd_h) {
    DMetric4 out;
    for (int a=0;a<4;++a) out.dg[a] = dmetric(F, x, a, h);
    return out;
}

// Potential gradient: (∂_0 V, ∂_1 V, ∂_2 V, ∂_3 V)
inline vec4 gradV(const IPotential& P, const vec4& x, real h = rslm::units::C().fd_h) {
    vec4 g;
    for (int a=0;a<4;++a) {
        vec4 xp = x, xm = x;
        xp.v[a] += h; xm.v[a] -= h;
        real vp = P.V(xp);
        real vm = P.V(xm);
        g.v[a] = (vp - vm) * (real(0.5)/h);
    }
    return g;
}

} // namespace rslm::deriv
