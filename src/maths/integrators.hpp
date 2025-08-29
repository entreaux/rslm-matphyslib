#pragma once
/**
 * RSLM Maths — integrators.hpp
 * ----------------------------
 * Simple, stable one-step integrators for geodesic motion with optional force.
 * We expose a velocity-Verlet–like step that keeps good energy behavior.
 */

#include "config.hpp"
#include "linalg.hpp"
#include "connection.hpp"
#include "deriv.hpp"
#include "field.hpp"
#include "quadform.hpp"
#include "units.hpp"
#include "trace.hpp"

namespace rslm::integ {

using rslm::cfg::real;
using rslm::linalg::vec4;
using rslm::conn::MetricPack;
using rslm::conn::Gamma;
using rslm::field::IMetricField;
using rslm::field::IPotential;

// Acceleration a^μ = - Γ^μ_{αβ} u^α u^β  +  f^μ, with f^μ = - g^{μν} ∂_ν V
inline vec4 accel(const MetricPack& M, const Gamma& G, const vec4& u, const IPotential* P, const vec4& x) {
    vec4 a{0,0,0,0};
    // - Γ term
    for (int mu=0; mu<4; ++mu) {
        real s = 0;
        for (int a_=0;a_<4;++a_) for (int b=0;b<4;++b) s += G.G[mu][a_][b] * u.v[a_]*u.v[b];
        a.v[mu] = -s;
    }
    // + forcing
    if (P) {
        vec4 gV = rslm::deriv::gradV(*P, x);
        for (int mu=0; mu<4; ++mu) {
            real s=0; for (int nu=0;nu<4;++nu) s += M.g_inv.m[mu][nu] * gV.v[nu];
            a.v[mu] += -s;
        }
    }
    return a;
}

// Project a timelike 4-velocity back onto the shell g(u,u) = -1 (c=1)
inline void renormalize_timelike(const rslm::linalg::mat4& g, vec4& u) {
    using rslm::quad::qform;
    real s = qform(g, u);
    // s should be -1; scale by 1/sqrt(-s)
    if (s >= real(0)) return; // not timelike; leave as-is (caller may decide)
    real k = real(1) / std::sqrt(-s);
    for (int i=0;i<4;++i) u.v[i] *= k;
}

// Velocity-Verlet style step (geometric-ish), small dtau advised.
inline void geodesic_step(const IMetricField& F, const IPotential* P,
                          vec4& x, vec4& u, real dtau) {
    MetricPack M = rslm::conn::prepare_metric(F, x);
    Gamma G = rslm::conn::christoffel(M);

    vec4 a0 = accel(M, G, u, P, x);

    // half-step velocity
    vec4 uh = u;
    for (int i=0;i<4;++i) uh.v[i] += real(0.5) * dtau * a0.v[i];

    // position update
    for (int i=0;i<4;++i) x.v[i] += dtau * uh.v[i];

    // recompute at new position
    MetricPack M1 = rslm::conn::prepare_metric(F, x);
    Gamma G1 = rslm::conn::christoffel(M1);
    vec4 a1 = accel(M1, G1, uh, P, x);

    // finish velocity
    for (int i=0;i<4;++i) u.v[i] = uh.v[i] + real(0.5) * dtau * a1.v[i];

    // keep timelike normalization (for stability)
    renormalize_timelike(M1.g, u);
}

inline void rk4_geodesic(const rslm::field::IMetricField& F,
                         rslm::linalg::vec4& x,
                         rslm::linalg::vec4& u,
                         rslm::cfg::real dtau,
                         const rslm::field::IPotential* P = nullptr) {
    // For now, route to the stable velocity-Verlet step.
    geodesic_step(F, P, x, u, dtau);
}

} // namespace rslm::integ
