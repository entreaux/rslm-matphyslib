#pragma once
/**
 * RSLM Maths — curvature.hpp
 * --------------------------
 * Riemann tensor R^μ_{ναβ}, Ricci R_{αβ}, and scalar R.
 * Derivatives of Γ are computed by finite differences over space-time.
 */

#include "config.hpp"
#include "linalg.hpp"
#include "connection.hpp"
#include "deriv.hpp"
#include "trace.hpp"
#include "units.hpp"

namespace rslm::curv {

using rslm::cfg::real;
using rslm::linalg::vec4;
using rslm::linalg::mat4;
using rslm::conn::Gamma;
using rslm::conn::MetricPack;
using rslm::field::IMetricField;

// ∂_a Γ^μ_{νβ} via central difference on the metric field (rebuild Γ at x±h e_a)
inline Gamma dGamma_dir(const IMetricField& F, const vec4& x, int a, real h = rslm::units::C().fd_h) {
    vec4 xp = x, xm = x; xp.v[a]+=h; xm.v[a]-=h;
    MetricPack Pp = rslm::conn::prepare_metric(F, xp);
    MetricPack Pm = rslm::conn::prepare_metric(F, xm);
    Gamma Gp = rslm::conn::christoffel(Pp);
    Gamma Gm = rslm::conn::christoffel(Pm);

    Gamma d{}; const real s = real(0.5)/h;
    for (int mu=0; mu<4; ++mu)
        for (int nu=0; nu<4; ++nu)
            for (int b=0; b<4; ++b)
                d.G[mu][nu][b] = (Gp.G[mu][nu][b] - Gm.G[mu][nu][b]) * s;
    return d;
}

struct Riemann {
    // R[mu][nu][a][b] = R^μ_{ναβ}
    real R[4][4][4][4]{};
};

inline Riemann riemann_at(const IMetricField& F, const vec4& x) {
    MetricPack M = rslm::conn::prepare_metric(F, x);
    Gamma G = rslm::conn::christoffel(M);

    // Precompute ∂_a Γ
    Gamma dG[4];
    for (int a=0;a<4;++a) dG[a] = dGamma_dir(F, x, a);

    Riemann out{};
    for (int mu=0; mu<4; ++mu)
    for (int nu=0; nu<4; ++nu)
    for (int a=0;  a<4; ++a)
    for (int b=0;  b<4; ++b) {
        real s = dG[a].G[mu][nu][b] - dG[b].G[mu][nu][a];
        // + Γ^μ_{σα} Γ^σ_{νβ} - Γ^μ_{σβ} Γ^σ_{να}
        for (int sig=0; sig<4; ++sig) {
            s += G.G[mu][sig][a] * G.G[sig][nu][b];
            s -= G.G[mu][sig][b] * G.G[sig][nu][a];
        }
        out.R[mu][nu][a][b] = s;
    }
    return out;
}

// Ricci: R_{αβ} = R^μ_{αμβ}
inline mat4 ricci(const Riemann& R) {
    mat4 Rc;
    for (int a=0;a<4;++a)
        for (int b=0;b<4;++b) {
            real s=0;
            for (int mu=0;mu<4;++mu) s += R.R[mu][a][mu][b];
            Rc.m[a][b]=s;
        }
    return Rc;
}

// Scalar: R = g^{αβ} R_{αβ}
inline real scalar(const mat4& g_inv, const mat4& Ric) {
    real s=0;
    for (int a=0;a<4;++a) for (int b=0;b<4;++b) s += g_inv.m[a][b] * Ric.m[a][b];
    return s;
}

// Frobenius norm ||R||_F
inline real frob_riemann(const Riemann& R) {
    long double S=0;
    for (int i=0;i<4;++i)
    for (int j=0;j<4;++j)
    for (int k=0;k<4;++k)
    for (int l=0;l<4;++l) {
        long double v = R.R[i][j][k][l];
        S += v*v;
    }
    return static_cast<real>(std::sqrt(S));
}

} // namespace rslm::curv
