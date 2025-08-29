#pragma once
/**
 * RSLM Maths — connection.hpp
 * ---------------------------
 * Christoffel symbols Γ^μ_{αβ} (Levi–Civita) from a metric field via
 * finite differences ∂_a g_{μν}.
 */

#include "config.hpp"
#include "linalg.hpp"
#include "quadform.hpp"
#include "deriv.hpp"
#include "trace.hpp"

namespace rslm::conn {

using rslm::cfg::real;
using rslm::linalg::vec4;
using rslm::linalg::mat4;
using rslm::linalg::sym4;
using rslm::deriv::DMetric4;
using rslm::field::IMetricField;

struct Gamma {
    // G[mu][a][b] = Γ^μ_{αβ}
    real G[4][4][4]{};
};

struct MetricPack {
    sym4 g;
    mat4 g_inv;
    DMetric4 dg; // ∂_a g
    bool inv_ok{false};
};

// Invert g with robust fallback tolerance
inline MetricPack prepare_metric(const IMetricField& F, const vec4& x) {
    MetricPack P;
    P.g = F.g(x);

    // inverse
    rslm::linalg::mat4 inv;
    rslm::cfg::real det=0, cond=0;
    P.inv_ok = rslm::linalg::inverse(P.g, inv, det, cond, rslm::cfg::real(1e-14));
    P.g_inv = inv;

    // ∂g
    P.dg = rslm::deriv::dmetric4(F, x);
    return P;
}

inline Gamma christoffel(const MetricPack& M) {
    Gamma out{};
    // Γ^μ_{αβ} = 1/2 g^{μν} ( ∂_α g_{νβ} + ∂_β g_{να} - ∂_ν g_{αβ} )
    for (int mu=0; mu<4; ++mu) {
        for (int a=0; a<4; ++a) {
            for (int b=0; b<4; ++b) {
                real s = 0;
                for (int nu=0; nu<4; ++nu) {
                    real term = M.dg.dg[a].m[nu][b] + M.dg.dg[b].m[nu][a] - M.dg.dg[nu].m[a][b];
                    s += M.g_inv.m[mu][nu] * term;
                }
                out.G[mu][a][b] = real(0.5) * s;
            }
        }
    }
    return out;
}

} // namespace rslm::conn
