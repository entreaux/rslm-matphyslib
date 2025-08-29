#pragma once
/**
 * RSLM Maths — physics/einstein_fit.hpp
 * -------------------------------------
 * Build Einstein tensor G_{μν} from curvature and compute the Frobenius
 * residual || G - κ T ||_F for diagnostics.
 */

#include <cmath>
#include <algorithm>

#include "config.hpp"
#include "linalg.hpp"
#include "curvature.hpp"
#include "connection.hpp"
#include "field.hpp"
#include "stress_energy.hpp"

namespace rslm::phys {

using rslm::cfg::real;
using rslm::linalg::mat4;
using rslm::linalg::sym4;
using rslm::field::IMetricField;

inline mat4 to_mat4(const sym4& g) {
    mat4 M{};
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) M.m[i][j] = g.m[i][j];
    return M;
}

inline sym4 to_sym4(const mat4& A) {
    sym4 S{};
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) S.m[i][j] = real(0.5)*(A.m[i][j] + A.m[j][i]);
    return S;
}

// Einstein tensor G = Ricci - 0.5 g R  (all at x)
inline sym4 einstein_at(const IMetricField& F, const rslm::linalg::vec4& x) {
    auto Rm = rslm::curv::riemann_at(F, x);
    auto Rc = rslm::curv::ricci(Rm);
    auto g  = F.g(x);
    real R  = rslm::curv::scalar(rslm::conn::prepare_metric(F, x).g_inv, Rc);
    sym4 G{};
    for (int i=0;i<4;++i)
        for (int j=0;j<4;++j)
            G.m[i][j] = Rc.m[i][j] - real(0.5) * g.m[i][j] * R;
    return G;
}

// Frobenius norm of a 4×4 (no metric weighting; pure componentwise)
inline real frob(const mat4& A) {
    long double s=0;
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) { long double v=A.m[i][j]; s += v*v; }
    return static_cast<real>(std::sqrt(s));
}

/** ||G - κT||_F at x */
inline real residual_norm(const IMetricField& F,
                          const std::vector<Event>& evs,
                          const rslm::linalg::vec4& x,
                          const TSParams& P)
{
    sym4 Gs = einstein_at(F, x);
    mat4  G = to_mat4(Gs);
    mat4  T = stress_energy_at(F, evs, x, P);
    // R = G - κ T
    mat4 R{};
    for (int i=0;i<4;++i) for (int j=0;j<4;++j)
        R.m[i][j] = G.m[i][j] - P.kappa * T.m[i][j];
    return frob(R);
}

} // namespace rslm::phys
