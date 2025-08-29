#pragma once
/**
 * RSLM Maths — audits/time_dilation.hpp
 * -------------------------------------
 * Compute time-dilation factor gamma from (g, u) without assuming Minkowski
 * coordinates. We infer a canonical future-directed timelike direction t by:
 *   - take the most negative eigenvalue/eigenvector of g (since signature -+++),
 *   - normalize t so that g(t,t) = -1,
 *   - set gamma = - g(u, t)  (component of u along the time axis in that frame),
 *   - spatial part w = u - gamma * t  (g-orthogonal to t),
 *   - v^2 = g(w,w) / gamma^2,  and v = sqrt(max(v^2, 0)).
 *
 * If u is properly normalized (g(u,u) = -1), identity holds: gamma = 1/sqrt(1-v^2).
 * This is a robust, frame-invariant way to report γ for diagnostics & windowing.
 */

#include <cmath>
#include <algorithm>
#include "config.hpp"
#include "linalg.hpp"
#include "quadform.hpp"
#include "metric.hpp"
#include "eigen_jacobi.hpp"  // symmetric 4×4 eigen (already in tree)

namespace rslm::audits {

using rslm::cfg::real;
using rslm::linalg::mat4;
using rslm::linalg::vec4;
using rslm::linalg::sym4;

struct TDReport {
    real gamma{1};
    real v_norm{0};   // |v| in [0,1)
    real q_u{0};      // g(u,u), should be -1 if normalized
};

// helper: pick timelike direction t with g(t,t) = -1
inline vec4 timelike_unit(const sym4& g) {
    // eigen-decompose g (symmetric): g = V diag(λ) V^T, V orthonormal (Euclidean)
    rslm::linalg::mat4 V;
    rslm::linalg::vec4 evals;
    rslm::eigen::jacobi_symmetric_4x4(g, V, evals); // function provided in tree

    // find most negative eigenvalue
    int k = 0;
    for (int i=1;i<4;++i) if (evals.v[i] < evals.v[k]) k = i;

    // unit Euclidean eigenvector (column k)
    vec4 e{ V.m[0][k], V.m[1][k], V.m[2][k], V.m[3][k] };

    // scale to timelike unit: g(e,e) = λ_k (since ||e||_2 = 1)
    real lam = evals.v[k];
    // lam should be < 0 for -+++; guard anyway
    real scale = real(1) / std::sqrt(std::max(real(1e-30), -lam));
    vec4 t{ e.v[0]*scale, e.v[1]*scale, e.v[2]*scale, e.v[3]*scale };
    return t; // now g(t,t) ≈ -1
}

inline TDReport time_dilation(const sym4& g, const vec4& u) {
    using rslm::quad::qform;
    TDReport R;
    R.q_u = qform(g, u);
    vec4 t = timelike_unit(g);

    // gamma = - g(u, t)
    real gut = real(0);
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) gut += u.v[i] * g.m[i][j] * t.v[j];
    R.gamma = -gut;

    // spatial projection w = u - gamma * t
    vec4 w{ u.v[0] - R.gamma*t.v[0],
            u.v[1] - R.gamma*t.v[1],
            u.v[2] - R.gamma*t.v[2],
            u.v[3] - R.gamma*t.v[3] };

    real gw = qform(g, w);               // should be >= 0
    real v2 = gw / std::max(real(1e-30), R.gamma*R.gamma);
    R.v_norm = std::sqrt(std::max(real(0), v2));
    return R;
}

} // namespace rslm::audits
