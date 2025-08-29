#pragma once
/**
 * RSLM Maths — tetrad.hpp
 * -----------------------
 * Build tetrads E (columns are frame vectors e_a) so that:
 *     Eᵀ g E = η   with η = diag(-1,+1,+1,+1)
 *
 * Also returns a positive-definite proxy:
 *     g_tilde = E I Eᵀ = E Eᵀ
 */

#include "config.hpp"
#include "linalg.hpp"
#include "metric.hpp"
#include "eigen_jacobi.hpp"
#include "trace.hpp"
#include <array>
#include <cmath>

namespace rslm::tetrad {

using rslm::cfg::real;
using rslm::linalg::mat4;
using rslm::linalg::sym4;
using rslm::linalg::vec4;

// Construct tetrad E and PD proxy g_tilde. Returns Frobenius error ||Eᵀ g E - η||_F.
inline real build_tetrad(const sym4& g_in, mat4& E_out, mat4& g_tilde_out, real eps = real(1e-12)) {
    // Ensure correct signature
    sym4 g = metric::project_signature(g_in, std::max(real(1e-12), eps));

    // Eigen-decomp g = Q Λ Qᵀ
    mat4 Q; vec4 lam;
    eigen::jacobi_symmetric_4x4(g, Q, lam, 64, real(1e-14));

    // Find negative index; if none, project_signature would have created one.
    int neg_idx = 0;
    for (int i=0;i<4;++i) if (lam.v[i] < 0) { neg_idx = i; break; }

    // Construct D^{-1} = diag(1/sqrt(|λ_i|))
    mat4 Dinv;
    for (int i=0;i<4;++i) {
        real a = std::sqrt(std::max(std::fabs(lam.v[i]), eps));
        Dinv.m[i][i] = real(1)/a;
        for (int j=0;j<4;++j) if (i!=j) Dinv.m[i][j] = 0;
    }

    // E0 = Q * Dinv
    mat4 E0 = rslm::linalg::mul(Q, Dinv);

    // Permute columns so that the negative eigenvalue maps to column 0 (time-like)
    mat4 E = E0;
    if (neg_idx != 0) {
        for (int r=0;r<4;++r) std::swap(E.m[r][0], E.m[r][neg_idx]);
    }

    // Verify orthonormality: M = Eᵀ g E - η
    mat4 ET = rslm::linalg::transpose(E);
    mat4 M  = rslm::linalg::mul(rslm::linalg::mul(ET, g), E);

    // Subtract η
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) {
        real eta = (i==j) ? (i==0 ? real(-1) : real(1)) : real(0);
        M.m[i][j] -= eta;
    }

    // Frobenius norm
    real fro = 0;
    for (int r=0;r<4;++r) for (int c=0;c<4;++c) fro += M.m[r][c]*M.m[r][c];
    fro = std::sqrt(fro);

    // PD proxy: g_tilde = E Eᵀ
    g_tilde_out = rslm::linalg::mul(E, rslm::linalg::transpose(E));
    E_out = E;

    TRACE_INFO("tetrad_fro_error", fro);
    return fro;
}

// Quick PD check for g_tilde: all eivals > 0 (within eps)
inline bool is_pd(const mat4& Gt, real eps = real(1e-12)) {
    mat4 Q; vec4 lam;
    eigen::jacobi_symmetric_4x4(Gt, Q, lam, 64, real(1e-14));
    for (int i=0;i<4;++i) if (lam.v[i] <= eps) return false;
    return true;
}

} // namespace rslm::tetrad
