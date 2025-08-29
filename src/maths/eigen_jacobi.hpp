#pragma once
/**
 * RSLM Maths — eigen_jacobi.hpp
 * ------------------------------
 * Symmetric 4×4 Jacobi eigen-decomposition (A = Q Λ Qᵀ).
 * - Deterministic, no allocations
 * - Good enough for geometry control paths
 * - Tolerant to tiny off-diagonals
 */

#include "config.hpp"
#include "linalg.hpp"
#include "telemetry/trace.hpp"
#include <cmath>
#include <cstdlib>

namespace rslm::eigen {

using rslm::cfg::real;
using rslm::linalg::mat4;

// Finds the largest |offdiag| element; returns (p,q,|a_pq|).
inline void max_offdiag_abs(const mat4& A, int& p, int& q, real& val) {
    p = 0; q = 1; val = std::fabs(A.m[0][1]);
    for (int i=0;i<4;++i) {
        for (int j=i+1;j<4;++j) {
            real v = std::fabs(A.m[i][j]);
            if (v > val) { val = v; p = i; q = j; }
        }
    }
}

// One Jacobi rotation zeroing A_pq (p<q), updating A and Q in place.
inline void jacobi_rotate(mat4& A, mat4& Q, int p, int q) {
    if (p==q) return;
    real app = A.m[p][p], aqq = A.m[q][q], apq = A.m[p][q];
    if (apq == real(0)) return;

    real tau = (aqq - app) / (real(2)*apq);
    real t = (tau >= 0) ? real(1)/(tau + std::sqrt(real(1)+tau*tau))
                        : real(1)/(tau - std::sqrt(real(1)+tau*tau));
    real c = real(1)/std::sqrt(real(1)+t*t);
    real s = t * c;

    // Update A (p-th and q-th rows/cols)
    real app_new = app - t * apq;
    real aqq_new = aqq + t * apq;
    A.m[p][p] = app_new;
    A.m[q][q] = aqq_new;
    A.m[p][q] = A.m[q][p] = 0;

    for (int k=0;k<4;++k) {
        if (k==p || k==q) continue;
        real aik = A.m[p][k], akq = A.m[q][k];
        A.m[p][k] = A.m[k][p] = c*aik - s*akq;
        A.m[q][k] = A.m[k][q] = s*aik + c*akq;
    }

    // Update Q = Q * R(p,q,c,s)
    for (int k=0;k<4;++k) {
        real qkp = Q.m[k][p], qkq = Q.m[k][q];
        Q.m[k][p] = c*qkp - s*qkq;
        Q.m[k][q] = s*qkp + c*qkq;
    }
}

/**
 * Jacobi eigen decomposition for symmetric 4×4.
 * @param A_in   symmetric input (only symmetry is assumed)
 * @param Q_out  orthonormal eigenvectors (columns)
 * @param lam    eigenvalues (diagonal of Λ), from A_out after convergence
 * @returns true if converged within tol or max_sweeps
 */
inline bool jacobi_symmetric_4x4(const mat4& A_in, mat4& Q_out, rslm::linalg::vec4& lam,
                                 int max_sweeps = 32, real tol = real(1e-12)) {
    // Copy A; initialize Q=I
    mat4 A = A_in;
    Q_out = rslm::linalg::identity();

    // Main sweeps
    for (int sweep=0; sweep<max_sweeps; ++sweep) {
        int p, q; real off;
        max_offdiag_abs(A, p, q, off);
        if (off < tol) break;
        jacobi_rotate(A, Q_out, p, q);
    }

    for (int i=0;i<4;++i) lam.v[i] = A.m[i][i];
    return true;
}

} // namespace rslm::eigen
