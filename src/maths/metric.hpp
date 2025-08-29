#pragma once
/**
 * RSLM Maths — metric.hpp
 * -----------------------
 * Metric utilities:
 *   - Minkowski η (−+++)
 *   - Construct g = Aᵀ η A        [guarantees Lorentzian signature]
 *   - Validate signature via eigen-inertia counting
 *   - Project a symmetric matrix to Lorentzian signature (nearest-ish)
 */

#include "config.hpp"
#include "linalg.hpp"
#include "eigen_jacobi.hpp"
#include "trace.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace rslm::metric {

using rslm::cfg::real;
using rslm::linalg::mat4;
using rslm::linalg::sym4;
using rslm::linalg::vec4;

// Minkowski η
inline sym4 minkowski_eta() { return rslm::linalg::minkowski_eta(); }

// Construct g = Aᵀ η A (ensures signature −+++ if A is invertible)
inline sym4 from_A(const mat4& A) {
    sym4 eta = minkowski_eta();
    mat4 tmp = rslm::linalg::mul(rslm::linalg::transpose(A), eta);
    mat4 g   = rslm::linalg::mul(tmp, A);
    return sym4(g);
}

// Count inertia (nneg, npos, nzero) using eigenvalues
inline void validate_signature(const sym4& g, int& nneg, int& npos, int& nzero, real eps = real(1e-10)) {
    mat4 Q; vec4 lam;
    eigen::jacobi_symmetric_4x4(g, Q, lam, 64, real(1e-14));
    nneg=npos=nzero=0;
    for (int i=0;i<4;++i) {
        if (lam.v[i] < -eps) ++nneg;
        else if (lam.v[i] > eps) ++npos;
        else ++nzero;
    }
    TRACE_INFO("signature_counts", "neg=" + std::to_string(nneg) + " pos=" + std::to_string(npos) + " zero=" + std::to_string(nzero));
}

// Project symmetric matrix to Lorentzian (−+++).
// Strategy: eigen-decomp g = Q Λ Qᵀ, then set signs to (−,+,+,+)
// keeping magnitudes |λ| but flooring by eps to avoid degeneracy.
// If there are 0 negatives, flip the smallest-|λ| positive to negative.
// If there are >1 negatives, keep the one with largest |λ| negative, flip others to positive.
inline sym4 project_signature(const sym4& g_in, real eps = real(1e-9)) {
    mat4 Q; vec4 lam;
    eigen::jacobi_symmetric_4x4(g_in, Q, lam, 64, real(1e-14));

    // Determine which index will be the single negative
    int neg_count=0, neg_idx=-1;
    std::array<real,4> absv{};
    for (int i=0;i<4;++i) {
        absv[i] = std::fabs(lam.v[i]);
        if (lam.v[i] < 0) { ++neg_count; if (neg_idx<0 || absv[i] > absv[neg_idx]) neg_idx=i; }
    }
    if (neg_count==0) {
        // flip the smallest |λ|
        int minabs = 0;
        for (int i=1;i<4;++i) if (absv[i] < absv[minabs]) minabs = i;
        neg_idx = minabs;
    } else {
        // keep the largest-magnitude negative; others become positive
    }

    // Build Λ' with desired signs and floor
    real L[4] = {0,0,0,0};
    for (int i=0;i<4;++i) {
        real a = std::max(absv[i], eps);
        L[i] = (i==neg_idx) ? -a : +a;
    }

    // Reconstruct Q Λ' Qᵀ
    mat4 D;
    for (int r=0;r<4;++r) for (int c=0;c<4;++c) D.m[r][c] = (r==c) ? L[r] : real(0);

    // g_proj = Q D Qᵀ
    mat4 tmp = rslm::linalg::mul(Q, D);
    mat4 gP  = rslm::linalg::mul(tmp, rslm::linalg::transpose(Q));

    sym4 out(gP);
    int nneg,npos,nzero;
    validate_signature(out, nneg, npos, nzero);
    if (!(nneg==1 && npos==3)) {
        TRACE_WARN("project_signature_warn", "inertia_after=" + std::to_string(nneg) + "," + std::to_string(npos) + "," + std::to_string(nzero));
    }
    return out;
}

} // namespace rslm::metric
