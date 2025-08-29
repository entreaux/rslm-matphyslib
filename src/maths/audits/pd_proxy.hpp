#pragma once
/**
 * RSLM Maths — audits/pd_proxy.hpp
 * --------------------------------
 * Build a simple positive-definite proxy metric \tilde g from a Lorentzian g.
 * We use \tilde g = g^T g (for symmetric g this is g^2), which is symmetric PD
 * as long as g is non-singular. This is adequate as an optimizer inner product.
 * Also includes a tiny 4×4 Cholesky check for PD verification.
 *
 * NOTE: This is *diagnostic* PD. For production we might switch to tetrad-based
 * \tilde g = e^T I e after we finalize our tetrad builder; the wrapper names are
 * kept generic so it’s easy to swap later.
 */

#include <cmath>
#include <algorithm>
#include "config.hpp"
#include "linalg.hpp"
#include "metric.hpp"

namespace rslm::audits {

using rslm::cfg::real;
using rslm::linalg::mat4;
using rslm::linalg::sym4;

// Build PD proxy: \tilde g = g^2
inline mat4 pd_proxy_square(const sym4& g) {
    // g is symmetric; \tilde g = g * g is symmetric and (semi)PD.
    return rslm::linalg::mul(g, g);
}

// Simple 4×4 Cholesky (A = L L^T). Returns false if not PD.
inline bool cholesky4(const mat4& A, mat4& L, real eps = real(0)) {
    L = mat4{};
    for (int i=0;i<4;++i) {
        for (int j=0;j<=i;++j) {
            real sum = A.m[i][j];
            for (int k=0;k<j;++k) sum -= L.m[i][k] * L.m[j][k];

            if (i == j) {
                if (sum <= eps) return false;
                L.m[i][j] = std::sqrt(sum);
            } else {
                L.m[i][j] = sum / L.m[j][j];
            }
        }
    }
    return true;
}

struct PDReport {
    bool ok{false};
    real min_diagL{0};
    real max_diagL{0};
};

// Check PD by Cholesky and return a tiny report (for logs)
inline PDReport check_pd(const mat4& A) {
    mat4 L;
    PDReport R;
    R.ok = cholesky4(A, L);
    if (R.ok) {
        R.min_diagL = L.m[0][0];
        R.max_diagL = L.m[0][0];
        for (int i=0;i<4;++i) {
            R.min_diagL = std::min(R.min_diagL, L.m[i][i]);
            R.max_diagL = std::max(R.max_diagL, L.m[i][i]);
        }
    }
    return R;
}

} // namespace rslm::audits
