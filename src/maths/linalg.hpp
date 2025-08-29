#pragma once
/**
 * RSLM Maths — linalg.hpp
 * -----------------------
 * Minimal core algebra for 4D spacetime:
 *   - vec4  : 4-vector (contravariant or covariant, depending on use)
 *   - mat4  : 4×4 matrix, row-major
 *   - sym4  : 4×4 symmetric matrix (constructed by symmetrization)
 *
 * Operations:
 *   - identity(), diag(a0,a1,a2,a3)
 *   - transpose(M), mul(M,N), mul(M,v)
 *   - det(M) with partial pivoting (tracks permutation sign)
 *   - inverse(M, out_inv, out_det, out_cond_inf) → bool success
 *   - norm_inf(M), condition number (∞-norm) via inv
 *
 * Notes:
 *   - Designed for robustness over performance (we’ll NEON-opt later).
 *   - Logging is throttled: it prints matrices/vectors only when
 *     RSLM_TRACE_HEAVY=1 is set, and even then it emits few samples.
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <ostream>
#include <string>
#include <algorithm>

#include "config.hpp"
#include "trace.hpp"

namespace rslm::linalg {

using rslm::cfg::real;

// ----------------------------------------------------------------------------
// Helpers: heavy-trace gate + throttled emitter
// ----------------------------------------------------------------------------
inline bool heavy_trace_enabled() {
    static int v = []{
        if (const char* e = std::getenv("RSLM_TRACE_HEAVY")) return (std::string(e) != "0") ? 1 : 0;
        return 0;
    }();
    return v != 0;
}
struct Throttle {
    std::uint64_t n{0};
    std::uint64_t first{4};
    std::uint64_t stride{1000};
    bool tick() { ++n; return (n <= first) || (stride && (n % stride == 0)); }
};

// ----------------------------------------------------------------------------
// vec4
// ----------------------------------------------------------------------------
struct vec4 {
    std::array<real,4> v{real(0),real(0),real(0),real(0)};

    vec4() = default;
    vec4(real t, real x, real y, real z) { v[0]=t; v[1]=x; v[2]=y; v[3]=z; }
    explicit vec4(std::initializer_list<real> a) {
        std::size_t i=0; for (real e : a) { if (i<4) v[i++]=e; }
        for (; i<4; ++i) v[i]=real(0);
    }

    real& operator[](std::size_t i)       { return v[i]; }
    const real& operator[](std::size_t i) const { return v[i]; }
};

// Stream pretty-printer (compact)
inline std::ostream& operator<<(std::ostream& os, const vec4& a) {
    os.setf(std::ios::fixed, std::ios::floatfield);
    os << "[" << double(a.v[0]) << ", " << double(a.v[1]) << ", "
       << double(a.v[2]) << ", " << double(a.v[3]) << "]";
    return os;
}

// ----------------------------------------------------------------------------
// mat4 (row-major)
// ----------------------------------------------------------------------------
struct mat4 {
    // m[r][c]
    std::array<std::array<real,4>,4> m{{
        {{real(0),real(0),real(0),real(0)}},
        {{real(0),real(0),real(0),real(0)}},
        {{real(0),real(0),real(0),real(0)}},
        {{real(0),real(0),real(0),real(0)}},
    }};

    mat4() = default;
    explicit mat4(real s) {
        for (int r=0;r<4;++r) for (int c=0;c<4;++c) m[r][c]=(r==c)?s:real(0);
    }
    explicit mat4(std::initializer_list<real> a) {
        int k=0;
        for (real e : a) { m[k/4][k%4]=e; if (++k==16) break; }
        for (;k<16;++k) m[k/4][k%4]=real(0);
    }

    real*       operator[](std::size_t r)       { return m[r].data(); }
    const real* operator[](std::size_t r) const { return m[r].data(); }
};

inline std::ostream& operator<<(std::ostream& os, const mat4& A) {
    os.setf(std::ios::fixed, std::ios::floatfield);
    os << "[";
    for (int r=0;r<4;++r) {
        os << (r==0?"[":" [");
        for (int c=0;c<4;++c) {
            os << double(A.m[r][c]);
            if (c<3) os << ", ";
        }
        os << (r<3?"],\n":"]");
    }
    os << "]";
    return os;
}

// Symmetric wrapper: builds 0.5*(M + M^T)
struct sym4 : mat4 {
    sym4() = default;
    explicit sym4(const mat4& A) {
        for (int r=0;r<4;++r) for (int c=0;c<4;++c) {
            m[r][c] = real(0.5)*(A.m[r][c] + A.m[c][r]);
        }
    }
    static sym4 from_diag(real a0, real a1, real a2, real a3) {
        sym4 S; for (int i=0;i<4;++i) for (int j=0;j<4;++j) S.m[i][j]=real(0);
        S.m[0][0]=a0; S.m[1][1]=a1; S.m[2][2]=a2; S.m[3][3]=a3;
        return S;
    }
};

// ----------------------------------------------------------------------------
// Constructors/utilities
// ----------------------------------------------------------------------------
inline mat4 identity() {
    return mat4{real(1),0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
}
inline mat4 diag(real a0, real a1, real a2, real a3) {
    mat4 A; for (int i=0;i<4;++i) for (int j=0;j<4;++j) A.m[i][j]=real(0);
    A.m[0][0]=a0; A.m[1][1]=a1; A.m[2][2]=a2; A.m[3][3]=a3; return A;
}
// Minkowski metric η with signature (-,+,+,+)
inline sym4 minkowski_eta() {
    return sym4::from_diag(real(-1), real(1), real(1), real(1));
}

inline mat4 transpose(const mat4& A) {
    mat4 T;
    for (int r=0;r<4;++r) for (int c=0;c<4;++c) T.m[r][c]=A.m[c][r];
    return T;
}

inline mat4 mul(const mat4& A, const mat4& B) {
    mat4 C;
    for (int r=0;r<4;++r) {
        for (int c=0;c<4;++c) {
            real s=0;
            for (int k=0;k<4;++k) s += A.m[r][k]*B.m[k][c];
            C.m[r][c]=s;
        }
    }
    return C;
}

inline vec4 mul(const mat4& A, const vec4& x) {
    vec4 y;
    for (int r=0;r<4;++r) {
        real s=0;
        for (int k=0;k<4;++k) s += A.m[r][k]*x.v[k];
        y.v[r]=s;
    }
    return y;
}

// ∞-norm (max row-sum of abs)
inline real norm_inf(const mat4& A) {
    real best = 0;
    for (int r=0;r<4;++r) {
        real s = 0;
        for (int c=0;c<4;++c) s += std::fabs(A.m[r][c]);
        if (s>best) best=s;
    }
    return best;
}

// ----------------------------------------------------------------------------
// Determinant with partial pivoting (also used during inversion)
// ----------------------------------------------------------------------------
inline real det(mat4 A) {
    Throttle th;
    if (heavy_trace_enabled() && th.tick()) TRACE_DEBUG("det_A", A);

    int sign = 1;
    real det = 1;
    for (int k=0;k<4;++k) {
        int piv = k;
        real amax = std::fabs(A.m[k][k]);
        for (int r=k+1;r<4;++r) {
            real v = std::fabs(A.m[r][k]);
            if (v > amax) { amax = v; piv = r; }
        }
        if (amax == real(0)) return real(0);
        if (piv != k) { std::swap(A.m[piv], A.m[k]); sign = -sign; }
        real akk = A.m[k][k];
        det *= akk;
        for (int r=k+1;r<4;++r) {
            real f = A.m[r][k] / akk;
            for (int c=k;c<4;++c) A.m[r][c] -= f * A.m[k][c];
        }
    }
    real out = det * real(sign);
    if (heavy_trace_enabled() && th.tick()) TRACE_DEBUG("det", out);
    return out;
}

// Inverse (Gauss–Jordan 4x4)
inline bool inverse(const mat4& A, mat4& Ainv, real& out_det, real& out_cond_inf, real eps = real(1e-14)) {
    Throttle th;
    if (heavy_trace_enabled() && th.tick()) TRACE_DEBUG("inv_A", A);

    real aug[4][8];
    for (int r=0;r<4;++r) {
        for (int c=0;c<4;++c) aug[r][c] = A.m[r][c];
        for (int c=0;c<4;++c) aug[r][4+c] = (r==c)?real(1):real(0);
    }

    int sign = 1;
    out_det = 1;

    for (int k=0;k<4;++k) {
        int piv = k;
        real amax = std::fabs(aug[k][k]);
        for (int r=k+1;r<4;++r) {
            real v = std::fabs(aug[r][k]);
            if (v > amax) { amax = v; piv = r; }
        }
        if (amax < eps) { return false; }

        if (piv != k) {
            for (int c=0;c<8;++c) std::swap(aug[piv][c], aug[k][c]);
            sign = -sign;
        }

        real akk = aug[k][k];
        out_det *= akk;

        real inv_akk = real(1) / akk;
        for (int c=0;c<8;++c) aug[k][c] *= inv_akk;

        for (int r=k+1;r<4;++r) {
            real f = aug[r][k];
            if (f==real(0)) continue;
            for (int c=0;c<8;++c) aug[r][c] -= f * aug[k][c];
        }
    }

    out_det *= real(sign);

    for (int k=3;k>=0;--k) {
        for (int r=0;r<k;++r) {
            real f = aug[r][k];
            if (f==real(0)) continue;
            for (int c=0;c<8;++c) aug[r][c] -= f * aug[k][c];
        }
    }

    for (int r=0;r<4;++r) for (int c=0;c<4;++c) Ainv.m[r][c] = aug[r][4+c];

    real nA = norm_inf(A);
    real nInv = norm_inf(Ainv);
    out_cond_inf = nA * nInv;

    if (heavy_trace_enabled() && th.tick()) {
        TRACE_DEBUG("inv_det", out_det);
        TRACE_DEBUG("inv_cond_inf", out_cond_inf);
        TRACE_DEBUG("Ainv", Ainv);
    }
    return true;
}

} // namespace rslm::linalg
