#pragma once
/**
 * RSLM Maths — quadform.hpp
 * -------------------------
 * Quadratic forms and index-raising/lowering helpers for 4D tensors.
 *
 * Conventions:
 * - Inputs are plain mat4/vec4; when used as a metric g, we expect signature (-,+,+,+)
 *   for Minkowski unless stated otherwise.
 */

#include "config.hpp"
#include "linalg.hpp"
#include "telemetry/trace.hpp"
#include <cmath>

namespace rslm::quad {

using rslm::cfg::real;
using rslm::linalg::mat4;
using rslm::linalg::vec4;

// vᵀ g v
inline real qform(const mat4& g, const vec4& v) {
    real s = 0;
    for (int r=0;r<4;++r) {
        real row = 0;
        for (int c=0;c<4;++c) row += g.m[r][c] * v.v[c];
        s += v.v[r] * row;
    }
    TRACE_DEBUG("qform", s);
    return s;
}

// Lower index: v_μ = g_{μν} v^ν
inline vec4 lower(const mat4& g, const vec4& v_contra) {
    vec4 res;
    for (int r=0;r<4;++r) {
        real s = 0;
        for (int c=0;c<4;++c) s += g.m[r][c] * v_contra.v[c];
        res.v[r] = s;
    }
    TRACE_DEBUG("lower", res);
    return res;
}

// Raise index: v^μ = g^{μν} v_ν   (requires g_inv)
inline vec4 raise(const mat4& g_inv, const vec4& v_cov) {
    vec4 res;
    for (int r=0;r<4;++r) {
        real s = 0;
        for (int c=0;c<4;++c) s += g_inv.m[r][c] * v_cov.v[c];
        res.v[r] = s;
    }
    TRACE_DEBUG("raise", res);
    return res;
}

} // namespace rslm::quad
