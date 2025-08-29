#pragma once
/**
 * RSLM Maths — field.hpp
 * ----------------------
 * Metric-field and potential-field interfaces + a couple of baseline fields.
 * The metric interface returns a symmetric 4×4 metric g_{μν}(x).
 */

#include "config.hpp"
#include "linalg.hpp"
#include "metric.hpp"
#include "trace.hpp"

namespace rslm::field {

using rslm::cfg::real;
using rslm::linalg::vec4;
using rslm::linalg::mat4;
using rslm::linalg::sym4;

// ---------------- Interfaces ----------------
struct IMetricField {
    virtual ~IMetricField() = default;
    virtual sym4 g(const vec4& x) const = 0;
};

struct IPotential {
    virtual ~IPotential() = default;
    virtual real V(const vec4& x) const = 0;            // scalar potential
};

// --------------- Baseline fields ------------
struct MinkowskiField final : IMetricField {
    sym4 g(const vec4&) const override { return rslm::linalg::minkowski_eta(); }
};

// Very small, smooth “bump” curvature for testing (keeps signature via projection).
// g(x) = Proj( η + ε * diag( -e^{-r^2}, e^{-r^2}, e^{-r^2}, e^{-r^2} ) )
struct GaussianBumpField final : IMetricField {
    real eps;
    explicit GaussianBumpField(real epsilon = real(1e-2)) : eps(epsilon) {}
    sym4 g(const vec4& x) const override {
        sym4 eta = rslm::linalg::minkowski_eta();
        real r2 = x.v[0]*x.v[0] + x.v[1]*x.v[1] + x.v[2]*x.v[2] + x.v[3]*x.v[3];
        real s  = std::exp(-r2);
        mat4 tweak = rslm::linalg::diag(-eps*s, eps*s, eps*s, eps*s);
        sym4 gtry( rslm::linalg::mul(rslm::linalg::identity(), rslm::linalg::identity()) );
        // g = η + tweak (symmetrized by sym4 ctor), then project to Lorentzian
        mat4 tmp = rslm::linalg::mul(rslm::linalg::identity(), rslm::linalg::identity()); (void)tmp;
        sym4 gsum( rslm::linalg::mul(rslm::linalg::identity(), rslm::linalg::identity()) );
        // build: copy eta then add tweak on-diagonal
        sym4 g = eta;
        for (int i=0;i<4;++i) g.m[i][i] += tweak.m[i][i];
        return rslm::metric::project_signature(g);
    }
};

// Zero potential (for pure geodesics)
struct ZeroPotential final : IPotential {
    real V(const vec4&) const override { return real(0); }
};

// Smooth quadratic test potential: V = 0.5 * k * (x^2 + y^2 + z^2)  (ignores t)
struct RadialPotential final : IPotential {
    real k;
    explicit RadialPotential(real k_) : k(k_) {}
    real V(const vec4& x) const override {
        return real(0.5) * k * (x.v[1]*x.v[1] + x.v[2]*x.v[2] + x.v[3]*x.v[3]);
    }
};

} // namespace rslm::field
