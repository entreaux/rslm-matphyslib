# rslm-matlib
RSLM MatPhysLib

A compact C++20 maths + physics library for Relativistic Spacetime Language Models (RSLM).
It is a learnable Lorentzian geometry toolbelt (metrics, connections, curvature, geodesics), physics terms (stress–energy, Einstein “fit” residual), and diagnostics (heatmaps, overlays, OBJ surfaces, CSV paths, rich logs).
CPU-only and tuned for Apple Silicon (M-series). Works great as an internal API for future RSLM iterations.

Highlights

Lorentzian core: (-+++) metrics, Christoffels, Riemann/Ricci/scalar curvature, timelike normalization, tetrads and a PD proxy metric for optimizers.

Geodesic integrators: stable Velocity-Verlet step with optional potential forcing, chain-safe helpers.

Physics terms: semantic stress–energy tensor builder and Einstein-residual fit utilities for diagnostics.

Telemetry: zero-dependency structured logger with run IDs, throttling, and plain-text logs.

Diagnostics suite: grid samplers, palettes, PPM export, OBJ surface export, 2D overlays, path masks.

Unified facade: include rslmmatlib.hpp to reach the public API without wading through submodules.

Module Map
src/maths/
  config.hpp            # central knobs & compile-time toggles
  units.hpp             # c, epsilons, finite-diff steps, dtau defaults
  numeric.hpp           # safe sqrt, saturating tanh, comparisons, etc.
  linalg.hpp            # tiny fixed 4D vectors/matrices, mat4/sym4 ops
  quadform.hpp          # g(u,u), mixed forms, raising/lowering
  metric.hpp            # signature checks, tetrads, PD proxy metric
  connection.hpp        # Γ (Christoffel), metric packs
  deriv.hpp             # finite differences on fields/potentials
  curvature.hpp         # Riemann, Ricci, scalar curvature
  field.hpp             # IMetricField, IPotential interfaces
  integrators.hpp       # velocity-Verlet geodesic step, helpers
  physics/
    stress_energy.hpp   # semantic T_{μν}(x) builder (RBF + energy/mass)
    einstein_fit.hpp    # G_{μν} - κ T_{μν} diagnostics & samplers
  diagnostics/
    grid.hpp            # grid generation & sampling utilities
    palette.hpp         # color maps (Thermal5, etc.)
    ppm.hpp             # PPM writer w/ optional mask overlay
    overlay.hpp         # path masks & compositing helpers
    export.hpp          # OBJ surface exporter, CSV path writer
  telemetry/
    logger.hpp/.cpp     # plain-text structured logger (run_id, levels)
    trace.hpp           # scope-based tracing + throttling
  audits/
    pd_proxy.hpp        # PD metric proxy & Cholesky checks
    time_dilation.hpp   # γ(v) helpers and sanity tests

rslmmatlib.hpp           # unified facade 


Configuration Knobs (config.hpp / units.hpp)

Centralized constants used throughout (good for quick hyper-tuning):
Units / numerics
constexpr real c = 1;
constexpr real g_eps = 1e-12; (metric tolerance)
constexpr real fd_h = 1e-4; (finite-difference step)
constexpr real dtau = 5e-1; (typical geodesic step

Telemetry
constexpr int flush_every = 200;
constexpr bool console_mirror_default = true;

Sampling / export
constexpr const char* out_dir = "out";

Overlay mask thickness defaults in diagnostics/overlay.hpp
Stress–energy / Einstein fit
RBF width sigma, eta stabilizer, κ coupling (see physics/*.hpp)
Tip: keep dtau small for stability, especially with strong curvature or forcing.

Design Notes
Signature-safe: All metrics are checked for the intended (-+++) signature; tetrads give a PD proxy for optimizer norms while keeping true Lorentzian geometry for retractions.
Stable geodesics: Velocity-Verlet-like step with timelike shell renormalization guards against drift.
Einstein “fit”: We do not solve Einstein’s equations; we minimize a residual (G_{μν} - κ T_{μν}) to use curvature as a capacity/saliency diagnostic, matching the RSLM blueprint.
