#pragma once
/**
 * RSLM Maths — diagnostics/ppm.hpp
 * --------------------------------
 * Save a Grid2D as a color PPM (ASCII P3). Pure text, easy to diff.
 * You can pass an optional overlay mask (same dims) where non-zero pixels are drawn black.
 */

#include <string>
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

#include "grid.hpp"
#include "palette.hpp"

namespace rslm::diag {

/**
 * Save heatmap using a palette functor P::map(double)->RGB.
 * If overlay.size()==nx*ny and overlay[i*ny+j]!=0, pixel is drawn black.
 * vmin/vmax: if equal or inverted, we auto-compute from data.
 */
template <typename Palette>
inline bool save_ppm(const Grid2D& G, const std::string& path,
                     double vmin = NAN, double vmax = NAN,
                     const std::vector<std::uint8_t>& overlay = {})
{
    // auto-range if needed
    if (!(vmin < vmax)) {
        vmin = G.val.empty() ? 0.0 : G.val[0];
        vmax = vmin;
        for (auto v : G.val) { if (v < vmin) vmin = v; if (v > vmax) vmax = v; }
        if (vmin == vmax) { vmin -= 1.0; vmax += 1.0; }
    }

    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "P3\n%zu %zu\n255\n", G.ny, G.nx); // width=ny, height=nx

    auto norm = [&](double v) {
        double t = (v - vmin) / (vmax - vmin);
        if (t < 0) t = 0; if (t > 1) t = 1;
        return t;
    };

    for (std::size_t i=0;i<G.nx;++i) {
        for (std::size_t j=0;j<G.ny;++j) {
            if (!overlay.empty() && overlay[i*G.ny + j]) {
                std::fprintf(f, "0 0 0 ");
            } else {
                double t = norm(G.at(i,j));
                RGB c = Palette::map(t);
                std::fprintf(f, "%d %d %d ", int(c.r), int(c.g), int(c.b));
            }
        }
        std::fprintf(f, "\n");
    }
    std::fclose(f);
    return true;
}

} // namespace rslm::diag
