#pragma once
/**
 * RSLM Maths — diagnostics/overlay.hpp
 * ------------------------------------
 * Build overlay masks (polyline) from a sequence of vec4 points on an XY slice.
 * Also export path CSV.
 */

#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "linalg.hpp"
#include "grid.hpp"

namespace rslm::diag {

using rslm::linalg::vec4;

inline bool save_path_csv(const std::vector<vec4>& path, const std::string& file) {
    std::FILE* f = std::fopen(file.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "t,x,y,z\n");
    for (auto& p : path) {
        std::fprintf(f, "%.17g,%.17g,%.17g,%.17g\n",
                     double(p.v[0]), double(p.v[1]), double(p.v[2]), double(p.v[3]));
    }
    std::fclose(f);
    return true;
}

/** Build an overlay mask (same size as grid) from a polyline of points on XY. */
inline std::vector<std::uint8_t> mask_from_path_xy(const Grid2D& G,
                                                   const std::vector<vec4>& path,
                                                   int thickness = 1)
{
    std::vector<std::uint8_t> M(G.nx * G.ny, 0);

    auto to_ij = [&](double x, double y) -> std::pair<int,int> {
        int j = int(std::round((x - G.x0) / G.dx));
        int i = int(std::round((y - G.y0) / G.dy));
        return {i, j};
    };

    auto mark = [&](int i, int j) {
        if (i<0 || j<0 || i>=int(G.nx) || j>=int(G.ny)) return;
        M[i*G.ny + j] = 255;
    };

    auto draw_segment = [&](int i0,int j0,int i1,int j1) {
        int di = std::abs(i1-i0), dj = std::abs(j1-j0);
        int si = (i0<i1)?1:-1, sj=(j0<j1)?1:-1;
        int err = di - dj;
        int i=i0, j=j0;
        while (true) {
            for (int ti=-thickness; ti<=thickness; ++ti)
                for (int tj=-thickness; tj<=thickness; ++tj)
                    mark(i+ti, j+tj);
            if (i==i1 && j==j1) break;
            int e2 = 2*err;
            if (e2 > -dj) { err -= dj; i += si; }
            if (e2 <  di) { err += di; j += sj; }
        }
    };

    for (std::size_t k=1;k<path.size();++k) {
        auto [i0,j0] = to_ij(path[k-1].v[1], path[k-1].v[2]);
        auto [i1,j1] = to_ij(path[k].v[1],   path[k].v[2]);
        draw_segment(i0,j0,i1,j1);
    }
    return M;
}

} // namespace rslm::diag
