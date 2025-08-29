#pragma once
/**
 * RSLM Maths — diagnostics/export.hpp
 * -----------------------------------
 * Export Grid2D to:
 *  - CSV   : x,y,value (one row per cell)
 *  - ASCII : matrix style (ny columns per line), good for quick diffs
 *  - OBJ   : 3D surface with z = scale * value, vertices laid on (x,y)
 *
 * All files are plain text. No JSON.
 */

#include <cstdio>
#include <string>
#include <filesystem>

#include "config.hpp"
#include "trace.hpp"
#include "grid.hpp"

namespace rslm::diag {

inline void ensure_dir(const std::filesystem::path& p) {
    std::error_code ec;
    std::filesystem::create_directories(p, ec);
}

inline bool save_csv(const Grid2D& G, const std::string& path) {
    ensure_dir(std::filesystem::path(path).parent_path());
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "x,y,value\n");
    for (std::size_t i=0;i<G.nx;++i) {
        for (std::size_t j=0;j<G.ny;++j) {
            double x = double(G.x0 + j*G.dx);
            double y = double(G.y0 + i*G.dy);
            double v = double(G.at(i,j));
            std::fprintf(f, "%.17g,%.17g,%.17g\n", x, y, v);
        }
    }
    std::fclose(f);
    TRACE_INFO("save_csv", path);
    return true;
}

inline bool save_ascii(const Grid2D& G, const std::string& path) {
    ensure_dir(std::filesystem::path(path).parent_path());
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    for (std::size_t i=0;i<G.nx;++i) {
        for (std::size_t j=0;j<G.ny;++j) {
            double v = double(G.at(i,j));
            std::fprintf(f, (j+1<G.ny) ? "%.6e " : "%.6e", v);
        }
        std::fprintf(f, "\n");
    }
    std::fclose(f);
    TRACE_INFO("save_ascii", path);
    return true;
}

/**
 * Export as a triangulated surface: vertices on (x,y), height = scale * value.
 * - Vertices: (x, y, z) with z = scale*value
 * - Faces: two triangles per quad
 */
inline bool save_obj_surface(const Grid2D& G, const std::string& path, double scale = 1.0) {
    ensure_dir(std::filesystem::path(path).parent_path());
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;

    std::fprintf(f, "# RSLM surface OBJ: z = %g * value\n", scale);
    // vertices
    for (std::size_t i=0;i<G.nx;++i) {
        for (std::size_t j=0;j<G.ny;++j) {
            double x = double(G.x0 + j*G.dx);
            double y = double(G.y0 + i*G.dy);
            double z = scale * double(G.at(i,j));
            std::fprintf(f, "v %.17g %.17g %.17g\n", x, y, z);
        }
    }
    // simple grid faces (1-based indices)
    auto vid = [&](std::size_t i, std::size_t j) { return int(i*G.ny + j + 1); };
    for (std::size_t i=0;i+1<G.nx;++i) {
        for (std::size_t j=0;j+1<G.ny;++j) {
            int v00 = vid(i,j),   v01 = vid(i,j+1);
            int v10 = vid(i+1,j), v11 = vid(i+1,j+1);
            std::fprintf(f, "f %d %d %d\n", v00, v01, v11);
            std::fprintf(f, "f %d %d %d\n", v00, v11, v10);
        }
    }
    std::fclose(f);
    TRACE_INFO("save_obj", path);
    return true;
}

} // namespace rslm::diag
