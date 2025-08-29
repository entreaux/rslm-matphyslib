#pragma once
/**
 * RSLM Maths — diagnostics/palette.hpp
 * ------------------------------------
 * Minimal color maps for heatmaps. We keep it tiny and dependency-free.
 * The map() function takes normalized t∈[0,1] and returns uint8 RGB.
 */

#include <cstdint>
#include <algorithm>

namespace rslm::diag {

struct RGB { std::uint8_t r, g, b; };

// Clamp helper
inline double clamp01(double x) { return x < 0 ? 0 : (x > 1 ? 1 : x); }

// Linear interpolate between two RGB
inline RGB lerp(const RGB& a, const RGB& b, double t) {
    auto L = [&](std::uint8_t A, std::uint8_t B) -> std::uint8_t {
        return static_cast<std::uint8_t>(A + (B - A) * t + 0.5);
    };
    return { L(a.r,b.r), L(a.g,b.g), L(a.b,b.b) };
}

/** “Thermal” 5-stop palette: blue → cyan → yellow → orange → red. */
struct Thermal5 {
    static RGB map(double t01) {
        const RGB stops[5] = {
            {  0,  32, 128}, // deep blue
            {  0, 180, 180}, // cyan
            {250, 250,  70}, // yellow
            {240, 150,  40}, // orange
            {200,  30,  30}, // red
        };
        t01 = clamp01(t01);
        const double seg = 1.0 / 4.0;
        int i = static_cast<int>(t01 / seg);
        if (i >= 4) i = 3;
        double lt = (t01 - i*seg) / seg;
        return lerp(stops[i], stops[i+1], lt);
    }
};

/** Grayscale palette (for quick comparisons). */
struct Gray {
    static RGB map(double t01) {
        std::uint8_t v = static_cast<std::uint8_t>(clamp01(t01) * 255.0 + 0.5);
        return {v, v, v};
    }
};

} // namespace rslm::diag
