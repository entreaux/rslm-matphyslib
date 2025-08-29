#pragma once
/**
 * RSLM Telemetry: Formatting helpers
 * ----------------------------------
 * - Plain-text escaping (for .txt logs)
 * - Scalar formatting with fixed precision
 * - Generic to_string fallback
 */

#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <cstdio>

namespace rslm::telemetry::format {

// Escape for our .txt key=value lines.
// - Replace tabs/newlines with visible escapes.
// - Quote " when used inside a quoted value.
inline std::string escape_txt(std::string_view s) {
    std::string out;
    out.reserve(s.size() + 16);
    for (char c : s) {
        switch (c) {
            case '\t': out += "\\t"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\"': out += "\\\""; break;
            default:   out += c; break;
        }
    }
    return out;
}

template <typename T>
inline std::string to_string_num(T v, int decimals = 6) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss << std::setprecision(decimals) << v;
    return oss.str();
}

template <typename T>
inline std::string to_string_generic(const T& v) {
    std::ostringstream oss;
    oss << v;
    return oss.str();
}

} // namespace rslm::telemetry::format
