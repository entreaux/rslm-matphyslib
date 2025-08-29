#pragma once
/**
 * Trace helpers (TXT logger)
 * - TRACE_SCOPE("name") emits scope_enter/scope_exit at Level::trace
 * - TRACE_VAR(x)        emits x at Level::trace
 * - TRACE_INFO/DEBUG/WARN/ERROR  key/value helpers at fixed levels
 */

#include "logger.hpp"
#include "format.hpp"

#define RSLM_LOC ::rslm::telemetry::Location{__FILE__, __LINE__, __func__}

namespace rslm::telemetry {

class TraceScope {
public:
    explicit TraceScope(const char* scope_name) : scope_(scope_name) {
        Logger::instance().log(Level::trace, RSLM_LOC, scope_, "scope_enter", scope_);
    }
    ~TraceScope() {
        Logger::instance().log(Level::trace, RSLM_LOC, scope_, "scope_exit", scope_);
    }
private:
    std::string scope_;
};

} // namespace rslm::telemetry

#define TRACE_SCOPE(name_literal) ::rslm::telemetry::TraceScope _rslm_scope_guard_(name_literal)

#define TRACE_VAR(var_expr) \
    do { \
        auto& _L = ::rslm::telemetry::Logger::instance(); \
        if (_L.enabled()) { \
            _L.log(::rslm::telemetry::Level::trace, RSLM_LOC, "", #var_expr, ::rslm::telemetry::format::to_string_generic((var_expr))); \
        } \
    } while (0)

#define TRACE_INFO(key_str, value_expr) \
    do { \
        auto& _L = ::rslm::telemetry::Logger::instance(); \
        if (_L.enabled()) { \
            _L.log(::rslm::telemetry::Level::info, RSLM_LOC, "", (key_str), ::rslm::telemetry::format::to_string_generic((value_expr))); \
        } \
    } while (0)

#define TRACE_DEBUG(key_str, value_expr) \
    do { \
        auto& _L = ::rslm::telemetry::Logger::instance(); \
        if (_L.enabled()) { \
            _L.log(::rslm::telemetry::Level::debug, RSLM_LOC, "", (key_str), ::rslm::telemetry::format::to_string_generic((value_expr))); \
        } \
    } while (0)

#define TRACE_WARN(key_str, value_expr) \
    do { \
        auto& _L = ::rslm::telemetry::Logger::instance(); \
        if (_L.enabled()) { \
            _L.log(::rslm::telemetry::Level::warn, RSLM_LOC, "", (key_str), ::rslm::telemetry::format::to_string_generic((value_expr))); \
        } \
    } while (0)

#define TRACE_ERROR(key_str, value_expr) \
    do { \
        auto& _L = ::rslm::telemetry::Logger::instance(); \
        if (_L.enabled()) { \
            _L.log(::rslm::telemetry::Level::error_, RSLM_LOC, "", (key_str), ::rslm::telemetry::format::to_string_generic((value_expr))); \
        } \
    } while (0)
