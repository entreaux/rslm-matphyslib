#pragma once
/**
 * RSLM Telemetry: Logger (header) — TXT output
 * --------------------------------------------
 * Thread-safe, line-oriented TXT logger with:
 *  - 128-bit RUN/session ID (UUID-like string)
 *  - Monotonic 64-bit log_id across the process
 *  - Console mirroring (human readable)
 *  - File logs at ./logs/<RUN>.txt
 *  - Levels: TRACE, DEBUG, INFO, WARN, ERROR
 *  - Location metadata: file, line, function
 *
 * Each line is tab-separated key=value pairs:
 *   ts=...    run_id=...   log_id=42   level=TRACE   file=...:line  func=...
 *   scope=... name=...     value_ck=...   value="escaped string"
 */

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>

namespace rslm::telemetry {

enum class Level : uint8_t { trace=0, debug=1, info=2, warn=3, error_=4 };

struct Location {
    std::string file;
    int line = 0;
    std::string func;
};

struct Record {
    Level level = Level::info;
    std::string ts_iso8601;
    std::string run_id;
    uint64_t log_id = 0;
    Location loc;
    std::string scope;
    std::string name;
    std::string value_str;
    uint64_t value_ck = 0;
};

class Logger {
public:
    static Logger& instance();

    void start_run(std::optional<std::string> forced_run_id = std::nullopt);
    void stop_run();

    void set_console_mirror(bool enabled);
    void set_min_level(Level lvl);

    bool enabled() const noexcept { return enabled_; }
    void set_enabled(bool e) { enabled_ = e; }

    void log(Level lvl, const Location& loc,
             std::string scope, std::string name, std::string value_str);

    std::string run_id() const { return run_id_; }
    uint64_t next_log_id() { return ++log_counter_; }
    
    void set_flush_every(std::uint64_t n) { std::scoped_lock lk(mtx_); flush_every_ = (n == 0 ? 1 : n); }


private:
    Logger();
    ~Logger();

    std::string gen_run_id_() const; // timestamp + 128-bit rand
    std::string now_iso8601_() const;
    static uint64_t fnv1a64_(std::string_view s);
    void ensure_log_dir_open_();
    void write_text_line_(const Record& rec);
    void write_console_line_(const Record& rec);
    
    std::uint64_t flush_every_ = 1;

private:
    std::atomic<uint64_t> log_counter_{0};
    std::string run_id_;
    std::ofstream file_;
    std::filesystem::path file_path_;
    std::mutex mtx_;
    bool console_mirror_ = true;
    bool enabled_ = true;
    Level min_level_ = Level::trace;
    bool started_ = false;
};

} // namespace rslm::telemetry
