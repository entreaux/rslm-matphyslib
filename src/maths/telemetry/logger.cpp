#include "logger.hpp"
#include "format.hpp"

#include <cstdio>
#include <cstdlib>
#include <ctime>

namespace rslm::telemetry {

static const char* level_name(Level L) {
    switch (L) {
        case Level::trace: return "trace";
        case Level::debug: return "debug";
        case Level::info:  return "info";
        case Level::warn:  return "warn";
        case Level::error_:return "error";
    }
    return "TRACE";
}

Logger& Logger::instance() {
    static Logger inst;
    return inst;
}

Logger::Logger() {
    if (const char* e = std::getenv("RSLM_TRACE")) {
        enabled_ = (std::string(e) != "0");
    }
}

Logger::~Logger() {
    stop_run();
}

void Logger::start_run(std::optional<std::string> forced_run_id) {
    std::scoped_lock lk(mtx_);
    if (started_) return;

    run_id_ = forced_run_id.value_or(gen_run_id_());
    ensure_log_dir_open_();
    started_ = true;

    // Write header directly to avoid re-entrant lock (no call to log()).
    Record rec;
    rec.level      = Level::info;
    rec.ts_iso8601 = now_iso8601_();
    rec.run_id     = run_id_;
    rec.log_id     = next_log_id();
    rec.loc        = Location{__FILE__, __LINE__, __func__};
    rec.scope      = "logger";
    rec.name       = "run_start";
    rec.value_str  = "run_id=" + run_id_;
    rec.value_ck   = fnv1a64_(rec.value_str);

    write_text_line_(rec);
    if (console_mirror_) write_console_line_(rec);
}

void Logger::stop_run() {
    std::scoped_lock lk(mtx_);
    if (!started_) return;

    Record rec;
    rec.level = Level::info;
    rec.ts_iso8601 = now_iso8601_();
    rec.run_id = run_id_;
    rec.log_id = next_log_id();
    rec.loc = Location{__FILE__, __LINE__, __func__};
    rec.scope = "logger";
    rec.name = "run_stop";
    rec.value_str = "closing";
    rec.value_ck = fnv1a64_(rec.value_str);

    write_text_line_(rec);
    if (console_mirror_) write_console_line_(rec);

    if (file_.is_open()) file_.close();
    started_ = false;
}

void Logger::set_console_mirror(bool enabled) {
    std::scoped_lock lk(mtx_);
    console_mirror_ = enabled;
}

void Logger::set_min_level(Level lvl) {
    std::scoped_lock lk(mtx_);
    min_level_ = lvl;
}

void Logger::log(Level lvl, const Location& loc,
                 std::string scope, std::string name, std::string value_str) {
    if (!enabled_) return;
    if (static_cast<int>(lvl) < static_cast<int>(min_level_)) return;

    std::scoped_lock lk(mtx_);
    if (!started_) {
        run_id_ = gen_run_id_();
        ensure_log_dir_open_();
        started_ = true;
    }

    Record rec;
    rec.level = lvl;
    rec.ts_iso8601 = now_iso8601_();
    rec.run_id = run_id_;
    rec.log_id = next_log_id();
    rec.loc = loc;
    rec.scope = std::move(scope);
    rec.name = std::move(name);
    rec.value_str = std::move(value_str);
    rec.value_ck = fnv1a64_(rec.value_str);

    write_text_line_(rec);
    if (console_mirror_) write_console_line_(rec);
}

std::string Logger::gen_run_id_() const {
    using namespace std;
    using namespace std::chrono;

    auto now = system_clock::now();
    auto secs = time_point_cast<std::chrono::seconds>(now);
    std::time_t t = system_clock::to_time_t(secs);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char ts[32];
    std::strftime(ts, sizeof(ts), "%Y-%m-%dT%H-%M-%S", &tm);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> U;
    uint64_t a = U(gen), b = U(gen);

    std::ostringstream oss;
    oss << ts << "_"
        << std::hex << std::setfill('0')
        << std::setw(16) << a << "-"
        << std::setw(16) << b;
    return oss.str();
}

std::string Logger::now_iso8601_() const {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto t = system_clock::to_time_t(now);
    auto us = duration_cast<microseconds>(now.time_since_epoch()) % 1'000'000;

    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[64];
    std::snprintf(buf, sizeof(buf),
                  "%04d-%02d-%02dT%02d:%02d:%02d.%06lldZ",
                  tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
                  tm.tm_hour, tm.tm_min, tm.tm_sec,
                  static_cast<long long>(us.count()));
    return std::string(buf);
}

uint64_t Logger::fnv1a64_(std::string_view s) {
    const uint64_t FNV_OFFSET = 14695981039346656037ull;
    const uint64_t FNV_PRIME  = 1099511628211ull;
    uint64_t h = FNV_OFFSET;
    for (unsigned char c : s) { h ^= c; h *= FNV_PRIME; }
    return h;
}

void Logger::ensure_log_dir_open_() {
    std::filesystem::path dir = std::filesystem::path("logs");
    std::error_code ec;
    if (!std::filesystem::exists(dir))
        std::filesystem::create_directories(dir, ec);

    file_path_ = dir / (run_id_ + ".txt");
    file_.open(file_path_, std::ios::out | std::ios::app);
}

void Logger::write_text_line_(const Record& r) {
    using format::escape_txt;

    // key=value pairs, tab-separated. value is always quoted and escaped.
    file_
        << "ts="     << r.ts_iso8601
        << "\trun_id=" << r.run_id
        << "\tlog_id=" << r.log_id
        << "\tlevel="  << level_name(r.level)
        << "\tfile="   << escape_txt(r.loc.file) << ":" << r.loc.line
        << "\tfunc="   << escape_txt(r.loc.func)
        << "\tscope="  << escape_txt(r.scope)
        << "\tname="   << escape_txt(r.name)
        << "\tvalue_ck=" << r.value_ck
        << "\tvalue=\""  << escape_txt(r.value_str) << "\""
        << "\n";
    if (flush_every_ == 1 || (r.log_id % flush_every_) == 0) {
        file_.flush();
    }
}

void Logger::write_console_line_(const Record& r) {
    std::fprintf(stdout,
        "[%s] (%s) #%llu %s %s:%d %s [%s] %s = \"%s\"\n",
        r.ts_iso8601.c_str(),
        r.run_id.c_str(),
        static_cast<unsigned long long>(r.log_id),
        level_name(r.level),
        r.loc.file.c_str(), r.loc.line,
        r.loc.func.c_str(),
        r.scope.c_str(),
        r.name.c_str(),
        r.value_str.c_str()
    );
}

} // namespace rslm::telemetry
