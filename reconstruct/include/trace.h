#pragma once

#include <chrono>
#include <iostream>
#include <mutex>
#include <sstream>

/* ============================================================
 *  Compile-time switches
 * ============================================================ */

// default: event on
#ifndef TRACE_EVENT_ENABLED
#define TRACE_EVENT_ENABLED 1
#endif

// default: duration on
#ifndef TRACE_DURATION_ENABLED
#define TRACE_DURATION_ENABLED 1
#endif

/* ============================================================
 *  Event trace
 * ============================================================ */

#if TRACE_EVENT_ENABLED

// 线程安全的输出锁
inline std::mutex& get_trace_mutex() {
    static std::mutex trace_mutex;
    return trace_mutex;
}

inline void trace_event(const char* name, int depth) {
    std::lock_guard<std::mutex> lock(get_trace_mutex());
    
    for (int i = 0; i < depth; ++i)
        std::cout << "  ";

    std::cout << name << std::endl;
}

#define TRACE_EVENT(name, depth) trace_event(name, depth)

#else

#define TRACE_EVENT(name, depth) do {} while (0)

#endif

/* ============================================================
 *  Duration trace with nesting (RAII)
 * ============================================================ */

#if TRACE_DURATION_ENABLED

inline thread_local int g_trace_depth = 0;

class TraceScope {
public:
    explicit TraceScope(const char* name)
        : name_(name),
          depth_(g_trace_depth++),
          start_(std::chrono::high_resolution_clock::now()) {}

    ~TraceScope() {
        auto end = std::chrono::high_resolution_clock::now();
        auto us =
            std::chrono::duration_cast<std::chrono::microseconds>(
                end - start_
            ).count();

        g_trace_depth--;

        // 使用stringstream构建完整输出，然后一次性输出确保原子性
        std::ostringstream oss;
        for (int i = 0; i < depth_; ++i)
            oss << "  ";
        oss << name_ << ": " << us << " us";

        {
            std::lock_guard<std::mutex> lock(get_trace_mutex());
            std::cout << oss.str() << std::endl;
        }
    }

    int depth() const { return depth_; }

private:
    const char* name_;
    int depth_;
    std::chrono::high_resolution_clock::time_point start_;
};

#define TRACE_SCOPE(name) \
    [[maybe_unused]] TraceScope _trace_scope_##__LINE__(name)

#else

#define TRACE_SCOPE(name) do {} while (0)

#endif

/* ============================================================
 *  Unified trace point
 * ============================================================ */

#define TRACE(name)            \
    TRACE_SCOPE(name);         \
    TRACE_EVENT(name,          \
        TRACE_DURATION_ENABLED ? g_trace_depth - 1 : 0)
