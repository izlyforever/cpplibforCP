#pragma once
#include <bits/stdc++.h>

class Timer final {
  std::chrono::high_resolution_clock::time_point start_;
  std::string name_;
 public:
  Timer(std::string name = {}) : start_(std::chrono::high_resolution_clock::now()), name_(name) {}
  ~Timer() {
    auto elapsedTime = std::chrono::high_resolution_clock::now() - start_;
    std::cerr << std::setprecision(3) << std::fixed << "[Time used: " <<
        name_ << "] " << elapsedTime.count() / 1'000'000.0 << "ms\n";
  }

  static uint64_t now_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::system_clock::now().time_since_epoch()).count();
  }

  static uint64_t tick_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now().time_since_epoch()).count();
  }

  static std::tm localTime(const std::time_t& time) {
    std::tm tm{};
#ifdef _WIN32
    ::localtime_s(&tm, &time);
#else
    ::localtime_r(&time, &tm);
  #endif
    return tm;
  }

  static std::tm localTime() {
    std::time_t now_t = ::time(nullptr);
    return localTime(now_t);
  }

  // Return local time (like 04:05:06.789), for instance: UTC+8 in China
  static std::string localTimeString() {
    std::tm time = localTime();
    std::stringstream sstream;
    sstream << std::put_time(&time, "%T");
    return sstream.str();
  }

  static std::string localFullTimeString() {
    std::tm time = localTime();
    std::stringstream sstream;
    sstream << std::put_time(&time, "%c %T");
    sstream << '.' << std::setfill('0') << std::setw(3) << now_ms() % 1000;
    return sstream.str();
  }
};
