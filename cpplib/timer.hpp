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
  static std::string nowTime() {
    time_t t = time(0);
    char buffer[9] = {0};
    strftime(buffer, 9, "%H:%M:%S", localtime(&t));
    return std::string(buffer);
  }
};