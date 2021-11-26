#pragma once
#include <bits/stdc++.h>

class Timer final {
  std::chrono::high_resolution_clock::time_point start_;
  std::string name_;
 public:
  Timer(std::string name = {}) : start_(std::chrono::high_resolution_clock::now()), name_(name) {}
  ~Timer() {
    auto elapsedTime = std::chrono::high_resolution_clock::now() - start_;
    std::cerr << "[Time used: " << name_ << "] " << elapsedTime.count() / 1'000'000.0 << "ms\n";
  }
};