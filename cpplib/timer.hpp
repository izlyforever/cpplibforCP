#pragma once
#include <bits/stdc++.h>

class Timer {
  std::chrono::steady_clock::time_point start_;
 public:
  Timer() : start_(std::chrono::steady_clock::now()) {}
  void show() {
    auto elapsedTime = std::chrono::steady_clock::now() - start_;
    std::cerr << "Time used: " << elapsedTime.count() / 1'000'000 << "ms\n";
  }
};
