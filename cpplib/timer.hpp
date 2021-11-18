#pragma once
#include <bits/stdc++.h>

class Timer {
  std::chrono::high_resolution_clock::time_point start;
 public:
  Timer() : start(std::chrono::high_resolution_clock::now()) {}
  void show(std::string s = {}) {
    auto elapsedTime = std::chrono::high_resolution_clock::now() - start;
    std::cerr << "[Time used: " << s << "] " << elapsedTime.count() / 1'000'000.0 << "ms\n";
  }
};
