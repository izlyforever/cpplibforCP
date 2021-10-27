#pragma once
#include <bits/stdc++.h>

class Timer {
  std::chrono::steady_clock::time_point start;
 public:
  Timer() : start(std::chrono::steady_clock::now()) {}
  void show(std::string s = {}) {
    auto elapsedTime = std::chrono::steady_clock::now() - start;
    std::cerr << "Time used[" << s << "]: " << elapsedTime.count() / 1000000 << "ms\n";
  }
};