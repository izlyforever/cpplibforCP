#pragma once
#include <bits/stdc++.h>

class Timer {
  std::chrono::steady_clock::time_point start;
 public:
  Timer() : start(std::chrono::steady_clock::now()) {}
  void show() {
    auto elapsedTime = std::chrono::steady_clock::now() - start;
    std::cout << "Time used: " << elapsedTime.count() / 1000000 << "ms\n";
  }
};

template<typename T>
void debug(std::vector<T> a){
  for (auto &i : a) std::cout << i << ' ';
  std::cout << std::endl;
}