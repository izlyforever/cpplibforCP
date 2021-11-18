#include <bits/stdc++.h>
#include "../../cpplib/math/mixed.hpp"

int main() {
  std::cin.tie(nullptr)->sync_with_stdio(false);

  std::vector<int> a{1, 2, 5, 6, 8, 17};
  KnuthShuffle(a);
  for (auto &x : a) std::cout << x << ' ';
  std::cout << '\n';

  auto b = uniformChoose(1e9, 10);
  for (auto &x : b) std::cout << x << ' ';
  std::cout << '\n'; 
  
  return 0;
}