// docs/test/math/poly1.cpp
#include <bits/stdc++.h>
using LL = long long;

#include "../../cpplib/math.hpp"

// https://www.spoj.com/problems/FACTMODP/en/
int main() {
  // freopen("in", "r", stdin);
  std::cin.tie(nullptr)->sync_with_stdio(false);
  int cas = 1;
  std::cin >> cas;
  while (cas--) {
    LL n, p;
    std::cin >> n >> p;
    std::cout << factorial(n, p) << '\n';
    // std::cout << factorialOrigin(n, p) << '\n';
  }
  return 0;
}