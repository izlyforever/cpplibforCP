// docs/test/math/gauss1.cpp
#include <bits/stdc++.h>
#define clog(x) std::clog << (#x) << " is " << (x) << '\n';
using LL = long long;
#include "../../cpplib/math/mixed.hpp"

int main() {
  //freopen("in", "r", stdin);
  std::cin.tie(nullptr)->sync_with_stdio(false);

  std::vector<std::vector<double>> A{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  std::vector<double> b{1.0, 2.0, 3.0};
  for (auto x : Gauss(std::move(A), std::move(b))) std::cout << x << '\n';

  return 0;
}