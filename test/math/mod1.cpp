#include <bits/stdc++.h>
#include "../../cpplib/math/mod.hpp"

constexpr int M = 998244353;
using mod = MInt<M>;

int main() {
  //freopen("in", "r", stdin);
  std::cin.tie(nullptr)->sync_with_stdio(false);
  int a, b;
  std::cin >> a >> b;

  mod a1(a), b1(b);
  std::cout << a1 + b1 << '\n';

  ModInt::setMod(M);
  ModInt a2(a), b2(b);
  std::cout << a2 - b2 << '\n';

  ModLL::setMod(M);
  ModInt a3(a), b3(b);
  std::cout << a3 * b3 << '\n';
  std::cout << a3 / b3 << '\n';
  return 0;
}