#pragma once
#include <bits/stdc++.h>

namespace FMT {
constexpr int M = 998244353, inv2 = (M + 1) / 2;
auto add = [](int& x, int y) {
  (x += y) >= M && (x -= M);
};
auto sub = [](int& x, int y) {
  (x -= y) < 0 && (x += M);
};
auto extend = [](int n) {
  int r = std::__lg(n);
  while ((1 << r) < n) ++r;
  return r;
};
auto FMTor = [](std::vector<int>& a, bool isRev) {
  int n = extend(a.size());
  a.resize(1 << n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < (1 << n); ++j) if ((j >> i) & 1) {
      if (isRev) sub(a[j], a[j ^ (1 << i)]);
      else add(a[j], a[j ^ (1 << i)]);
    }
  }
};
auto FMTand = [](std::vector<int>& a, bool isRev) {
  int n = extend(a.size());
  a.resize(1 << n);
  for (int i = 0; i < n; ++i) {
    for (int j = (1 << n) - 1; j >= 0; --j) if ((j >> i) & 1) {
      if (isRev) sub(a[j ^ (1 << i)], a[j]);
      else add(a[j ^ (1 << i)], a[j]);
    }
  }
};
auto FMTxor = [](std::vector<int>& a, bool isRev) {
  int n = extend(a.size());
  a.resize(1 << n);
  for (int i = 0; i < n; ++i) {
    for (int j = (1 << n) - 1; j >= 0; --j) if ((j >> i) & 1) {
      int u = a[j], v = a[j ^ (1 << i)];
      a[j] = (v - u + M) % M;
      a[j ^ (1 << i)] = (u + v) % M;
    }
    if (isRev) for (auto& x : a) x = 1LL * inv2 * x % M;
  }
};
auto fun = [](std::function<void(std::vector<int>& , bool)> f, std::vector<int> a, std::vector<int> b) {
  int n = extend(std::max(a.size(), b.size()));
  a.resize(1 << n); b.resize(1 << n);
  f(a, 0); f(b, 0);
  std::vector<int> c(1 << n);
  for (int i = 0; i < (1 << n); ++i) c[i] = 1LL * a[i] * b[i] % M;
  f(c, 1);
  return c;
};
auto Or = [](std::vector<int> a, std::vector<int> b) {
  return fun(FMTor, a, b);
};
auto And = [](std::vector<int> a, std::vector<int> b) {
  return fun(FMTand, a, b);
};
auto Xor = [](std::vector<int> a, std::vector<int> b) {
  return fun(FMTxor, a, b);
};
// i = j | k and j & k = 0
auto OrAnd = [](std::vector<int> a, std::vector<int> b) {
  int n = extend(std::max(a.size(), b.size()));
  a.resize(1 << n); b.resize(1 << n);
  std::vector<std::vector<int>> sa(n + 1, std::vector<int>(1 << n));
  auto sb = sa, sc = sa;
  for (int i = 0; i < (1 << n); ++i) sa[__builtin_popcount(i)][i] = a[i];
  for (int i = 0; i < (1 << n); ++i) sb[__builtin_popcount(i)][i] = b[i];
  for (int i = 0; i <= n; ++i) {
    FMTor(sa[i], 0);FMTor(sb[i], 0);
  }
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= i; ++j) {
      for (int k = 0; k < (1 << n); ++k) {
        add(sc[i][k], 1LL * sa[j][k] * sb[i - j][k] % M);
      }
    }
    FMTor(sc[i], 1);
  }
  std::vector<int> c(1 << n);
  for (int i = 0; i < (1 << n); ++i) c[i] = sc[__builtin_popcount(i)][i];
  return c;
};
} // namespace FMT
// https://www.luogu.com.cn/problem/P6097
