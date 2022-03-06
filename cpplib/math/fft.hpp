#pragma once
#include <bits/stdc++.h>

namespace FFT {
const double PI = std::acos(-1);
using C = std::complex<double>;
std::vector<int> rev;
std::vector<C> roots{C(0, 0), C(1, 0)};
void dft(std::vector<C>& a) {
  int n = (int)a.size();
  if ((int)rev.size() != n) {
    int k = __builtin_ctz(n) - 1;
    rev.resize(n);
    for (int i = 0; i < n; ++i) {
      rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
    }
  }
  if ((int)roots.size() < n) {
    int k = __builtin_ctz(roots.size());
    roots.resize(n);
    while ((1 << k) < n) {
      C e = std::polar(1.0, PI / (1 << k));
      for (int i = 1 << (k - 1); i < (1 << k); ++i) {
        roots[2 * i] = roots[i];
        roots[2 * i + 1] = roots[i] * e;
      }
      ++k;
    }
  }
  for (int i = 0; i < n; ++i) if (rev[i] < i) {
    std::swap(a[i], a[rev[i]]);
  }
  for (int k = 1; k < n; k *= 2) {
    for (int i = 0; i < n; i += 2 * k) {
      for (int j = 0; j < k; ++j) {
        auto u = a[i + j], v = a[i + j + k] * roots[k + j];
        a[i + j] = u + v;
        a[i + j + k] = u - v;
      }
    }
  }
}
void idft(std::vector<C>& a) {
  int n = (int)a.size();
  std::reverse(a.begin() + 1, a.end());
  dft(a);
  for (auto& x : a) x /= n;
}
} // namespace FFT
