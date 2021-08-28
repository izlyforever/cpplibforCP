#pragma once
#include <bits/stdc++.h>

// assume M is NTT-friendly and 3 is a primitive root of N
template<int M>
class NTT {
  std::vector<int> rev;
  std::vector<MInt<M>> roots{0, 1};
 public:
  static inline const MInt<M> g = 3;
  void dft(std::vector<MInt<M>> &a) {
    int n = a.size();
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
        auto e = pow(g, (M - 1) >> (k + 1));
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
  void idft(std::vector<MInt<M>> &a) {
    int n = a.size();
    std::reverse(a.begin() + 1, a.end());
    dft(a);
    auto inv = pow(MInt<M>(n), M - 2);
    for (auto &x : a) x *= inv;
  }
};
