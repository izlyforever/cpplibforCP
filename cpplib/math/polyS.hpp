#pragma once
#include <bits/stdc++.h>
using LL = long long;

class PolyS : public std::vector<int> {
  static inline std::vector<int> rev_, roots_{0, 1};
  using ULL = unsigned long long;
  static unsigned powMod(unsigned x, unsigned n) {
    static const unsigned m = 998244353U;
    static const unsigned mr = 998244351U;
    static const unsigned m1 = 301989884U;
    static const unsigned m1inv = 232013824U;
    unsigned xx = (ULL(x) << 32) % m, rr = m1;
    while (n) {
      if (n & 1) {
        ULL t = ULL(rr) * xx;
        rr = (t + ULL(unsigned(t) * mr) * m) >> 32;
      }
      ULL t = ULL(xx) * xx;
      xx = (t + ULL(unsigned(t) * mr) * m) >> 32;
      n >>= 1;
    }
    return ULL(rr) * m1inv % m;
  }
  void dft() {
    int n = (int)size();
    if ((int)rev_.size() != n) {
      int k = __builtin_ctz(n) - 1;
      rev_.resize(n);
      for (int i = 0; i < n; ++i) {
        rev_[i] = rev_[i >> 1] >> 1 | (i & 1) << k;
      }
    }
    if ((int)roots_.size() < n) {
      int k = __builtin_ctz(roots_.size());
      roots_.resize(n);
      while ((1 << k) < n) {
        int e = powMod(G, (M - 1) >> (k + 1));
        for (int i = 1 << (k - 1); i < (1 << k); ++i) {
          roots_[2 * i] = roots_[i];
          roots_[2 * i + 1] = 1LL * roots_[i] * e % M;
        }
        ++k;
      }
    }
    for (int i = 0; i < n; ++i) if (rev_[i] < i) {
      std::swap((*this)[i], (*this)[rev_[i]]);
    }
    for (int k = 1; k < n; k *= 2) {
      for (int i = 0; i < n; i += 2 * k) {
        for (int j = 0; j < k; ++j) {
          int u = (*this)[i + j];
          int v = 1LL * (*this)[i + j + k] * roots_[k + j] % M;
          int x = u + v, y = u - v;
          if (x >= M) x -= M;
          if (y < 0) y += M;
          (*this)[i + j] = x;
          (*this)[i + j + k] = y;
        }
      }
    }
  }
  void idft() {
    int n = (int)size();
    std::reverse(begin() + 1, end());
    dft();
    // not that n is power of 2, and M = 1 + c 2^x 
    int invN = M - (M - 1) / n;
    for (int i = 0; i < n; ++i) {
      (*this)[i] = 1LL * (*this)[i] * invN % M;
    }
  }
  void standard() {
    while (!empty() && !back()) pop_back();
  }
  void reverse() {
    std::reverse(begin(), end());
    standard();
  }
 public:
  static inline constexpr int G = 3, M = 998244353; // 1 +  2^23 * 7 * 17
  PolyS() {}
  PolyS(const int& x) : std::vector<int>{x} { standard();}
  PolyS(const std::vector<int>& a) : std::vector<int>{a} { standard();}
  PolyS(std::vector<int>&& a) : std::vector<int>(std::move(a)) { standard();}
  int at(int id) const {
    if (id < 0 || id >= (int)size()) return 0;
    return (*this)[id];
  }
  PolyS operator-() const {
    auto A = (*this);
    for (auto& x : A) x = (x == 0 ? 0 : M - x);
    return A;
  }
  PolyS mulXn(int n) const {
    auto A = *this;
    if (!A.empty()) A.insert(A.begin(), n, 0);
    return A;
  }
  PolyS modXn(int n) const {
    if (n >= (int)size()) return *this;
    return PolyS({begin(), begin() + n});
  }
  PolyS modXnR(int n) {
    this->resize(n);
    return PolyS(std::move(*this));
  }
  PolyS divXn(int n) const {
    if ((int)size() <= n) return PolyS();
    return PolyS({begin() + n, end()});
  }
  PolyS& operator+=(const PolyS& rhs) {
    if ((int)size() < (int)rhs.size()) resize(rhs.size());
    for (int i = 0, rs = rhs.size(); i < rs; ++i) {
      if (((*this)[i] += rhs[i]) >= M) (*this)[i] -= M;
    }
    standard();
    return *this;
  }
  PolyS& operator-=(const PolyS& rhs) {
    if (size() < rhs.size()) resize(rhs.size());
    for (int i = 0, rs = rhs.size(); i < rs; ++i) {
      if (((*this)[i] -= rhs[i]) < 0) (*this)[i] += M;
    }
    standard();
    return *this;
  }
  PolyS& operator*=(PolyS&& rhs) {
    int n = (int)size(), m = rhs.size(), tot = std::max(1, n + m - 1);
    int sz = 1 << std::__lg(tot * 2 - 1);
    resize(sz);
    rhs.resize(sz);
    dft();
    rhs.dft();
    for (int i = 0; i < sz; ++i) {
      (*this)[i] = 1LL * (*this)[i] * rhs[i] % M;
    }
    idft();
    standard();
    return *this;
  }
  PolyS& operator/=(PolyS&& rhs) {
    int n = (int)size(), m = rhs.size();
    if (n < m) return (*this) = PolyS();
    reverse();
    rhs.reverse();
    (*this) *= rhs.inv(n - m + 1);
    resize(n - m + 1);
    reverse();
    return *this;
  }
  PolyS& operator*=(const PolyS& rhs) {
    return (*this) *= PolyS(rhs);
  }
  PolyS& operator/=(const PolyS& rhs) {
    return (*this) /= PolyS(rhs);
  }
  PolyS& operator%=(const PolyS& rhs) {
    return (*this) -= (*this) / rhs * rhs;
  }
  PolyS operator+(const PolyS& rhs) const {
    return PolyS(*this) += rhs;
  }
  PolyS operator-(const PolyS& rhs) const {
    return PolyS(*this) -= rhs;
  }
  PolyS operator*(const PolyS& rhs) const {
    return PolyS(*this) *= rhs;
  }
  PolyS operator/(const PolyS& rhs) const {
    return PolyS(*this) /= rhs;
  }
  PolyS operator%(const PolyS& rhs) const {
    return PolyS(*this) %= rhs;
  }
  PolyS powModPoly(int n, PolyS p) const {
    PolyS r(1), x(*this);
    while (n) {
      if (n&1) (r *= x) %= p;
      n >>= 1; (x *= x) %= p;
    }
    return r;
  }
  int inner(const PolyS& rhs) const {
    int r = 0, n = std::min(size(), rhs.size());
    for (int i = 0; i < n; ++i) {
      r = (r + 1LL * (*this)[i] * rhs[i]) % M;
    }
    return r;
  }
  PolyS derivation() const {
    if (empty()) return PolyS();
    int n = (int)size();
    std::vector<int> r(n - 1);
    for (int i = 1; i < n; ++i) r[i - 1] =  1LL * (*this)[i] * i % M;
    return PolyS(r);
  }
  PolyS integral() const {
    if (empty()) return PolyS();
    int n = (int)size();
    std::vector<int> r(n + 1), inv(n + 1, 1);
    for (int i = 2; i <= n; ++i) inv[i] = 1LL * (M - M / i) * inv[M % i] % M;
    for (int i = 0; i < n; ++i) r[i + 1] = 1LL * (*this)[i] * inv[i + 1] % M;
    return PolyS(r);
  }
  PolyS inv(int n) const {
    // assert((*this)[0] != 0);
    PolyS x(powMod((*this)[0], M - 2));
    int k = 1;
    while (k < n) {
      k *= 2;
      x *= (PolyS(2) - this->modXn(k) * x).modXnR(k);
    }
    return x.modXnR(n);
  }
  // assume a[0] = 1
  PolyS log(int n) const {
    return (derivation() * inv(n)).integral().modXnR(n);
  }
  // assume a[0] = 0
  PolyS exp(int n) const {
    PolyS x(1);
    int k = 1;
    while (k < n) {
      k *= 2;
      x = (x * (PolyS(1) - x.log(k) + this->modXn(k))).modXnR(k);
    }
    return x.modXnR(n);
  }
  // assume a[0] = 1;
  PolyS sqrt(int n) const {
    const static int inv2 = (M + 1) / 2;
    PolyS x(1);
    int k = 1;
    while (k < n) {
      k *= 2;
      x += this->modXn(k) * x.inv(k);
      x = x.modXnR(k) * inv2;
    }
    return x.modXnR(n);
  }
  int eval(int x) const {
    int r = 0, t = 1;
    for (int i = 0, n = (int)size(); i < n; ++i) {
      r = (r + 1LL * (*this)[i] * t) % M;
      t = 1LL * t * x % M;
    }
    return r;
  }
  // transpose convolution {\rm MULT}(F(x),G(x))=\sum_{i\ge0}(\sum_{j\ge 0}f_{i+j}g_j)x^i
  PolyS mulT(PolyS&& rhs) const {
    if (rhs.size() == 0) return PolyS();
    int n = rhs.size();
    std::reverse(rhs.begin(), rhs.end());
    return ((*this) * rhs).divXn(n - 1);
  }
  // multi-evaluation(new tech)
  std::vector<int> evals(const std::vector<int>& x) const {
    if (size() == 0) return std::vector<int>(x.size());
    int n = (int)x.size();
    std::vector<int> ans(n);
    std::vector<PolyS> g(4 * n);
    std::function<void(int, int, int)> build = [&](int l, int r, int p) {
      if (r - l == 1) {
        // g[p] = std::vector<int>{1, x[i] ? M - x[l] : 0};
        g[p] = PolyS({1, x[l] ? M - x[l] : 0});
      } else {
        int m = (l + r) / 2;
        build(l, m, 2 * p);
        build(m, r, 2 * p + 1);
        g[p] = g[2 * p] * g[2 * p + 1];
      }
    };
    build(0, n, 1);
    std::function<void(int, int, int, const PolyS& )> solve = [&](int l, int r, int p, const PolyS& f) {
      if (r - l == 1) {
        ans[l] = f.at(0);
      } else {
        int m = (l + r) / 2;
        solve(l, m, 2 * p, f.mulT(std::move(g[2 * p + 1])).modXnR(m - l));
        solve(m, r, 2 * p + 1, f.mulT(std::move(g[2 * p])).modXnR(r - m));
      }
    };
    solve(0, n, 1, mulT(g[1].inv(size())).modXnR(n));
    return ans;
  } // https://www.luogu.com.cn/problem/P5050
}; // https://www.luogu.com.cn/training/3015#information
