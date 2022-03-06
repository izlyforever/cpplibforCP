// Main reference: https://www.luogu.com.cn/blog/command-block/sheng-cheng-han-shuo-za-tan
#pragma once
#include <bits/stdc++.h>
#include "basic.hpp"
#include "mod.hpp"
using LL = long long;

// many function will fail for the case n > mod
// using valT = decltype(T::a)::value_type;template<typename T, typename valT>
template<typename T, typename valT, typename enable = ModT<valT>>
class Poly : public T {
  static inline const valT INV2 = (valT::mod() + 1) / 2;
  static inline const int MAXN = 1e6 + 2;  // assume size(a) < MAXN
  static inline const auto& BINOM = BinomModp<valT>::Instance(MAXN);
 public:
  using T::T;
  // never use it if valT = MINT<M>, this method is awesome
  static void setMod(LL p, int n = MAXN) {
    valT::setMod(p);
    BinomModp<valT>::Instance(std::min(LL(n), p));
  }
  Poly operator-() const {
    auto A = *this;
    for (auto& x : A) x = -x;
    return A;
  }
  Poly mulXn(int n) const {
    auto A = (*this);
    A.insert(A.begin(), n, 0);
    return A;
  }
  Poly modXn(int n) const {
    if (n >= (int)this->size()) return (*this);
    return Poly({this->begin(), this->begin() + n});
  }
  Poly modXnR(int n) {
    this->resize(n);
    return Poly(std::forward<Poly>(*this));
  }
  Poly divXn(int n) const {
    if ((int)this->size() <= n) return Poly();
    return Poly({this->begin() + n, this->end()});
  }
  Poly& operator+=(const Poly& rhs) {
    if (this->size() < rhs.size()) this->resize(rhs.size());
    for (int i = 0, nr = rhs.size(); i < nr; ++i) (*this)[i] += rhs[i];
    this->standard();
    return *this;
  }
  Poly& operator-=(const Poly& rhs) {
    if (this->size() < rhs.size()) this->resize(rhs.size());
    for (int i = 0, nr = rhs.size(); i < nr; ++i) (*this)[i] -= rhs[i];
    this->standard();
    return *this;
  }
  Poly operator+(const Poly& rhs) const {
    return Poly(*this) += rhs;
  }
  Poly operator-(const Poly& rhs) const {
    return Poly(*this) -= rhs;
  }
  Poly operator*(const Poly& rhs) const {
    return this->mul(rhs);
  }
  Poly& operator*=(const Poly& rhs) {
    return (*this) = (*this) * rhs;
  }
  // assume a[0] \neq 0
  Poly inv(int n) const {
    // assert((*this)[0] != 0);
    Poly x((*this)[0].inv());
    int k = 1;
    while (k < n) {
      k *= 2;
      x *= (Poly(2) - this->modXn(k) * x).modXnR(k);
    }
    return x.modXnR(n);
  }
  Poly& operator/=(Poly rhs) {
    int n = (int)this->size(), m = rhs.size();
    if (n < m) return (*this) = Poly();
    this->reverse();
    rhs.reverse();
    (*this) *= rhs.inv(n - m + 1);
    this->resize(n - m + 1);
    this->reverse();
    return *this;
  }
  Poly operator/(Poly rhs) const {
    return Poly(*this) /= rhs;
  }
  Poly& operator%=(const Poly& rhs) {
    return *this -= (*this) / rhs * rhs;
  }
  Poly operator%(const Poly& rhs) const {
    return Poly(*this) %= rhs;
  }
  Poly powModPoly(LL n, const Poly& p) const {
    Poly r(1), x(*this);
    while (n) {
      if (n&1) r = r * x % p;
      n >>= 1; x = x * x % p;
    }
    return r;
  }
  valT inner(const Poly& rhs) const {
    valT r(0);
    int n = (int)std::min(this->size(), rhs.size());
    for (int i = 0; i < n; ++i) r += (*this)[i] * rhs[i];
    return r;
  }
  Poly derivation() const {
    if (this->empty()) return Poly();
    int n = (int)this->size();
    std::vector<valT> r(n - 1);
    for (int i = 1; i < n; ++i) r[i - 1] = (*this)[i] * valT(i);
    return Poly(std::forward<Poly>(r));
  }
  Poly integral() const {
    if (this->empty()) return Poly();
    int n = (int)this->size();
    std::vector<valT> r(n + 1);
    for (int i = 0; i < n; ++i) r[i + 1] = (*this)[i] * BINOM.inv_[i + 1];
    return Poly(std::forward<Poly>(r));
  }
  // assume a[0] = 1
  Poly log(int n) const {
    return (derivation() * inv(n)).integral().modXnR(n);
  }
  // assume a[0] = 0
  Poly exp(int n) const {
    Poly x(1);
    int k = 1;
    while (k < n) {
      k *= 2;
      x = (x * (Poly(1) - x.log(k) + this->modXn(k))).modXnR(k);
    }
    return x.modXnR(n);
  }
  Poly sin(int n) const {
    // 3 should be primitive root of valT::mod() = 4 x + 1
    static const valT COMPLEXI = pow(valT(3), (valT::mod() - 1) / 4);
    auto A = *this;
    for (auto& x : A) x *= COMPLEXI;
    A = A.exp(n);
    A -= A.inv(n);
    auto m = -COMPLEXI * INV2;
    for (auto& x : A) x *= m;
    return A;
  }
  Poly cos(int n) const {
    // 3 should be primitive root of valT::mod() = 4 x + 1
    static const valT COMPLEXI = pow(valT(3), (valT::mod() - 1) / 4);
    auto A = *this;
    for (auto& x : A) x *= COMPLEXI;
    A = A.exp(n);
    A += A.inv(n);
    for (auto& x : A) x *= INV2;
    return A;
  }
  // assume a[0] = 0
  Poly asin(int n) const { // unkown for acos
    auto D = this -> derivation();
    auto A = (Poly(1) - (*this) * (*this)).modXnR(n - 1);
    D = (D * A.sqrt(n - 1).inv(n - 1)).modXnR(n - 1);
    return D.integral();
  }
  // assume a[0] = 0
  Poly atan(int n) const {
    auto D = this -> derivation();
    auto A = (Poly(1) + (*this) * (*this)).modXnR(n - 1);
    D = (D * A.inv(n - 1)).modXnR(n - 1);
    return D.integral();
  }
  // assum a[0] = 1
  Poly sqrt(int n) const {
    Poly x(1);
    int k = 1;
    while (k < n) {
      k *= 2;
      x += this->modXn(k) * x.inv(k);
      x = x.modXnR(k) * INV2;
    }
    return x.modXnR(n);
  }
  // compose poly common algorithm F(A(x)) in $O(n^2)$, however Brent-Kung algorithm with $(n \log n)^{1.5}$ may be slower.
  Poly compose(Poly A, int n) const {
    A = A.modXnR(n);
    int sn = std::sqrt(n);
    std::vector<Poly> h(sn);
    h[0] = {1};
    for (int i = 1; i < sn; ++i) h[i] = (h[i - 1] * A).modXnR(n);
    auto B = (h.back() * A).modXnR(n);
    Poly now{1}, ans;
    for (int i = 0; i < n; i += sn) {
      Poly sm;
      for (int j = 0; j < sn && i + j < n; ++j) {
        auto m = this->at(i + j);
        auto tmp = h[j];
        for (auto& x : tmp) x *= m;
        sm += tmp;
      }
      ans += (now * sm).modXnR(n);
      now = (now * B).modXnR(n);
    }
    return ans;
  }
  // compose inverse (assmue a[0] = 0, a[1] \neq 0) based on Langrange: $[x^n]G(x) = \frac{1}{n}[x^{n - 1}](\frac{x}{F(x)})^n$
  Poly composeInv(int n) const {
    auto A = this->divXn(1).inv(n - 1);
    int sn = std::sqrt(n);
    std::vector<Poly> h(sn + 1);
    h[0] = {1};
    for (int i = 1; i <= sn; ++i) h[i] = (h[i - 1] * A).modXnR(n - 1);
    std::vector<valT> ans(n), inv(n);
    auto M = valT::mod();
    inv[1] = 1;
    for (int i = 2; i < n; ++i) inv[i] = inv[M % i] * valT(M - M / i);
    Poly now{1};
    for (int i = 0; i < n; i += sn) {
      for (int j = 1; j <= sn && i + j < n; ++j) {
        valT tmp;
        auto& sg = h[j];
        for (int k = 0, sk = i + j - 1; k <= sk; ++k) {
          tmp += now.at(k) * sg.at(sk - k);
        }
        ans[i + j] = tmp * inv[i + j];
      }
      now = (now * h.back()).modXnR(n - 1);
    }
    return Poly(std::forward<Poly>(ans));
  }
  Poly toFallingPowForm() {
    int n = (int)this->size();
    std::vector<valT> x(n);
    for (int i = 0; i < n; ++i) x[i] = i;
    auto y = this->evals(x);
    auto tmp = BINOM.ifac_;
    for (int i = 1; i < n; i += 2) tmp[i] = -tmp[i];
    Poly A = Poly(std::forward<Poly>(y)) * Poly(std::forward<Poly>(tmp));
    return A.modXnR(n);
  }
  Poly fromFallingPowForm() {
    int n = (int)this->size();
    Poly A = ((*this) * Poly(BINOM.ifac_)).modXnR(n);
    std::vector<valT> x(n), y = A;
    for (int i = 0; i < n; ++i) x[i] = i;
    y.resize(n);
    for (int i = 0; i < n; ++i) y[i] *= BINOM.fac_[i];
    return Lagrange(x, y);
  }
  valT eval(valT x) const {
    valT r(0), t(1);
    for (int i = 0, n = (int)this->size(); i < n; ++i) {
      r += (*this)[i] * t;
      t *= x;
    }
    return r;
  }
  // transpose {\rm MULT}(F(x),G(x))=\sum_{i\ge0}(\sum_{j\ge 0}f_{i+j}g_j)x^i
  Poly mulT(Poly&& rhs) const {
    if (rhs.size() == 0) return Poly();
    int n = rhs.size();
    rhs.reverse();
    return ((*this) * rhs).divXn(n - 1);
  }
  // multi-evaluation(new tech)
  std::vector<valT> evals(const std::vector<valT>& x) const {
    if ((int)this->size() == 0) return std::vector<valT>(x.size());
    int n = (int)x.size();
    std::vector<valT> ans(n);
    std::vector<Poly> g(4 * n);
    std::function<void(int, int, int)> build = [&](int l, int r, int p) {
      if (r - l == 1) {
        g[p] = Poly({1, -x[l]});
      } else {
        int m = (l + r) / 2;
        build(l, m, 2 * p);
        build(m, r, 2 * p + 1);
        g[p] = g[2 * p] * g[2 * p + 1];
      }
    };
    build(0, n, 1);
    std::function<void(int, int, int, const Poly& )> solve = [&](int l, int r, int p, const Poly& f) {
      if (r - l == 1) {
        ans[l] = f.at(0);
      } else {
        int m = (l + r) / 2;
        solve(l, m, 2 * p, f.mulT(std::move(g[2 * p + 1])).modXnR(m - l));
        solve(m, r, 2 * p + 1, f.mulT(std::move(g[2 * p])).modXnR(r - m));
      }
    };
    solve(0, n, 1, mulT(g[1].inv((int)this->size())).modXnR(n));
    return ans;
  } // https://www.luogu.com.cn/problem/P5050

  // \sum_{i = 0}^{n - 1} a_i / (1 - b_i x)
  static std::vector<valT> sumFraction(std::vector<valT> a, std::vector<valT> b, int N) {
    std::function<std::pair<Poly, Poly>(int, int)> solve = [&](int l, int r) -> std::pair<Poly, Poly> {
      if (r - l == 1) return {Poly(a[l]), Poly({1, - b[l]})};
      int m = (l + r) / 2;
      auto [pl, ql] = solve(l, m);
      auto [pr, qr] = solve(m, r);
      return {pl * qr + pr * ql, ql * qr};
    };
    auto [p, q] = solve(0, a.size());
    p *= q.inv(N);
    std::vector<valT> ans = p;
    ans.resize(N);
    return ans;
  } // https://codeforces.com/gym/102978/problem/D

  // compute $h(m), \cdots, h(m + cnt - 1)$ accroding to $h(0), h(1), \cdots, h(d)$
  static std::vector<valT> valToVal(std::vector<valT> h, valT m, int cnt) { // m > h.size()
    int d = (int)h.size() - 1;
    for (int i = 0; i <= d; ++i) {
      h[i] *= BINOM.ifac_[i] * BINOM.ifac_[d - i];
      if ((d - i) & 1) h[i] = -h[i];
    }
    std::vector<valT> f(d + cnt);
    auto now = m - valT(d);
    for (int i = 0; i < d + cnt; ++i) {
      f[i] = now.inv();
      ++now;
    }
    h = Poly(f) * Poly(h);
    h.resize(d + cnt);
    h = std::vector<valT>(h.begin() + d, h.end());
    now = 1;
    for (int i = 0; i <= d; ++i) now *= m - valT::raw(i);
    h[0] *= now;
    for (int i = 1; i < cnt; ++i) {
      now *= m + valT::raw(i);
      now *= (m + valT(i - d - 1)).inv();
      h[i] *= now;
    }
    return h;
  }; // https://www.luogu.com.cn/problem/P5667
  
  static Poly Lagrange(std::vector<valT> x, std::vector<valT> y) {
    std::function<Poly(int l, int r)> mulP = [&](int l, int r) -> Poly {
      if (r - l == 1) return Poly({-x[l], 1});
      int m = (l + r) / 2;
      return mulP(l, m) * mulP(m, r);
    };
    int n = (int)x.size();
    auto A = mulP(0, n).derivation();
    auto z = A.evals(x);
    for (int i = 0; i < n; ++i) y[i] /= z[i];
    std::function<std::pair<Poly, Poly>(int, int)> solve = [&](int l, int r) -> std::pair<Poly, Poly> {
      if (r - l == 1) {
        return {Poly(y[l]), Poly({-x[l], 1})};
      }
      int m = (l + r) / 2;
      auto [pl, ql] = solve(l, m);
      auto [pr, qr] = solve(m, r);
      return {pl * qr + pr * ql, ql * qr};
    };
    auto [p, q] = solve(0, x.size());
    return p;
  }

  // $a_n = \sum_{i = 1}^{k} f_i a_{n - i}$: https://oi-wiki.org/math/linear-recurrence/
  // find n-th term of The recursive formula for the constant coefficient of order k in $O(k \log k \log n)$
  static valT linearRecursion(const std::vector<valT>& a, std::vector<valT> f, LL n) {
    if (n < (int)a.size()) return a[n];
    int m = (int)f.size();
    std::reverse(f.begin(), f.end());
    std::vector<valT> g(m);
    g.emplace_back(1);
    Poly A = Poly({0, 1}), p = Poly(std::forward<Poly>(g)) - Poly(std::forward<Poly>(f));
    Poly R = A.powModPoly(n, p);
    return R.inner(a);
  } // https://www.luogu.com.cn/problem/P4723

  // ans[i] = 1^i + 2^i + ... + (n - 1)^i, 0 < i < k
  static std::vector<valT> prefixPowSum(int n, int k) {
    // Poly Numerator = Poly({0, n}).exp(k + 1).divXn(1);
    // Poly denominator  = Poly({0, 1}).exp(k + 1).divXn(1);
    std::vector<valT> a(k), b(k);
    for (int i = 0; i < k; ++i) a[i] = b[i] = BINOM.ifac_[i + 1];
    valT cur = 1;
    for (int i = 0; i < k; ++i) a[i] *= (cur *= valT::raw(n));
    auto Numerator = Poly(std::forward<Poly>(a)), denominator = Poly(std::forward<Poly>(b));

    auto f = (Numerator * denominator.inv(k)).modXnR(k) - Poly(1);
    auto ans = f;
    ans.resize(k);
    valT now(1);
    for (int i = 2; i < k; ++i) {
      now *= valT(i);
      ans[i] *= now;
    }
    return ans;
  }
  // $\prod_{i = 0}^{n - 1} (x + i)$ in $O(n \log n)$
  static Poly prod(int n) {
    std::function<Poly(int)> solve = [&](int n) -> Poly {
      if (n == 1) return Poly({0, 1});
      int k = n / 2;
      auto A = solve(k);
      std::vector<valT> tmp(k + 1);
      valT now{1};
      for (int i = 0; i <= k; ++i) {
        tmp[i] = now * BINOM.ifac_[i];
        now *= k;
      }
      auto B = A;
      for (int i = 0, nb = (int)B.size(); i < nb; ++i) B[i] *= BINOM.fac_[i];
      B = B.mulT(tmp).modXnR(k + 1);
      for (int i = 0, nb = (int)B.size(); i < nb; ++i) B[i] *= BINOM.ifac_[i];
      A *= B;
      if (2 * k != n) {
        B = A;
        for (auto& x : B) x *= valT::raw(n - 1);
        A = A.mulXn(1) + B;
      }
      return A;
    };
    return solve(n);
  }
  static Poly prodS(int n) { // $O(n \log^2 n)$
    std::function<Poly(int l, int r)> solve = [&](int l, int r) -> Poly {
      if (r - l == 1) return Poly({l, 1});
      int m = (l + r) / 2;
      return solve(l, m) * solve(m, r);
    };
    return solve(0, n);
  }
  // Stirling number
  static std::vector<valT> stirling1row(int n) {
    std::vector<valT> B = prod(n);
    B.resize(n + 1);
    return B;
  }
  std::vector<valT> stirling1col(int n, int k) {
    if (k > n)  return std::vector<valT>(n + 1);
    auto B = Poly({1, -1}).log(n + 2 - k).divXn(1);
    B = (-B).log(n + 1 - k);
    for (auto& x : B) x *= valT::raw(k);
    std::vector<valT> ans = B.exp(n + 1 - k).mulXn(k);
    ans.resize(n + 1);
    auto ifacK = BINOM.ifac_[k];
    for (int i = 0; i <= n; ++i) ans[i] *= BINOM.fac_[i] * ifacK;
    return ans;
  }
  std::vector<valT> stirling2row(int n) {
    auto tmp = BINOM.ifac_, a = BINOM.ifac_;
    for (int i = 1; i <= n; i += 2) tmp[i] = -tmp[i];
    for (int i = 0; i <= n; ++i) a[i] *= pow(valT::raw(i), n);
    std::vector<valT> ans = Poly(std::forward<Poly>(a)) * Poly(std::forward<Poly>(tmp));
    ans.resize(n + 1);
    return ans;
  }
  std::vector<valT> stirling2col(int n, int k) {
    if (k > n)  return std::vector<valT>(n + 1);
    auto A = Poly(BINOM.ifac_);
    A = A.divXn(1).modXnR(n + 1 - k);
    A = A.log(n + 1 - k);
    for (auto& x : A) x *= valT::raw(k);
    A = A.exp(n + 1 - k).mulXnR(k);
    std::vector<valT> ans = A;
    ans.resize(n + 1);
    auto ifacK = BINOM.ifac_[k];
    for (int i = 0; i <= n; ++i) ans[i] *= BINOM.fac_[i] * ifacK;
    return ans;
  }
}; // https://www.luogu.com.cn/training/3015#information


template<typename valT, typename enable = ModT<valT>>
class PolyBase : public std::vector<valT> {
 protected:
  void standard() {
    while (!this->empty() && this->back() == valT(0)) this->pop_back();
  }
  void reverse() {
    std::reverse(this->begin(), this->end());
    this->standard();
  }
 public:
  PolyBase() {}
  PolyBase(const valT& x) : std::vector<valT>{x} { standard();}
  PolyBase(const std::vector<valT>& a) : std::vector<valT>{a} { standard();}
  PolyBase(std::vector<valT>&& a) : std::vector<valT>{std::move(a)} { standard();}
  valT at(int id) const {
    if (id < 0 || id >= (int)this->size()) return 0;
    return (*this)[id];
  }
};


// There will be `PolyNTT`, `PolyFFT`, `PolyFFTDynamic`, `PolyMFT` provided to suit for different module $M$.

// - PolyMFT: $M > \text{INT_MAX}$
// - PolyFFTDynamic: else if $M$ is uncertain.
// - PolyNTT: else if $M$ is fixed NTT-friendly, such as $M = 998244353$
// - PolyFFT: else
// - PolyOrigin for testing
