#pragma once
#include <bits/stdc++.h>
#include "mod.hpp"
using LL = long long;

// n choose k: 1 stand for choosen
void GospersHack(int n, int k) {
  int cur = (1 << k) - 1;
  while (!(cur >> n)) {
    // do something below

    // do something above
    int lb = __builtin_ctz(cur);
    int r = cur + (1 << lb);
    cur = ((r ^ cur) >> (2 + lb)) | r;
  }
}
void GospersHackS(int n, int k) {
  int cur = (1 << k) - 1;
  while (!(cur >> n)) {
    // do something below

    // do something above
    int lb = cur & -cur;
    int r = cur + lb;
    cur = (((r ^ cur) >> 2) / lb) | r;
  }
}

// |(i, j) : 0 < i, j < n, i + j = d|
int twoSumCount(int n, int d) {
  return std::max(0, std::min(d - 1, 2 * n - d - 1));
}
// https://atcoder.jp/contests/abc220/tasks/abc220_e

// knuth Shuffle: you may use std::shuffle instead
template<typename T>
void KnuthShuffle(std::vector<T>& a) {
  std::mt19937 rnd(std::chrono::steady_clock::now().time_since_epoch().count());
  for (int i = (int)a.size() - 1; i > 0; --i) {
    std::swap(a[i], a[rnd() % (i + 1)]);
  }
}
// choose m distinct elements from [0, n) with equal probability
std::vector<int> uniformChoose(int n, int m) {
  std::mt19937 rnd(std::chrono::steady_clock::now().time_since_epoch().count());
  std::map<int, int> mp;
  std::vector<int> ans;
  ans.reserve(m);
  for (int i = 0; i < m; ++i) {
    int x = rnd() % (n - i);
    ans.emplace_back(mp.count(x) ? mp[x] : x);
    mp[x] = mp.count(i) ? mp[i] : i;
  }
  return ans;
}

// Fibonacci[n] % M
int Fib(int n, int M) {
  int a = 0, b = 1, c = 1, d = 0;
  while (n) {
    if (n & 1) {
      int x = 1LL * a * c % M;
      a = (1LL * a * d + 1LL * b * c + x) % M;
      b = (1LL * b * d + x) % M;
    }
    n >>= 1;
    int x = 1LL * c * c % M;
    c = (2LL * c * d + x) % M;
    d = (1LL * d * d + x) % M;
  }
  return a;
}

// $\sum_{i = 0}^{n - 1} \lfloor \frac{a \cdot i + b}{m} \rfloor$ in  $O(\log n)$
LL floorSum(int n, int m, int a, int b) {
  LL r = 0;
  if (a >= m) {
    r += 1LL * a / m * (n - 1) * n / 2;
    a %= m;
  }
  if (b >= m) {
    r += 1LL * b / m * n;
    b %= m;
  }
  int yMax = (1LL * a * n + b) / m;
  if (yMax == 0) return r;
  r += 1LL * (n - 1) * yMax;
  r -= floorSum(yMax, a, m, m - b - 1);
  return r;
}
// https://atcoder.jp/contests/practice2/tasks/practice2_c

// $\sum_{\sum c_i x_i = m} \frac{(\sum x_i)!}{\prod (x_i !)}$
int sumNum(const std::vector<int>& c, int m, int M) {
  std::vector<int> dp(m + 1);
  dp[0] = 1;
  for (int i = 1; i <= m; ++i) {
    for (auto x : c) if (x <= i) {
      dp[i] += dp[i - x];
      if (dp[i] >= M) dp[i] -= M;
    }
  }
  return dp[m];
}

// count min time: every time --n or ++m s.t. $n | m$. $O(\sqrt{m})$
int decInc(int n, int m) {
  if (n > m) return n - m;
  int ans = n - 1;
  // let n be i, then we need n - i + (m % i ? i - m % i : 0) <= n - m % i
  // Note that m % i = m - (m / i) * i
  for (int i = 1; i <= n; i = m / (m / i) + 1) {
    ans = std::min(ans, n - m % i);
  }
  // sepcial judge for the case m % i == 0
  for (int i = 1; i <= n && i * i <= m; ++i) if (m % i == 0) {
    ans = std::min(ans, n - i);
    if (m / i <= n) ans = std::min(ans, n - m / i);
  }
  return ans;
}

// finds min x s.t. L <= (A * x) % M <= R (or -1 if it does not exist)
int FirstInRange(int a, int m, int l, int r) { // 0 <= L <= R < M, 0 < a < M
  if (l == 0) return 0;
  std::function<std::pair<int, int>(int, int)> dfs = [&](int a, int m) -> std::pair<int, int> {
    if (a == 0) return std::pair(-1, 0);
    if (a > m - a) {
      std::tie(l, r) = std::pair(m - r, m - l);
      auto [x, y] = dfs(m - a, m);
      return std::pair(x, x - y - 1);
    }
    int k = (l + a - 1) / a;
    if (k * a <= r) return std::pair(k, 0);
    --k;
    l -= k * a;
    r -= k * a;
    int z = (m + a - 1) / a;
    auto [x, y] = dfs(z * a - m, a);
    return x == -1 ? std::pair(-1, -1) : std::pair(k + x * z - y, x);
  };
  return dfs(a, m).first;
}

template<typename T> // use ModInt, MInt, ModLL
class RecSeq : public std::vector<T> {
  // x^n = c_0 + c_1 x + \cdots c_{n - 1} x^{n - 1}
  static inline std::vector<T> c_;
public:
  using std::vector<T>::vector;
  static void setRec(std::vector<T> c) { c_ = std::move(c);}
	RecSeq operator*(const RecSeq& A) const {
    int n = (int)this->size(), m = (int)A.size();
		RecSeq R(n + m - 1);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
        R[i + j] += (*this)[i] * A[j];
      }
    }
		for (int i = n + m - 2, cn = c_.size(); i >= cn; --i) {
      for (int j = 0; j < cn; ++j) {
        R[i - cn + j] += R[i] * c_[j];
      }
		}
    if (R.size() > c_.size()) {
      R.resize(c_.size());
      R.shrink_to_fit();
    }
		return R;
	}
  friend RecSeq pow(RecSeq A, int n) {
    RecSeq R{1};
    while (n) {
      if (n & 1) R = R * A;
      n >>= 1;   A = A * A;
    }
    return R;
  }
};

// Gauss-Jordan Elimination $Ax = b$, float version
std::vector<double> Gauss(std::vector<std::vector<double>> A, std::vector<double> b) {
  int n = (int)A.size(), m = (int)A[0].size();
  assert(m == (int)b.size());
  std::vector<double> x(m);
  std::vector<int> p(m);
  std::iota(p.begin(), p.end(), 0);
  const double eps = 1e-12;
  auto findNonZero = [&](int i) { // find nonzero element with max $Abs(A[j][i])$
    int ans = i;
    for (int row = i; row < n; ++row) if (fabs(A[row][i]) > fabs(A[ans][i])) ans = row;
    return fabs(A[ans][i]) > eps ? ans : n;
  };
  auto triangleGauss = [&](int sz) { // A[i][i] = 1
    std::vector<double> x(sz);
    for (int i = sz - 1; i >=0; --i) {
      x[i] = b[i];
      for (int row = 0; row < i; ++row) b[row] -= A[row][i] * x[i];
    }
    x.resize(A[0].size());
    return x;
  };
  int sz = n;
  for (int i = 0, row = 0; i < n; ++i) {
    while (i < m) {
      row = findNonZero(i);
      if (row != n) break;
      for (int j = 0; j < n; ++j) A[j][i] = A[j][m - 1];
      std::swap(p[i], p[--m]);
    }
    if (i == m) {
      for (int row = m; row < n; ++row) if (fabs(b[row]) > eps) {
        // std::cout << "\nNo answer\n";
        return std::vector<double>();
      }
      sz = i;
      break;
    }
    if (row != i) {
      std::swap(A[row], A[i]);
      std::swap(b[row], b[i]);
    }
    b[i] /= A[i][i];
    for (int j = m - 1; j >= i; --j) A[i][j] /= A[i][i];
    for (int row = i + 1; row < n; ++row) {
      b[row] -= A[row][i] * b[i];
      for (int j = m - 1; j >= i; --j) {
        A[row][j] -= A[row][i] * A[i][j];
      }
    }
  }
  // if (sz != A[0].size()) std::cout << "\nInfinite answer\n";
  auto xt = triangleGauss(sz);
  for (int t = 0, na = A[0].size(); t < na; ++t) x[p[t]] = xt[t];
  return x;
}

// xor version of Gauss-Jordan Elimination with 0-1-matrix A
template <typename T>
std::vector<T> GaussXor(std::vector<std::vector<T>> A, std::vector<T> b) {
  int n = (int)A.size(), m = (int)A[0].size();
  std::vector<T> x(m);
  std::vector<int> p(m);
  std::iota(p.begin(), p.end(), 0);
  auto triangleGauss = [&](int sz) { // A[i][i] = 1
    std::vector<T> x(sz);
    for (int i = sz - 1; i >=0; --i) {
      x[i] = b[i];
      for (int row = 0; row < i; ++row) if (A[row][i]) {
        b[row] ^= x[i];
      }
    }
    x.resize(m);
    return x;
  };
  auto findNonZero = [&](int i) {
    for (int row = i; row < n; ++row) if (A[row][i]) return row;
    return n;
  };
  int sz = n;
  for (int i = 0, row = 0; i < n; ++i) {
    while (i < m) {
      row = findNonZero(i);
      if (row != n) break;
      for (int j = 0; j < n; ++j) A[j][i] = A[j][m - 1];
      std::swap(p[i], p[--m]);
    }
    if (i == m) {
      for (int row = m; row < n; ++row) if (b[row]) {
        // std::cout << "\nNo answer\n";
        return {};
      }
      sz = i;
      break;
    }
    if (row != i) {
      std::swap(A[i], A[row]);
      std::swap(b[i], b[row]);
    }
    for (int row = i + 1; row < n; ++row) if (A[row][i]) {
      b[row] ^= b[i];
      for (int j = m - 1; j >= i; --j) A[row][j] ^= A[i][j];
    }
  }
  // if (sz != m) std::cout << "\nInfinite answer\n";
  auto xt = triangleGauss(sz);
  for (int t = 0; t < m; ++t) x[p[t]] = xt[t];
  return x;
}

// mod version of Gauss-Jordan Elimination
template<typename valT, typename enable = ModT<valT>>
std::vector<valT> GaussModp(std::vector<std::vector<valT>> A, std::vector<valT> b) {
  int n = (int)A.size(), m = (int)A[0].size();
  assert(m == (int)b.size());
  std::vector<valT> x(m);
  std::vector<int> p(m);
  std::iota(p.begin(), p.end(), 0);
  auto triangleGauss = [&](int sz) { // A[i][i] = 1
    std::vector<valT> x(sz);
    for (int i = sz - 1; i >=0; --i) {
      x[i] = b[i];
      for (int row = 0; row < i; ++row) b[row] -= A[row][i] * x[i];
    }
    x.resize(m);
    return x;
  };
  auto findNonZero = [&](int i) {
    for (int row = i; row < n; ++row) if (A[row][i]) return row;
    return n;
  };
  int sz = n;
  for (int i = 0, row = 0; i < n; ++i) {
    while (i < m) {
      row = findNonZero(i);
      if (row != n) break;
      for (int j = 0; j < n; ++j) A[j][i] = A[j][m - 1];
      std::swap(p[i], p[--m]);
    }
    if (i == m) {
      for (int row = m; row < n; ++row) if (b[row]) {
        // std::cout << "\nNo answer\n";
        return {};
      }
      sz = i;
      break;
    }
    if (row != i) {
      std::swap(A[i], A[row]);
      std::swap(b[i], b[row]);
    }
    auto invA = A[i][i].inv();
    b[i] *= invA;
    for (int j = m - 1; j >= i; --j) A[i][j] *= invA;
    for (int row = i + 1; row < n; ++row) {
      b[row] -= A[row][i] * b[i];
      for (int j = m - 1; j >= i; --j) A[row][j] -= A[row][i] * A[i][j];
    }
  }
  // if (sz != m) std::cout << "\nInfinite answer\n";
  auto xt = triangleGauss(sz);
  for (int t = 0; t < m; ++t) x[p[t]] = xt[t];
  return x;
}

// Simplex algorithm
using VD = std::vector<double>;
const double eps = 1e-10;
const double inf = 1e10;
// make sure that A = (I, A') and b >= 0, compute max cx
VD simplexCore(VD c, std::vector<VD> A, VD b) {
  int n = (int)A.size(), m = (int)c.size();
  std::vector<int> p(m);
  std::iota(p.begin(), p.end(), 0);
  for (int i = 0; i < n; ++i) A[i].emplace_back(b[i]);
  c.emplace_back(0);
  A.emplace_back(c);
  for (int j = n; j <= m; ++j) {
    for (int i = 0; i < n; ++i) {
      A[n][j] -= A[n][i] * A[i][j];
    }
  }
  auto check = [&]() -> bool {
    for (int j = n; j < m; ++j) if (A[n][j] > eps) {
      bool flag = false;
      for (int i = 0; i < n; ++i) if (A[i][j] > eps) {
        flag = true;
        break;
      }
      if (!flag) return false;
    }
    return true;
  };
  while (1) {
    int ch = std::max_element(A[n].begin() + n, A[n].begin() + m) - A[n].begin(), hc;
    if (A[n][ch] < eps) break;
    assert(check()); // otherwise unbounded, no max solution
    double theta = DBL_MAX;
    for (int i = 0; i < n; ++i) if (A[i][ch] > eps && A[i].back() / A[i][ch] < theta) {
      theta = A[i].back() / A[i][ch];
      hc = i;
    }
    std::swap(p[ch], p[hc]);
    double tmp = 1 / A[hc][ch];
    for (int j = n; j <= m; ++j) A[hc][j] *= tmp;
    for (int i = 0; i <= n; ++i) if (i != hc) {
      for (int j = n; j <= m; ++j) if (j != ch) {
        A[i][j] -= A[i][ch] * A[hc][j];
      }
    }
    for (int i = 0; i <= n; ++i) A[i][ch] *= -tmp;
    A[hc][ch] = tmp;
  }
  VD x(m);
  for (int i = 0; i < n; ++i) x[p[i]] = A[i].back();
  // watch(-A.back().back()); // max_val
  return x; // point Corresponds to max_val
}
// compute max cx, with Aqx = bq and Alq x <= blq, end of 0 can be ommit in A and Aq
VD simplex(VD c, std::vector<VD> Aq, VD bq, std::vector<VD> Alq, VD blq) {
  assert(Aq.size() == bq.size());
  assert(Alq.size() == blq.size());
  int n = Aq.size() + Alq.size();
  int m = (int)c.size();
  for (int i = 0, nb = bq.size(); i < nb; ++i) if (bq[i] < -eps) {
    for (auto& x : Aq[i]) x = -x;
    bq[i] = -bq[i];
  }
  for (int i = 0, nb = blq.size(); i < nb; ++i) if (blq[i] < -eps) {
    for (auto& x : Alq[i]) x = -x;
    ++m;
  }
  std::vector<VD> A(n, VD(n + m));
  VD f(n + m), b(n);
  int now = n + c.size();
  for (int i = 0; i < n; ++i) A[i][i] = 1;
  for (int i = 0, na = Aq.size(); i < na; ++i) {
    for (int j = 0; j < na; ++j) A[i][n + j] = Aq[i][j];
    b[i] = bq[i];
    f[i] = -inf;
  }
  for (int i = 0, na = Alq.size(); i < na; ++i) {
    for (int j = 0; j < na; ++j) A[i + Aq.size()][n + j] = Alq[i][j];
    if (blq[i] < -eps) {
      A[i + Aq.size()][now++] = -1;
      f[i + Aq.size()] = -inf;
    }
    b[i + Aq.size()] = fabs(blq[i]);
  }
  for (int i = 0, nc = (int)c.size(); i < nc; ++i) f[n + i] = c[i];
  auto x = simplexCore(f, A, b);
  return VD(x.begin() + n, x.begin() + n + c.size());
}


// Polynomial multiplication with arbitrary modulus $O(n^{\log_2 3})$
using VL = std::vector<LL>;
VL Karatsuba(VL a, VL b, LL p) {
  if (a.size() < b.size()) std::swap(a, b);
  auto mulS = [&](VL a, VL b) {
    int n = (int)a.size(), m = (int)b.size(), sz = n + m - 1;
    std::vector<__int128_t> c(sz);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
        c[i + j] += a[i] * b[j];
      }
    }
    VL r(sz);
    for (int i = 0; i < sz; ++i) r[i] = c[i] % p;
    return r;
  };
  constexpr int N = 65;
  std::function<VL(VL, VL, int)> mul = [&](VL a, VL b, int n) -> VL {
    if (n < N) return mulS(a, b);
    int n2 = n / 2, n1 = n - 1;
    VL a2 = VL(a.begin() + n2, a.end());
    VL b2 = VL(b.begin() + n2, b.end());
    a.resize(n2); b.resize(n2);
    VL ap = a2, bp = b2;
    for (int i = 0; i < n2; ++i) if ((ap[i] += a[i]) >= p) ap[i] -= p;
    for (int i = 0; i < n2; ++i) if ((bp[i] += b[i]) >= p) bp[i] -= p;
    VL ab = mul(a, b, n2);
    VL a2b = mul(ap, bp, n2);
    VL a2b2 = mul(a2, b2, n2);
    for (int i = 0; i < n1; ++i) {
      if ((a2b[i] -= ab[i]) < 0) a2b[i] += p;
      if ((a2b[i] -= a2b2[i]) < 0) a2b[i] += p;
    }
    auto r = ab;
    r.emplace_back(0);
    r.insert(r.end(), a2b2.begin(), a2b2.end());
    for (int i = 0; i < n1; ++i) if ((r[i + n2] += a2b[i]) >= p) r[i + n2] -= p;
    return r;
  };
  int n = (int)a.size(), m = (int)b.size(), tot = std::max(1, n + m - 1);
  if (m < N || n / m * 2 > m) return mulS(a, b);
  int sz = 1 << std::__lg(tot * 2 - 1);
  while (sz < n) sz *= 2;
  a.resize(sz), b.resize(sz);
  auto r = mul(a, b, sz);
  r.resize(tot);
  return r;
} // https://www.luogu.com.cn/problem/P4245


// // Polynomial multiplication with arbitrary modulus $O(n^{\log_2 3})$  Parallel version
// using VL = std::vector<LL>;
// VL KaratsubaParallel(VL a, VL b, LL p) {
//   if (a.size() < b.size()) std::swap(a, b);
//   auto mulS = [&](VL a, VL b) {
//     int n = (int)a.size(), m = (int)b.size(), sz = n + m - 1;
//     std::vector<__int128_t> c(sz);
//     for (int i = 0; i < n; ++i) {
//       for (int j = 0; j < m; ++j) {
//         c[i + j] += a[i] * b[j];
//       }
//     }
//     VL r(sz);
//     for (int i = 0; i < sz; ++i) r[i] = c[i] % p;
//     return r;
//   };
//   constexpr int N = 65;
//   std::function<VL(VL, VL, int)> mul = [&](VL a, VL b, int n) -> VL {
//     if (n < N) return mulS(a, b);
//     int n2 = n / 2, n1 = n - 1;
//     VL a2 = VL(a.begin() + n2, a.end());
//     VL b2 = VL(b.begin() + n2, b.end());
//     a.resize(n2); b.resize(n2);
//     VL ap = a2, bp = b2;
//     for (int i = 0; i < n2; ++i) if ((ap[i] += a[i]) >= p) ap[i] -= p;
//     for (int i = 0; i < n2; ++i) if ((bp[i] += b[i]) >= p) bp[i] -= p;
//     std::future<VL> abThread = std::async(mul, a, b, n2);
//     VL a2b = mul(ap, bp, n2);
//     VL ab = abThread.get();
//     VL a2b2 = mul(a2, b2, n2);
//     for (int i = 0; i < n1; ++i) {
//       if ((a2b[i] -= ab[i]) < 0) a2b[i] += p;
//       if ((a2b[i] -= a2b2[i]) < 0) a2b[i] += p;
//     }
//     auto r = ab;
//     r.emplace_back(0);
//     r.insert(r.end(), a2b2.begin(), a2b2.end());
//     for (int i = 0; i < n1; ++i) if ((r[i + n2] += a2b[i]) >= p) r[i + n2] -= p;
//     return r;
//   };
//   int n = (int)a.size(), m = (int)b.size(), tot = std::max(1, n + m - 1);
//   if (m < N || n / m * 8 > m) return mulS(a, b);
//   int sz = 1 << std::__lg(tot * 2 - 1);
//   a.resize(sz), b.resize(sz);
//   auto r = mul(a, b, sz);
//   r.resize(tot);
//   return r;
// } // you may need: g++ -lpthread -o main main.cpp -O2 -std=c++17
// // https://www.luogu.com.cn/problem/P4245, the above code will CE since the platform don't support

// Segment DP: f_{l, r} = \min_{l \leq k < r} f_{l, k} + f_{k + 1, r} + w(l, r) \qquad (1 \leq l < r \leq n)
template<typename T>
std::vector<std::vector<T>> quadrangleItvDp(std::vector<std::vector<T>> w, int n) {
  std::vector<std::vector<T>> f(n + 1, std::vector<T>(n + 1)), mf(n + 1, std::vector<T>(n + 1));
  for (int i = 1; i < n; ++i) {
    f[i][i + 1] = w[i][i + 1];
    mf[i][i + 1] = i;
  }
  auto const inf = std::numeric_limits<T>::max() / 2;
  for (int len = 2; len < n; ++len) {
    for (int l = 1, r = len + 1; r <= n; ++l, ++r) {
      f[l][r] = inf;
      for (int k = mf[l][r - 1]; k <= mf[l + 1][r]; ++k) {
        if (f[l][r] > f[l][k] + f[k + 1][r]) {
          f[l][r] = f[l][k] + f[k + 1][r];
          mf[l][r] = k;
        }
      }
      f[l][r] += w[l][r];
    }
  }
  return f;
}

// roll DP: f_{i, j} = \min_{k < j} f_{i - 1, k} + w(k + 1, j) \quad (1 \leq i \leq n, 1 \leq j \leq m)
template<typename T>
std::vector<std::vector<T>> quadrangleRollDp(std::vector<std::vector<T>> w, int n, int m) {
// w is a (n + 1, n + 1) matrix, the answer is a (m + 1, n + 1) matrix
  std::vector<std::vector<T>> f(m + 1, std::vector<T>(n + 1)), mf(m + 1,  std::vector<T>(n + 2));
  auto const inf = std::numeric_limits<T>::max() / 2;
  for (int i = 1; i < n; ++i) f[0][i] = inf;
  for (int i = 1; i <= m; ++i) {
    mf[i][n + 1] = n;
    for (int j = n; j > 0; --j) {
      f[i][j] = inf;
      for (int k = std::max(i - 1, mf[i - 1][j]); k < j && k <= mf[i][j + 1]; ++k) {
        if (f[i][j] > f[i - 1][k] + w[k + 1][j]) {
          f[i][j] = f[i - 1][k] + w[k + 1][j];
          mf[i][j] = k;
        }
      }
    }
  }
  return f;
}
// https://www.luogu.com.cn/problem/P4767

class PalindromeNumber {
  std::vector<LL> ten{1};
  PalindromeNumber() {
    for (int i = 0; i < 18; ++i) ten.emplace_back(ten.back() * 10);
  }
  int len(LL n) {
    return std::upper_bound(ten.begin(), ten.end(), n) - ten.begin();
  }
  LL revBit(LL n) {
    LL r = 0;
    while (n) {
      r = r * 10 + n % 10;
      n /= 10;
    }
    return r;
  }
 public:
  PalindromeNumber(const PalindromeNumber& A) = delete;
  static PalindromeNumber& Instance() {
    static PalindromeNumber instance;
    return instance;
  }
  LL nthPalindrome(int k) {
    int i = 1;
    while (1) {
      if (k <= 9 * ten[i - 1]) {
        LL ans = ten[i - 1] + k - 1;
        return ten[i - 1] * ans + revBit(ans / 10);
      }
      k -= 9 * ten[i - 1];
      if (k <= 9 * ten[i - 1]) {
        LL ans = ten[i - 1] + k - 1;
        return ten[i] * ans + revBit(ans);
      }
      k -= 9 * ten[i - 1];
      ++i;
    }
  }
  int Palindrome(LL n) {  // numbers of Palindrome < n
    int x = len(n), x2 = x / 2, ans = 0;
    LL tmp = n / ten[x2];
    LL now = tmp * ten[x2] + revBit(x & 1 ? tmp / 10 : tmp);
    if (now >= n) --ans;
    return ans += tmp + ten[x2] - 1;
  }
  LL solve(LL n, int k) { return nthPalindrome(k + Palindrome(n)); }
};
// https://ac.nowcoder.com/acm/contest/11191/C

using ULL = unsigned long long;
// 40% faster
unsigned fastPowMod998244353(unsigned x, unsigned n) {
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

// 18% faster
unsigned fastPowMod1000000007(unsigned x, unsigned n) {
  static const unsigned m = 1000000007U;
  static const unsigned mr = 2226617417U;
  static const unsigned m1 = 294967268U;
  static const unsigned m1inv = 518424770U;
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

// 18% faster
unsigned fastPowMod1000000009(unsigned x, unsigned n) {
  static const unsigned m = 1000000009U;
  static const unsigned mr = 737024967U;
  static const unsigned m1 = 294967260U;
  static const unsigned m1inv = 171601999U;
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
// asm may helps https://www.rieselprime.de/ziki/Montgomery_multiplication


/*
// Definition for a Node.
class Node {
public:
  int val;
  Node* next;
  Node* random;
  Node(int _val) {
    val = _val;
    next = NULL;
    random = NULL;
  }
};
*/
template<typename Node>
Node* copyRandomList(Node* head) {
  if (head == nullptr) return nullptr;
  // insert Node
  auto* now = head;
  while (now) {
    auto p = new Node(now->val);
    auto* nxt = now->next;
    now->next = p;
    p->next = nxt;
    now = nxt;
  }
  // update random
  now = head;
  while (now) {
    auto* cur = now->next;
    if (now->random) {
      cur->random = now->random->next;
    } else {
      cur->random = nullptr;
    }
    now = cur->next;
  }
  // get ans and make head back to origin
  now = head;
  auto* ans = head->next;
  while (now) {
    auto* cur = now->next;
    now->next = cur->next;
    now = now->next;
    if (now) {
      cur->next = now->next;
    } else {
      cur->next = nullptr;
    }
  }
  return ans;
}
// https://leetcode-cn.com/submissions/detail/261932599/
