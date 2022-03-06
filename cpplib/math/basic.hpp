#pragma once
#include <bits/stdc++.h>
#include "mod.hpp"
#include "../template.hpp"
using LL = long long;

int powMod(int x, int n, int M) {
  int r = 1;
  while (n) {
    if (n&1) r = 1LL * r * x % M;
    n >>= 1; x = 1LL * x * x % M;
  }
  return r;
}

template<typename U, typename T, typename T2, typename enable = TwiceT<T, T2>>
T powModT(T x, U n, T M) {
  T r = 1;
  while (n) {
    if (n&1) r = T2(r) * x % M;
    n >>= 1; x = T2(x) * x % M;
  }
  return r;
}

template<typename T>
T floor(T a, T n) {
  if (n < 0) {
    n = -n;
    a = -a;
  }
  return a < 0 ? (a - n + 1) / n : a / n;
}
template<typename T>
T ceil(T a, T n) {
  if (n < 0) {
    n = -n;
    a = -a;
  }
  return a < 0 ? a / n : (a + n - 1) / n;
}

// never mixed it with cin and cout, useful for int128
template<typename T, typename enable = IntegerT<T>>
class FastIO {
 public:
  static T read(){
    T x = 0;
    bool negative = false;
    // you may use buffer instead for speed
    char ch = getchar();
    while (ch < '0' || ch > '9'){
      if (ch == '-') negative = true;
      ch = getchar();
    }
    while (ch >= '0' && ch <= '9') {
      x = x * 10 + ch - '0';
      ch = getchar();
    }
    return negative ?  -x : x;
  }
  static void print(T x){
    if (x < 0) {
      putchar('-');
      x = -x;
    }
    printCore(x);
  }
 private:
  static void printCore(T x){
    if (x > 9) printCore(x / 10);
    putchar(x % 10 + '0');
  }
};

// slightly faster than std::gcd
LL gcd(LL a, LL b) {
  if (!a || !b) return a | b;
  unsigned shift = __builtin_ctzll(a | b);
  a >>= __builtin_ctzll(a);
  do {
    b >>= __builtin_ctzll(b);
    if (a > b) std::swap(a, b);
    b -= a;
  } while (b);
  return a << shift;
}
// https://cp-algorithms.com/algebra/euclid-algorithm.html

// ax + by = gcd(a,b)
template<typename T, typename enable = SignedT<T>>
std::tuple<T, T, T> exGcd(T a, T b) {
  if (b == 0) return {a, 1, 0};
  auto [d, y, x] = exGcd(b, a % b);
  return {d, x, y - a / b * x};
}
// |x| <= max(1, b), |y| <= a ===> |y - a / b * x| <= a % b + a / b * b = a

// Chinese remainder theorem: x = ai mod mi, m_i > 0, 0 <= a_i < m_i
std::pair<LL, LL> crt2(LL a1, LL m1, LL a2, LL m2) {
  auto [d, t1, t2] = exGcd(m1, m2);
  assert((a2 - a1) % d == 0);
  LL m = m1 / d * m2;
  LL ans = (a1 + (a2 - a1) / d * t1 % m2 * m1) % m;
  return {ans < 0 ? ans + m: ans, m};
}

std::pair<LL, LL> crt(const std::vector<std::pair<LL, LL>>& A) {
  auto ans = A[0];
  for (int i = 1, na = (int)A.size(); i < na; ++i) {
    ans = crt2(ans.first, ans.second, A[i].first, A[i].second);
  }
  return ans;
}
// https://www.luogu.com.cn/problem/P1495

// O(N \log N) smalleset prime factor(may be faster)
std::vector<int> spfS(int N) {
  std::vector<int> sp(N);
  std::iota(sp.begin(), sp.end(), 0);
  for(int i = 2; i * i < N; ++i) if(sp[i] == i) {
    for(int j = i * i; j < N; j += i) if(sp[j] == j) {
      sp[j] = i;
    }
  }
  return sp;
}

// $O(N)$ smallest prime factor
std::vector<int> spf(int N) {
  std::vector<int> sp(N), p{0, 2};
  p.reserve(N);
  for (int i = 2; i < N; i += 2) sp[i] = 2;
  for (int i = 1; i < N; i += 2) sp[i] = i;
  for (int i = 3; i < N; i += 2) {
    if (sp[i] == i) p.emplace_back(i);
    for (int j = 2, np = (int)p.size(); j < np && p[j] <= sp[i] && i * p[j] < N; ++j) {
      sp[i * p[j]] = p[j]; // Note that sp[x] is assigned only once foreach x
    }
  }
  return sp;
}

// O(N) none square factor
std::vector<int> nsf(int N) {
  std::vector<int> ans(N);
  std::iota(ans.begin(), ans.end(), 0);
  for (int i = 1; i < N; ++i) if (ans[i] == i) {
    for (int j = 2; i * j * j < N; ++j) {
      ans[i * j * j] = i;
    }
  }
  return ans;
}

// O(N) none square factor
std::vector<int> nsfS(int N) {
  auto sp = spf(N);
  std::vector<int> ans(N);
  ans[1] = 1;
  for (int i = 2; i < N; ++i) {
    int si = i / sp[i];
    ans[i] = si % sp[i] == 0 ? ans[si / sp[i]] : ans[si] * sp[i];
  }
  return ans;
}

class Binom {
  static inline constexpr int N = 65;
  LL C_[N][N];
  Binom() {
    for (int i = 0; i < N; ++i) C_[i][0] = C_[i][i] = 1;
    for (int i = 1; i < N; ++i) {
      for (int j = 1; j < i; ++j) {
        C_[i][j] = C_[i - 1][j] + C_[i - 1][j - 1];
      }
    }
  }
 public:
  Binom(const Binom&) = delete;
  static Binom& Instance() {
    static Binom instance_;
    return instance_;
  }
  LL operator()(int m, int n) const {
    assert(n < N && m < N);
    return C_[m][n];
  }
};

template<typename valT, typename enable = ModT<valT>>
class BinomModp {
  static inline constexpr int N = 1e6 + 2;
  BinomModp() { fac_.reserve(N), ifac_.reserve(N), inv_.reserve(N);}
  void init(int n) {
    assert(n <= valT::mod());
    fac_[0] = 1;
    for (int i = 1; i < n; ++i) fac_[i] = fac_[i - 1] * valT::raw(i);
    ifac_[n - 1] = fac_[n - 1].inv();
    for (int i = n - 1; i > 0; --i) ifac_[i - 1] = ifac_[i] * valT::raw(i);
    for (int i = 1; i < n; ++i) inv_[i] = ifac_[i] * fac_[i - 1];
  }
 public:
  std::vector<valT> fac_, ifac_, inv_;
  BinomModp(const BinomModp&) = delete;
  static BinomModp& Instance(int N = 0) {
    static BinomModp instance_;
    if (N) instance_.init(N);
    return instance_;
  }
  valT binom(int n, int k) const {
    if (n < 0 || n < k) return valT(0);
    return fac_[n] * ifac_[k] * ifac_[n - k];
  }
  // M is a small prime number in this case
  valT lucas(int n, int k) const {
    valT r(1);
    const int M = valT::mod();
    while (n && k) {
      int np = n % M, kp = k % M;
      if (np < kp) return valT(0);
      r *= binom(np, kp);
      n /= M; k /= M;
    }
    return r;
  }
};

// Calculate f(m) where f is the Lagrange interpolation on $f(0), f(1), \cdots, f(n - 1)$
template<typename valT, typename enable = ModT<valT>>
valT Lagrange(const std::vector<valT>& f, int m) {
  int n = (int)f.size();
  if (m < n) return f[m];
  auto& B = BinomModp<valT>::Instance(n);
  std::vector<valT> AP(n), BP(n);
  AP[0] = BP[n - 1] = valT(1);
  for (int i = 1; i < n; ++i) AP[i] = AP[i - 1] * valT::raw(m + 1 - i);
  for (int i = n - 2; ~i; --i) BP[i] = BP[i + 1] * valT::raw(m - 1 - i);
  valT ans = 0;
  for (int i = 0; i < n; ++i) {
    valT x = f[i] * AP[i] * BP[i] * B.ifac_[i] * B.ifac_[n - 1 - i];
    ans += (n - 1 - i) & 1 ? -x : x;
  }
  return ans;
}
// Lagrange theorem $f(x) =  \sum_{i = 0}^{n - 1} f_i \prod_{j \neq i} \frac{x - j}{i - j}$
// Simplies $f(m) = \sum_{i = 0}^{n - 1} (-1)^{n - 1 - i} f_i \binom{m}{i} \binom{m - i - 1}{n - 1 - i}$

// Calculate powSum in $O(k)$ based on Lagrange interpolation
template<typename valT, typename enable = ModT<valT>>
valT powSum(int n, int k, const std::vector<int>& sp) {
  if (k == 0) return valT(n);
  std::vector<valT> f(k + 2);
  f[1] = valT(1);
  for (int i = 2, nf = (int)f.size(); i < nf; ++i) {
    if (sp[i] == i) f[i] = pow(valT(i), k);
    else f[i] = f[sp[i]] * f[i / sp[i]];
  }
  for (int i = 1, nf = (int)f.size(); i < nf; ++i) f[i] += f[i - 1];
  return Lagrange(f, n);
}
// https://codeforces.com/problemset/problem/622/F


template<typename valT, typename enable = ModT<valT>>
class Matrix {
  static inline constexpr int N = 1003;
  int n_;
 public:
  valT a_[N][N];
  Matrix() {}
  Matrix(int n, valT x = 0): n_(n) {
    all(0);
    for (int i = 0; i < n_; ++i) {
      a_[i][i] = x;
    }
  }
  void all(valT x) {
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        a_[i][j] = x;
      }
    }
  }
  Matrix& operator+=(const Matrix& rhs) {
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        a_[i][j] += rhs.a_[i][j];
      }
    }
    return (*this);
  }
  Matrix operator+(const Matrix& rhs) const {
    return Matrix(this) += rhs;
  }
  Matrix& operator-=(const Matrix& rhs) {
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        a_[i][j] -= rhs.a_[i][j];
      }
    }
    return (*this);
  }
  Matrix operator-(const Matrix& rhs) const {
    return Matrix(this) -= rhs;
  }
  Matrix operator*(const Matrix& rhs) const {
    Matrix R(n_);
    for (int i = 0; i < n_; ++i) {
      for (int k = 0; k < n_; ++k) {
        for (int j = 0; j < n_; ++j) {
          R.a_[i][j] += a_[i][k] * rhs.a_[k][j];
        }
      }
    }
    return R;
  }
  Matrix operator*=(const Matrix& rhs) {
    return (*this) = (*this) * rhs;
  }
  void print() {
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        std::cout << a_[i][j] << ' ';
      }
      std::cout << '\n';
    }
  }
  friend Matrix pow(Matrix A, int n) {
    Matrix R(A.n_, valT(1));
    while (n) {
      if (n&1) R = R * A;
      n >>= 1; A = A * A;
    }
    return R;
  }
};

// You must change mid to be random one or you will be hack
template<typename T> // don't use it
void quickSort(std::vector<T>& a) {
  std::function<void(int, int)> qSort = [&](int l, int r) {
    int i = l, j = r;
    auto mid = a[(l + r) / 2];
    while (i <= j) {
      while (i <= j && a[i] < mid) ++i;
      while (j >= i && a[j] > mid) --j;
      if (i <= j) {
        std::swap(a[i], a[j]);
        ++i;
        --j;
      }
    }
    if (i < r) qSort(i, r);
    if (l < j) qSort(l, j);
  };
  qSort(0, a.size() - 1);
}

class MexS {
  static inline const int B = 64; // submit use 64bit
  using ULL = unsigned long long;
  std::vector<LL> cnt_;
  std::vector<std::vector<ULL>> a_;
  int ans_;
 public:
  // the answer is at most n
  MexS(int n) : cnt_(n + 1), ans_(-1) {
    int x = cnt_.size();
    while (x > B) {
      a_.emplace_back(std::vector<ULL>((x + B - 1) / B, -1ULL));
      x /= B;
    }
    a_.emplace_back(std::vector<ULL>{-1ULL});
  }
  void insert(int id) {
    if (id < 0 || id >= cnt_.size() || cnt_[id]++) return;
    if (id == ans_) ans_ = -1;
    for (auto &x : a_) {
      int tid = id / B;
      x[tid] ^= 1ULL << id - tid * B;
      if (x[tid]) return;
      id = tid;
    }
  }
  void erase(int id) { // make sure there is an element in this set
    if (id < 0 || id >= cnt_.size() || --cnt_[id]) return;
    if (id <= ans_) ans_ = id;
    for (auto &x : a_) {
      int tid = id / B;
      x[tid] ^= 1ULL << id - tid * B;
      if (x[tid] != -1ULL) return;
      id = tid;
    }
  }
  int solve() {
    if (ans_ == -1) {
      ans_ = 0;
      for (auto it = a_.crbegin(); it != a_.crend(); ++it) {
        ans_ = ans_ * B + __builtin_ctzll((*it)[ans_]);
      }
    }
    return ans_;
  }
};

class MEX {
  // B may need to be bigger
  static inline constexpr int B = 20;
  std::array<std::map<int, int>, B> mp;
  std::map<int, int> S;
 public:
  void insert(int x) {
    if (S[x]++) return;
    int mask = 0;
    for (int i = B - 1; i >= 0; --i) {
      mask |= 1 << i;
      ++mp[i][x & mask];
    }
  }
  void erase(int x) {
    if (--S[x]) return;
    S.erase(x);
    int mask = 0;
    for (int i = B - 1; i >= 0; --i) {
      mask |= 1 << i;
      --mp[i][x & mask];
    }
  }
  // find mex(a_i ^ x)
  int solve(int x = 0) {
    int mask = 0, r = 0;
    for (int i = B - 1; i >= 0; --i) {
      mask |= x & (1 << i);
      if (mp[i][mask] == (1 << i)) {
        mask ^= 1 << i;
        r |= 1 << i;
      }
    }
    return r;
  }
};


// transform vector<int> to vector<valT>
template<typename valT, typename enable = ModT<valT>>
std::vector<valT> trans(const std::vector<int>& a) {
  int n = (int)a.size();
  std::vector<valT> ans(n);
  for (int i = 0; i < n; ++i) ans[i] = valT(a[i]);
  return ans;
}

// Shortest recursive relational formula: https://cmwqf.github.io/2020/07/18/%E6%B5%85%E8%B0%88Berlekamp-Massey%E7%AE%97%E6%B3%95/
template<typename valT, typename enable = ModT<valT>>
static std::vector<valT> BerlekampMassey(const std::vector<valT>& a) {
  std::vector<valT> ans, lst;
  valT delta = 0;
  for (int i = 0, w = -1, n = (int)a.size(); i < n; ++i) {
    valT t = 0;
    for (int j = 0, na = ans.size(); j < na; ++j) {
      t += ans[j] * a[i - 1 - j];
    }
    if (t == a[i]) continue;
    // first time ans fail
    if (w == -1) {
      w = i; delta = a[i];
      ans.emplace_back(0);
      continue;
    }
    auto now = ans;
    auto mul = (a[i] - t) / delta;
    if (i - w + lst.size() > ans.size()) ans.resize(i - w + lst.size());
    ans[i - w - 1] += mul;
    for (int j = 0, lj = lst.size(); j < lj; ++j) ans[i - w + j] -= mul * lst[j];
    if ((int)now.size() - i < (int)lst.size() - w) {
      w = i; delta = a[i] - t; std::swap(now, lst);
    }
  }
  return ans;
}
// https://www.luogu.com.cn/problem/P5487
