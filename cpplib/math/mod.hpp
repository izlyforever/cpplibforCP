#pragma once
#include <bits/stdc++.h>
using LL = long long;

// You should setMod before use it
class ModInt {
  static inline int M = 998244353;
  int n_;
// |x| <= max(1, b), |y| <= a ===> |y - a / b * x| <= a % b + a / b * b = a 
  static std::tuple<int, int, int> exGcdInternal(int a, int b) {
    if (b == 0) return {a, 1, 0};
    auto [d, y, x] = exGcdInternal(b, a % b);
    return {d, x, y - a / b * x};
  }
  static int inv(int a) {
    auto [d, x, y] = exGcdInternal(a, M);
    // assert(d == 1);
    return x < 0 ? x + M : x;
  }
  // assum M is prime
  static int invP(int x) {
    return x == 1 ? x : 1LL * (M - M / x) * invP(M % x) % M;
  }
 public:
  template<typename T>
   operator T() const {
    return static_cast<T>(n_);
  }
  static void setMod(int m) {
    assert(M == m);
  }
  static int mod() {
    return M;
  }
  // assume 0 <= x < M
  static ModInt raw(int x) {
    ModInt A;
    A.n_ = x;
    return A;
  }
  ModInt() { n_ = 0;}
  ModInt(const int& x) : n_(x % M) {
    if (n_ < 0) n_ += M;
  }
  ModInt(const LL& x) : n_(x % M) {
    if (n_ < 0) n_ += M;
  }
  ModInt operator-() const {
    return n_ == 0 ? *this : raw(M - n_);
  }
  ModInt& operator++() {
    if (++n_ == M) n_ = 0;
    return *this;
  }
  ModInt& operator--() {
    if (n_-- == 0) n_ += M;
    return *this;
  }
  ModInt& operator+=(const ModInt& A) {
    n_ += A.n_;
    if (n_ >= M) n_ -= M;
    return (*this);
  }
  ModInt& operator-=(const ModInt& A) {
    n_ -= A.n_;
    if (n_ < 0) n_ += M;
    return (*this);
  }
  ModInt& operator*=(const ModInt& A) {
    n_ = 1LL * n_ * A.n_ % M;
    return (*this);
  }
  ModInt& operator/=(const ModInt& A) {
    return (*this) *= A.inv();
  }
  ModInt operator+(const ModInt& A) const {
    return ModInt(*this) += A;
  }
  ModInt operator-(const ModInt& A) const {
    return ModInt(*this) -= A;
  }
  ModInt operator*(const ModInt& A) const {
    return ModInt(*this) *= A;
  }
  ModInt operator/(const ModInt& A) const {
    return ModInt(*this) /= A;
  }
  ModInt operator<<(int x) const {
    static constexpr int bits = 32;
    LL r = n_;
    while (x > bits) {
      x -= bits;
      r <<= bits;
      r %= M;
    }
    r <<= x;
    return ModInt(r);
  }
  ModInt& operator<<=(int x) {
    return (*this) = (*this) << x;
  }
  bool operator==(const ModInt& A) const {
    return n_ == A.n_;
  }
  bool operator!=(const ModInt& A) const {
    return n_ != A.n_;
  }
  ModInt inv() const {
    return inv(n_);
  }
  ModInt invP() const {
    return invP(n_);
  }
  friend ModInt pow(ModInt A, int n) {
    ModInt R(1);
    while (n) {
      if (n & 1) R *= A;
      n >>= 1;   A *= A;
    }
    return R;
  }
  friend std::istream& operator>>(std::istream& in, ModInt& A) {
    LL x;
    in >> x;
    A = ModInt(x);
    return in;
  }
  friend std::ostream& operator<<(std::ostream& out, const ModInt& A) {
    out << A.n_;
    return out;
  }
};

// You should setMod before use it
class ModLL {
  static inline LL M = 998244353;
  LL n_;
// |x| <= max(1, b), |y| <= a ===> |y - a / b * x| <= a % b + a / b * b = a 
  static std::tuple<LL, LL, LL> exGcdInternal(LL a, LL b) {
    if (b == 0) return {a, 1, 0};
    auto [d, y, x] = exGcdInternal(b, a % b);
    return {d, x, y - a / b * x};
  }
  static LL inv(LL a) {
    auto [d, x, y] = exGcdInternal(a, M);
    // assert(d == 1);
    return x < 0 ? x + M : x;
  }
  // assum M is prime
  static LL invP(LL x) {
    return x == 1 ? x : __int128_t(M - M / x) * invP(M % x) % M;
  }
 public:
  template<typename T>
   operator T() const {
    return static_cast<T>(n_);
  }
  static void setMod(LL m) {
    M = m;
  }
  static LL mod() {
    return M;
  }
  // assume 0 <= x < M
  static ModLL raw(LL x) {
    ModLL A;
    A.n_ = x;
    return A;
  }
  ModLL() { n_ = 0;}
  ModLL(const int& x) : n_(x % M) {
    if (n_ < 0) n_ += M;
  }
  ModLL(const LL& x) : n_(x % M) {
    if (n_ < 0) n_ += M;
  }
  ModLL(const __int128_t& x) : n_(x % M) {
    if (n_ < 0) n_ += M;
  }
  ModLL operator-() const {
    return n_ == 0 ? *this : raw(M - n_);
  }
  ModLL& operator++() {
    if (++n_ == M) n_ = 0;
    return *this;
  }
  ModLL& operator--() {
    if (n_-- == 0) n_ += M;
    return *this;
  }
  ModLL& operator+=(const ModLL& A) {
    n_ += A.n_;
    if (n_ >= M) n_ -= M;
    return (*this);
  }
  ModLL& operator-=(const ModLL& A) {
    n_ -= A.n_;
    if (n_ < 0) n_ += M;
    return (*this);
  }
  ModLL& operator*=(const ModLL& A) {
    n_ = __int128_t(n_) * A.n_ % M;
    return (*this);
  }
  ModLL& operator/=(const ModLL& A) {
    return (*this) *= A.inv();
  }
  ModLL operator+(const ModLL& A) const {
    return ModLL(*this) += A;
  }
  ModLL operator-(const ModLL& A) const {
    return ModLL(*this) -= A;
  }
  ModLL operator*(const ModLL& A) const {
    return ModLL(*this) *= A;
  }
  ModLL operator/(const ModLL& A) const {
    return ModLL(*this) /= A;
  }
  ModLL operator<<(int x) const {
    static constexpr int bits = 64;
    __int128_t r = n_;
    while (x > bits) {
      x -= bits;
      r <<= bits;
      r %= M;
    }
    r <<= x;
    return ModLL(r);
  }
  ModLL& operator<<=(int x) {
    return (*this) = (*this) << x;
  }
  bool operator==(const ModLL& A) const {
    return n_ == A.n_;
  }
  bool operator!=(const ModLL& A) const {
    return n_ != A.n_;
  }
  ModLL inv() const {
    return inv(n_);
  }
  ModLL invP() const {
    return invP(n_);
  }
  friend ModLL pow(ModLL A, LL n) {
    ModLL R(1);
    while (n) {
      if (n & 1) R *= A;
      n >>= 1;   A *= A;
    }
    return R;
  }
  friend std::istream& operator>>(std::istream& in, ModLL& A) {
    LL x;
    in >> x;
    A = ModLL(x);
    return in;
  }
  friend std::ostream& operator<<(std::ostream& out, const ModLL& A) {
    out << A.n_;
    return out;
  }
};


template<int N>
class MInt {
  static inline constexpr int M = N;
  int n_;
  // |x| <= max(1, b), |y| <= a ===> |y - a / b * x| <= a % b + a / b * b = a 
  static std::tuple<int, int, int> exGcdInternal(int a, int b) {
    if (b == 0) return {a, 1, 0};
    auto [d, y, x] = exGcdInternal(b, a % b);
    return {d, x, y - a / b * x};
  }
  static int inv(int a) {
    auto [d, x, y] = exGcdInternal(a, M);
    // assert(d == 1);
    return x < 0 ? x + M : x;
  }
  // assum M is prime
  static int invP(int x) {
    return x == 1 ? x : 1LL * (M - M / x) * invP(M % x) % M;
  }
 public:
  template<typename T>
   operator T() const {
    return static_cast<T>(n_);
  }
  static void setMod(int m) {
    assert(M == m);
  }
  static constexpr int mod() {
    return M;
  }
  // assume 0 <= x < M
  static MInt raw(int x) {
    MInt A;
    A.n_ = x;
    return A;
  }
  MInt() { n_ = 0;}
  MInt(const int& x) : n_(x % M) {
    if (n_ < 0) n_ += M;
  }
  MInt(const LL& x) : n_(x % M) {
    if (n_ < 0) n_ += M;
  }
  MInt operator-() const {
    return n_ == 0 ? *this : raw(M - n_);
  }
  MInt& operator++() {
    if (++n_ == M) n_ = 0;
    return *this;
  }
  MInt& operator--() {
    if (n_-- == 0) n_ += M;
    return *this;
  }
  MInt& operator+=(const MInt& A) {
    n_ += A.n_;
    if (n_ >= M) n_ -= M;
    return (*this);
  }
  MInt& operator-=(const MInt& A) {
    n_ -= A.n_;
    if (n_ < 0) n_ += M;
    return (*this);
  }
  MInt& operator*=(const MInt& A) {
    n_ = 1LL * n_ * A.n_ % M;
    return (*this);
  }
  MInt& operator/=(const MInt& A) {
    return (*this) *= A.inv();
  }
  MInt operator+(const MInt& A) const {
    return MInt(*this) += A;
  }
  MInt operator-(const MInt& A) const {
    return MInt(*this) -= A;
  }
  MInt operator*(const MInt& A) const {
    return MInt(*this) *= A;
  }
  MInt operator/(const MInt& A) const {
    return MInt(*this) /= A;
  }
  MInt operator<<(int x) const {
    static constexpr int bits = 32;
    LL r = n_;
    while (x > bits) {
      x -= bits;
      r <<= bits;
      r %= M;
    }
    r <<= x;
    return MInt(r);
  }
  MInt& operator<<=(int x) {
    return (*this) = (*this) << x;
  }
  bool operator==(const MInt& A) const {
    return n_ == A.n_;
  }
  bool operator!=(const MInt& A) const {
    return n_ != A.n_;
  }
  MInt inv() const {
    return inv(n_);
  }
  MInt invP() const {
    return invP(n_);
  }
  friend MInt pow(MInt A, int n) {
    MInt R(1);
    while (n) {
      if (n & 1) R *= A;
      n >>= 1;   A *= A;
    }
    return R;
  }
  friend std::istream& operator>>(std::istream& in, MInt& A) {
    LL x;
    in >> x;
    A = MInt(x);
    return in;
  }
  friend std::ostream& operator<<(std::ostream& out, const MInt& A) {
    out << A.n_;
    return out;
  }
};

// valT
template<class T>
struct is_MInt : std::false_type {};

template<int M>
struct is_MInt<MInt<M>> : std::true_type {};

template<class T>
inline constexpr bool is_mint_v = is_MInt<T>::value;

template<typename T>
using ModT = std::enable_if_t<std::is_same_v<ModLL, T> || std::is_same_v<ModInt, T> || is_mint_v<T>>;
