#pragma once
#include <bits/stdc++.h>
#include "../template.hpp"

// Bit Tree Maximal version
template<typename T>
struct BitreeMax {
  std::vector<T> s_;
  BitreeMax() {}
  BitreeMax(int n) : s_(n + 1, std::numeric_limits<T>::min()) {}
  static int lowbit(int n) { return n & (-n); }
  void modify(int id, T p) {
    int ns = s_.size();
    while (id < ns) {
      s_[id] = std::max(s_[id], p);
      id += lowbit(id);
    }
  }
  // cal maximal value in [1, id]
  T max(int id) {
    T r = std::numeric_limits<T>::min();
    while (id) {
      r = std::max(r, s_[id]);
      id -= lowbit(id);
    }
    return r;
  }
};

// Bit Tree Minimal version
template<typename T, typename enable = IntegerT<T>>
struct BitreeMin {
  std::vector<T> s_;
  BitreeMin() {}
  BitreeMin(int n) : s_(n + 1, std::numeric_limits<T>::max()) {}
  static int lowbit(int n) { return n & (-n); }
  void modify(int id, T p) {
    int ns = s_.size();
    while (id < ns) {
      s_[id] = std::min(s_[id], p);
      id += lowbit(id);
    }
  }
  // cal minimal value in [1, id]
  T min(int id) {
    T r = std::numeric_limits<T>::max();
    while (id) {
      r = std::min(r, s_[id]);
      id -= lowbit(id);
    }
    return r;
  }
};
// see https://codeforces.com/contest/1635/submission/147077087 for more elegent impl

template<typename T, typename enable = IntegerT<T>>
struct Bitree {
  std::vector<T> s_;
  Bitree() {}
  Bitree(int n) : s_(n + 1) {}
  // https://codeforces.com/blog/entry/59305
  Bitree(const std::vector<T>& a) : s_(a.size() + 1) {
    int n = a.size();
    for (int i = 0; i < n; ++i) s_[i + 1] = s_[i] + a[i];
    for (int i = n; i > 0; --i) s_[i] -= s_[i - lowbit(i)];
  }
  std::vector<T> getOrigin() const {
    auto a = s_;
    int n = s_.size() - 1;
    for (int i = 1; i <= n; ++i) a[i] += a[i - lowbit(i)];
    std::vector<T> ans(n);
    for (int i = n - 1; i >= 0; --i) ans[i] = a[i + 1] - a[i];
    return ans;
  }
  static int lowbit(int n) { return n & (-n); }
  void add(int id, T p) {
    int ns = s_.size();
    while (id < ns) {
      s_[id] += p;
      id += lowbit(id);
    }
  }
  T sum(int id) {
    T r = 0;
    while (id) {
      r += s_[id];
      id -= lowbit(id);
    }
    return r;
  }
  T sum(int l, int r) { return sum(r) - sum(l - 1); }
  T at(int id) { return sum(id, id);}
  // val must monic increase
  int lower_bound(T x) {
    int l = 1, r = s_.size() - 1;
    while (l <= r) {
      int m = (l + r) / 2;
      if (at(m) >= x) r = m - 1;
      else l = m + 1;
    }
    return l;
  }
  // val must monic increase
  int upper_bound(T x) {
    int l = 1, r = s_.size() - 1;
    while (l <= r) {
      int m = (l + r) / 2;
      if (at(m) > x) r = m - 1;
      else l = m + 1;
    }
    return l;
  }
  // find minimal index s.t. sum(id) >= x, sum must be increased
  int search(T val) {
    T sum = 0;
    int id = 0;
    for (int i = std::__lg(s_.size()); ~i; --i) {
      if (id + (1 << i) < (int)s_.size() && sum + s_[id + (1 << i)] < val) {
        id += (1 << i);
        sum += s_[id];
      }
    }
    return ++id;
  }
};

template<typename T, typename enable = IntegerT<T>>
class BitreePlus {
  int n_;
  // c[i] = a[i] - a[i - 1], b_i = (i - 1) * c_i
  Bitree<T> B, C;
  void add(int id, T p) {
    C.add(id, p);
    B.add(id, (id - 1) * p);
  }
 public:
  BitreePlus() {}
  BitreePlus(int n) : n_(n), B(n), C(n) {}
  void add(int l, int r, T p) {
    add(l, p);
    add(r + 1, -p);
  }
  T sum(int id) { return id * C.sum(id) - B.sum(id); }
  T sum(int l, int r) { return sum(r) - sum(l - 1); }
  T at(int id) { return sum(id, id); }
  // val must monic increase
  int lower_bound(T x) {
    int l = 1, r = n_;
    while (l <= r) {
      int m = (l + r) / 2;
      if (at(m) >= x) r = m - 1;
      else l = m + 1;
    }
    return l;
  }
  // val must monic increase
  int upper_bound(T x) {
    int l = 1, r = n_;
    while (l <= r) {
      int m = (l + r) / 2;
      if (at(m) > x) r = m - 1;
      else l = m + 1;
    }
    return l;
  }
  // find minimal index s.t. sum(id) >= x, sum must be increased
  int search(T val) {
    T sumB = 0, sumC = 0;
    int id = 0;
    for (int i = std::__lg(n_); ~i; --i)
      if (int idi = id + (1 << i); idi <= n_) {
        if (idi * (sumC + C.s[idi]) - B.s[idi] - sumB < val) {
          id = idi;
          sumB += B.s[id];
          sumC += C.s[id];
        }
      }
    return ++id;
  }
};

struct Bitree2M {
  static inline int lowbit(int n) { return n & (-n); }
  Bitree2M() {}
  Bitree2M(int n, int m) : n_(n), m_(m), mp_(n + 1) {}
  void add(int x, int y, int p) {
    for (int i = x; i <= n_; i += lowbit(i)) {
      auto &a = mp_[i];
      for (int j = y; j <= m_; j += lowbit(j)) {
        a[j] += p;
      }
    }
  }
  int sum(int x, int y) {
    x = std::min(x, n_);
    y = std::min(y, m_);
    int ans = 0;
    for (int i = x; i > 0; i -= lowbit(i)) {
      auto &a = mp_[i];
      if (a.empty()) continue;
      for (int j = y; j > 0; j -= lowbit(j)) {
        if (a.count(j)) {
          ans += a[j];
        }
      }
    }
    return ans;
  }
  int sum(int x1, int y1, int x2, int y2) {
    --x1; --y1;
    return sum(x2, y2) - sum(x2, y1) - sum(x1, y2) + sum(x1, y1);
  }
  int n_, m_;
  std::vector<std::map<int, int>> mp_;
};

struct Bitree2V {
  static inline int lowbit(int n) { return n & (-n); }
  Bitree2V() {}
  Bitree2V(int n, int m) : n_(n), m_(m), mp_(n + 1, std::vector<int>(m + 1)) {}
  void add(int x, int y, int p) {
    for (int i = x; i <= n_; i += lowbit(i)) {
      auto &a = mp_[i];
      for (int j = y; j <= m_; j += lowbit(j)) {
        a[j] += p;
      }
    }
  }
  int sum(int x, int y) {
    x = std::min(x, n_);
    y = std::min(y, m_);
    int ans = 0;
    for (int i = x; i > 0; i -= lowbit(i)) {
      auto &a = mp_[i];
      for (int j = y; j > 0; j -= lowbit(j)) {
        ans += a[j];
      }
    }
    return ans;
  }
  int sum(int x1, int y1, int x2, int y2) {
    --x1; --y1;
    return sum(x2, y2) - sum(x2, y1) - sum(x1, y2) + sum(x1, y1);
  }
  int n_, m_;
  std::vector<std::vector<int>> mp_;
};
