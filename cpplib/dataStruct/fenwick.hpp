#pragma once
#include <bits/stdc++.h>
#include "../template.hpp"

// Bit Tree Mininal version
template<typename T, typename enable = IntegerT<T>>
struct BitreeMin {
  std::vector<T> s_;
  BitreeMin() {}
  BitreeMin(int n) : s_(n + 1, std::numeric_limits<T>::max()) {}
  int lowbit(int n) { return n & (-n); }
  void modify(int id, T p) {
    int ns = s_.size();
    while (id < ns) {
      s_[id] = std::min(s_[id], p);
      id += lowbit(id);
    }
  }
  // cal minial value in [1, id]
  T min(int id) {
    T r = std::numeric_limits<T>::max();
    while (id) {
      r = std::min(r, s_[id]);
      id -= lowbit(id);
    }
    return r;
  }
};

template<typename T, typename enable = IntegerT<T>>
struct Bitree {
  std::vector<T> s_;
  Bitree() {}
  Bitree(int n) : s_(n + 1) {}
  int lowbit(int n) { return n & (-n); }
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