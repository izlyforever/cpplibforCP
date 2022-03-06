#pragma once
#include <bits/stdc++.h>

// Second Block abs version(online)
class BlockAbs {
  int l_, r_; // fa_ \in [l_, r_]
  int f_, d_; // x \in [l_, r_] has real value f_ x - d_, where f_ = 1, -1
  std::vector<int> fa_;
  int find(int x) {
    return x == fa_[x] ? x : fa_[x] = find(fa_[x]);
  }
  void merge(int x, int y) { // merge x to y
    fa_[find(x)] = find(y);
  }
 public:
  BlockAbs(int mx) : l_(0), r_(mx), f_(1), d_(0), fa_(mx + 1) {
    std::iota(fa_.begin(), fa_.end(), 0);
  }
  void add(int x) { // |fi - d_ - x| = | i - (x + d_) f_|
    x = (x + d_) * f_;
    if ((l_ + r_) < 2 * x) {
      f_ = -1;
      d_ = -x;
      if (x < r_) {
        for (int i = r_; i > x; --i) merge(i, 2 * x - i);
        r_ = x;
      }
    } else {
      f_ = 1;
      d_ = x;
      if (x > l_) {
        for (int i = l_; i < x; ++i) merge(i, 2 * x - i);
        l_ = x;
      }
    }
  }
  int query(int x) {
    return find(x) * f_ - d_;
  }
};
// https://codeforces.com/gym/103104/problem/K

// Second Block minus version(online, space optim can be done if offline)
class BlockMinus {
  std::vector<int> fa_, sz_, a_;
  int l_, delta_, mx_; // real value x - delta_
  int find(int x) {
    return x == fa_[x] ? x : fa_[x] = find(fa_[x]);
  }
  void merge(int x, int y) { // merge x to y
    x = find(x); y = find(y);
    if (x == y) return;
    fa_[x] = y;
    sz_[y] += sz_[x];
    sz_[x] = 0;
  }
  void modifyPart(int ql, int qr, int x) {
    for (int i = ql; i < qr; ++i) {
      a_[i] = find(a_[i]);
      if (a_[i] - delta_ > x) {
        --sz_[a_[i]];
        a_[i] = find(a_[i] - x);
        ++sz_[a_[i]];
      }
    }
  }
  void modifyAll(int x) {
    if (x < mx_ - delta_ - x) {
      for (int i = delta_ + 1; i <= x + delta_; ++i) merge(i, i + x);
      delta_ += x;
    } else {
      for (int i = mx_; i > x + delta_; --i) merge(i, i - x);
      mx_ = x + delta_;
    }
  }
  int queryPart(int ql, int qr, int x) {
    int ans = 0;
    for (int i = ql; i < qr; ++i) {
      if (find(a_[i]) - delta_ == x) ++ans;
    }
    return ans;
  }
  int queryAll(int x) {
    x += delta_;
    if (x > mx_ || find(x) != x) return 0;
    return sz_[x];
  }
 public:
  void init(const std::vector<int>& a, int l, int r) {
    l_ = l, delta_ = 0;
    a_ = {a.begin() + l, a.begin() + r};
    mx_ = *std::max_element(a_.begin(), a_.end());
    fa_.resize(mx_ + 1); std::iota(fa_.begin(), fa_.end(), 0);
    sz_.resize(mx_ + 1); for (auto x : a_) ++sz_[x];
  }
  void modify(int ql, int qr, int x) {
    if (x >= mx_ - delta_) return;
    if (qr - ql == (int)a_.size()) modifyAll(x);
    else modifyPart(ql - l_, qr - l_, x);
  }
  int query(int ql, int qr, int x) {
    if (qr - ql == (int)a_.size()) return queryAll(x);
    return queryPart(ql - l_, qr - l_, x);
  }
};
// https://codeforces.com/contest/896/problem/E
