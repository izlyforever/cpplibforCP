#pragma once
#include <bits/stdc++.h>
#include "../template.hpp"

// Returns the original value corresponding to the array value after discretization
template<typename T, typename enable = IntegerT<T>>
std::vector<T> discrete(std::vector<T>& a) {
  auto r = a;
  std::sort(r.begin(), r.end());
  r.erase(std::unique(r.begin(), r.end()), r.end());
  for (auto& x : a) {
    x = std::lower_bound(r.begin(), r.end(), x) - r.begin();
  }
  return r;
}

// make sure that a[i].first <= a[i].second
void disjointInterval(std::vector<std::pair<int, int>>& a) {
  if (a.size() <= 1) return;
  std::vector<std::pair<int, int>> b;
  std::sort(a.begin(), a.end());
  int l = a[0].first, r = a[0].second;
  for (int i = 1, n_ = (int)a.size(); i < n_; ++i) {
    if (a[i].first <= r) {
      r = std::max(r, a[i].second);
    } else {
      b.emplace_back(l, r);
      l = a[i].first, r = a[i].second;
    }
  }
  b.emplace_back(l, r);
  std::swap(a, b);
}

template<typename T>
class Internal {
  std::set<std::pair<T, T>> S;
 public:
  void clear() {
    S.clear();
  }
  int size() const {
    return S.size();
  }
  bool contain(T x) {
    return contain(x, x + 1);
  }
  bool contain(T l, T r) {
    if (l >= r) return false;
    auto it = S.lower_bound({l + 1, l});
    return it != S.begin() && std::prev(it)->second >= r;
  }
  bool intersect(T l, T r) {
    if (l >= r) return false;
    auto it = S.lower_bound({l + 1, l});
    return (it != S.end() && it->first < r) || (it != S.begin() && std::prev(it)->second > l);
  }
  void add(T l, T r) {
    if (l >= r) return;
    auto it = S.lower_bound({l + 1, l});
    if (it != S.begin()) {
      auto pit = std::prev(it);
      if (pit->second >= l) {
        l = pit->first;
        r = std::max(r, pit->second);
        S.erase(pit);
      }
    }
    while (it != S.end() && it->first <= r) {
      r = std::max(r, it->second);
      it = S.erase(it);
    }
    S.emplace(l, r);
  }
  void sub(T l, T r) {
    if (l >= r) return;
    auto it = S.lower_bound({l + 1, l});
    if (it != S.begin()) {
      auto pit = std::prev(it);
      if (pit->second > l) {
        int pl = pit->first;
        int pr = pit->second;
        S.erase(pit);
        if (pl < l) S.insert(pl, l);
        if (r > pr) S.insert(r, pr);
      }
    }
    while (it != S.end() && it->first < r) {
      if (it->second > r) S.insert(r, it->second);
      it = S.erase(it);
    }
  }
};

template<typename T, typename enable = IntegerT<T>>
class RingBuffer {
  int m_, id_;
  std::vector<T> a_;
 public:
  RingBuffer(int m) : m_(m), id_(0), a_(m, -1) {};
  T getCurrent() const {
    return a_[id_];
  }
  void insert(T x) {
    a_[id_++] = x;
    if (id_ == m_) id_ = 0;
  }
};
// https://codeforces.com/gym/103274/problem/G
