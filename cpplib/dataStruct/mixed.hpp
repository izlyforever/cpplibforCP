#pragma once
#include <bits/stdc++.h>
#include "basic.hpp"
#include "fenwick.hpp"

using LL = long long;

// a will becomes next lexicographical order of a, satisfies $-1 < a_0 < a_1 < \cdots, a_{n - 1} < mx$
bool nextBinom(std::vector<int>& a, int mx) {
  int n = (int)a.size(), i = 1;
  while (i <= n && a[n - i] == mx - i) ++i;
  if (i > n) return false;
  ++a[n - i];
  for (int j = 1; j < i; ++j) {
    a[n - j] = a[n - i] + (i - j);
  }
  return true;
}

// total number binom{mx}{n}
void bruteForceBinom(int n, int mx) {
  std::vector<int> a(n);
  std::iota(a.begin(), a.end(), 0);
  do {
    // do something
    for (auto x : a) std::cout << x << ' ';
    std::cout << '\n';
  } while (nextBinom(a, mx));
}

// Error Correction Code: O(n_ m_ + k_^k_ n_)
class ECC {
  std::vector<std::vector<int>> a_;   // origin data: n rows, m cols.
  int k_;                             // Maximum number of differences allowed
  std::vector<std::vector<int>> bad_; // difference with current answer
  int n_, m_, mxId_;
  void updateMxId(int i) {
    if (bad_[i].size() > bad_[mxId_].size()) mxId_ = i;
  }
  bool dfs(int c) {  // remain time that current answer can change
    auto bd = bad_[mxId_];
    if ((int)bd.size() <= k_) return true;
    if ((int)bd.size() - k_ > c) return false;
    // Note that bd is O(k_) instead of O(m_)
    std::vector<int> f(bd.size() - k_);
    iota(f.begin(), f.end(), 0);
    int tMxId = mxId_;
    do {
      mxId_ = tMxId;
      std::queue<int> tmp;
      for (auto x : f) {
        tmp.push(r_[bd[x]]);
        for (int i = 0; i < n_; ++i) {
          if (a_[i][bd[x]] == r_[bd[x]]) {
            bad_[i].emplace_back(bd[x]);
          }
          if (a_[i][bd[x]] == a_[mxId_][bd[x]]) {
            bad_[i].erase(
              std::find(bad_[i].begin(), bad_[i].end(), bd[x]));
          }
        }
        r_[bd[x]] = a_[mxId_][bd[x]];
      }
      for (int i = 0; i < n_; ++i) updateMxId(i);
      if (dfs(c - f.size())) return true;
      for (auto x : f) {
        for (int i = 0; i < n_; ++i) {
          if (a_[i][bd[x]] == r_[bd[x]]) {
            bad_[i].emplace_back(bd[x]);
          }
          if (a_[i][bd[x]] == tmp.front()) {
            bad_[i].erase(
              std::find(bad_[i].begin(), bad_[i].end(), bd[x]));
          }
        }
        r_[bd[x]] = tmp.front();
        tmp.pop();
      }
    } while (nextBinom(f, bd.size()));
    return false;
  }
 public:
  std::vector<int> r_;   // m_ cols vector, current answer
  ECC(std::vector<std::vector<int>> a) : a_(a), r_(a[0]) {
    n_ = a_.size(); m_ = r_.size();
    bad_.resize(n_);
    mxId_ = 0;
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < m_; ++j)
        if (a_[i][j] != r_[j]) {
          bad_[i].emplace_back(j);
        }
      updateMxId(i);
    }
  }
  void setK(int k) { k_ = k; }
  bool solve() { return dfs(k_); }
};


// length of longest increasing subsquence
int LIS(std::vector<int>& a) {
  std::vector<int> b;
  for (auto x : a) {
    if (auto it = std::lower_bound(b.begin(), b.end(), x); it == b.end()) {
      b.emplace_back(x);
    } else {
      *it = x;
    }
  }
  return b.size();
}
// length of longest non-decreasing subsquence
int LNDS(std::vector<int>& a) {
  std::vector<int> b;
  for (auto x : a) {
    if (auto it = std::upper_bound(b.begin(), b.end(), x); it == b.end()) {
      b.emplace_back(x);
    } else {
      *it = x;
    }
  }
  return b.size();
}

// longest non-decreasing subsquence, return choosen index
template<typename T>
auto LISP(const std::vector<T>& a) {
  int n = a.size();
  std::vector<T> b;
  b.reserve(n);
  // p[i] means the preview number of choosen i
  std::vector<int> p(n);
  std::iota(p.begin(), p.end(), 0);
  // q[i] means current b[i] is a[q[i]]
  std::vector<int> q;
  q.reserve(n);
  for (int i = 0; i < n; ++i) {
    auto it = std::upper_bound(b.begin(), b.end(), a[i]);
    if (it == b.end()) {
      if (!q.empty()) p[i] = q.back();
      b.emplace_back(a[i]);
      q.emplace_back(i);
    } else {
      *it = a[i];
      auto t = it - b.begin();
      q[t] = i;
      if (t > 0) p[i] = q[t - 1];
    }
  }
  std::stack<int> c;
  int now = q.back();
  c.push(now);
  while (now != p[now]) {
    now = p[now];
    c.push(now);
  }
  return c;
}

// monicDeque: index of every max element of SubInterval of length m
std::vector<int> monicDequeMax(std::vector<int>& a, int m) {
  std::vector<int> r;
  std::deque<int> Q;
  for (int i = 0, na = (int)a.size(); i < na; ++i) {
    if (!Q.empty() && i - Q.front() >= m) Q.pop_front();
    // change > to < if you want monicDequeMin
    while (!Q.empty() && a[i] > a[Q.back()]) Q.pop_back();
    Q.push_back(i);
    if (i >= m - 1) r.emplace_back(Q.front());
  }
  return r;
}

// f is index of a such that $a_{f_0} < a_{f_1} < a_{f_m}$
std::vector<int> monicStack(const std::vector<int>& a) {
  int n = (int)a.size();
  std::vector<int> f(n);
  std::stack<int> S;
  for (int i = 0; i < n; ++i) {
    while (!S.empty() && a[S.top()] < a[i]) {
      f[S.top()] = i;
      S.pop();
    }
    S.push(i);
  }
  return f;
}
// https://www.luogu.com.cn/problem/P5788

// Cartesian Tree
struct cNode {
  int id, val, par, ch[2];
  void init(int _id, int _val, int _par) {
    id = _id, val = _val, par = _par, ch[0] = ch[1] = 0;
  }
};
int cartesian_build(std::vector<cNode>& tree, int n) {
  for (int i = 1; i <= n; ++i) {
    int k = i - 1;
    while (tree[k].val < tree[i].val) k = tree[k].par;
    tree[i].ch[0] = tree[k].ch[1];
    tree[k].ch[1] = i;
    tree[i].par = k;
    tree[tree[i].ch[0]].par = i;
  }
  return tree[0].ch[1];
}
// https://codeforces.com/contest/1490/problem/D


// we can only use the ideal of merge sort
LL inverseOrderCount(std::vector<int> a) {
  discrete(a);
  Bitree<int> A(*std::max_element(a.begin(), a.end()) + 1);
  LL ans = 0;
  for (int i = (int)a.size() - 1; i >= 0; --i) {
    ans += A.sum(a[i]);
    A.add(a[i] + 1, 1);
  }
  return ans;
}
// https://codeforces.com/contest/1602/problem/E
