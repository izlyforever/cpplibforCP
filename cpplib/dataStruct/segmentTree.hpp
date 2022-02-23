#pragma once
#include <bits/stdc++.h>
using LL = long long;

// all you need is provide data in Info and operator +
struct Info {
  friend Info operator+(const Info& A, const Info& B) {
    Info C;
    // impl
    return C;
  }
};

template<class Info>
class SegmentTreeInfo {
  int n_;
  std::vector<Info> info_;
  void pull(int p) {
    info_[p] = info_[p << 1] + info_[p << 1 | 1];
  }
  void modify(const Info& v, int x, int l, int r, int p) {
    if (r - l == 1) {
      info_[p] = v;
    } else {
      int m = (l + r) / 2;
      if (x < m) modify(v, x, l, m, p << 1);
      else modify(v, x, m, r, p << 1 | 1);
      pull(p);
    }
  }
 public:
  SegmentTreeInfo(int n) : n_(n), info_(4 << std::__lg(n)) {}
  SegmentTreeInfo(const std::vector<Info>& a) : SegmentTreeInfo(a.size()) {
    std::function<void(int, int, int)> build = [&](int l, int r, int p) {
      if (r - l == 1) {
        info_[p] = a[l];
      } else {
        int m = (l + r) / 2;
        build(l, m, p << 1);
        build(m, r, p << 1 | 1);
        pull(p);
      }
    };
    build(0, n_, 1);
  }
  void modify(const Info& v, int x) {
    modify(v, x, 0, n_, 1);
  }
};

// sum version is simple, the class blow is an example
class SegmentTree {
  int n_;
  std::vector<LL> sm, tag;
  void pull(int p) { sm[p] = sm[p << 1] + sm[p << 1 | 1]; }
  void pushTag(LL x, int l, int r, int p) {
    tag[p] += x;
    sm[p] += x * (r - l);
  }
  void push(int l, int r, int p) {
    if (tag[p]) {
      int m = (l + r) / 2;
      pushTag(tag[p], l, m, p << 1);
      pushTag(tag[p], m, r, p << 1 | 1);
      tag[p] = 0;
    }
  }
  void rangeAdd(int L, int R, LL x, int l, int r, int p) {
    if (L <= l && R >= r) {
      pushTag(x, l, r, p);
      return;
    }
    push(l, r, p);
    int m = (l + r) / 2;
    if (L < m) rangeAdd(L, R, x, l, m, p << 1);
    if (R > m) rangeAdd(L, R, x, m, r, p << 1 | 1);
    pull(p);
  }
  // you should implement it to meet for needs
  LL query(int L, int R, int l, int r, int p) {
    if (L <= l && R >= r) return sm[p];
    push(l, r, p);
    LL ans = 0;
    int m = (l + r) / 2;
    if (L < m) ans += query(L, R, l, m, p << 1);
    if (R > m) ans += query(L, R, m, r, p << 1 | 1);
    return ans;
  }
  void resize() {
    tag.resize(4 * n_);
    sm.resize(4 * n_);
  }
 public:
  SegmentTree(int n) : n_(n) { resize(); }
  SegmentTree(const std::vector<LL>& a) {
    n_ = (int)a.size();
    resize();
    std::function<void(int, int, int)> build = [&](int l, int r, int p) {
      if (r - l == 1) {
        sm[p] = a[l];
        return;
      }
      int m = (l + r) / 2;
      build(l, m, p << 1);
      build(m, r, p << 1 | 1);
      pull(p);
    };
    build(0, n_, 1);
  }
  void add(int L, int R, LL v) { rangeAdd(--L, R, v, 0, n_, 1); }
  LL query(int L, int R) { return query(--L, R, 0, n_, 1); }
};
// https://www.luogu.com.cn/problem/P3372

// There are two versions: sum and min/max
class SegmentTreeAddCountMin {
  int n_;
  std::vector<int> mx, tag;
  void pull(int p) {
    mx[p] = std::min(mx[p << 1], mx[p << 1 | 1]);
  }
  void pushTag(int x, int l, int r, int p) {
    tag[p] += x;
    mx[p] += x;
  }
  void push(int l, int r, int p) {
    if (tag[p]) {
      int m = (l + r) / 2;
      pushTag(tag[p], l, m, p << 1);
      pushTag(tag[p], m, r, p << 1 | 1);
      tag[p] = 0;
    }
  }
  void rangeAdd(int L, int R, int x, int l, int r, int p) {
    if (L <= l && R >= r) {
      pushTag(x, l, r, p);
      return;
    }
    push(l, r, p);
    int m = (l + r) / 2;
    if (L < m) rangeAdd(L, R, x, l, m, p << 1);
    if (R > m) rangeAdd(L, R, x, m, r, p << 1 | 1);
    pull(p);
  }
  // you should implement it to meet for needs
  int query(int L, int R, int l, int r, int p) {
    if (L <= l && R >= r) return mx[p];
    push(l, r, p);
    LL ans = 0;
    int m = (l + r) / 2;
    if (L < m) ans += query(L, R, l, m, p << 1);
    if (R > m) ans += query(L, R, m, r, p << 1 | 1);
    return ans;
  }
  void resize() {
    tag.resize(4 * n_);
    mx.resize(4 * n_);
  }
 public:
  SegmentTreeAddCountMin(int n) : n_(n) { resize(); }
  SegmentTreeAddCountMin(const std::vector<int>& a) {
    n_ = (int)a.size();
    resize();
    std::function<void(int, int, int)> build = [&](int l, int r, int p) {
      if (r - l == 1) {
        mx[p] = a[l];
        return;
      }
      int m = (l + r) / 2;
      build(l, m, p << 1);
      build(m, r, p << 1 | 1);
      pull(p);
    };
    build(0, n_, 1);
  }
  void add(int L, int R, int v) { rangeAdd(L, R, v, 0, n_, 1); }
  int query(int L, int R) { return query(L, R, 0, n_, 1); }
};
// min/max version slightly hard: you should record min/max, and second min/max value: https://codeforces.com/gym/102992/problem/J


// Persistent Segment Tree, reference: https://zhuanlan.zhihu.com/p/250565583
class PstSegTree {
  struct Node {
    int l, r;
    LL val;
  };
  int n_;
  std::vector<int> root_;  // version number
  std::vector<Node> tree_;
  void pushUp(int p) {
    tree_[p].val = tree_[tree_[p].l].val + tree_[tree_[p].r].val;
  }
  void update(int pos, int val, int l, int r, int p, int q) {
    tree_[q] = tree_[p];
    if (r - l == 1) {
      tree_[q].val = val;
    } else {
      int m = (l + r) / 2;
      if (pos < m)
        update(pos, val, l, m, tree_[p].l, tree_[q].l = newNode());
      else
        update(pos, val, m, r, tree_[p].r, tree_[q].r = newNode());
      pushUp(q);
    }
  }
  LL query(int L, int R, int l, int r, int p) {
    if (L <= l && R >= r) return tree_[p].val;
    int m = (l + r) / 2;
    LL ans = 0;
    if (L < m) ans += query(L, R, l, m, tree_[p].l);
    if (R > m) ans += query(L, R, m, r, tree_[p].r);
    return ans;
  }
 public:
  int newNode() {
    int sz_ = tree_.size();
    tree_.emplace_back(Node());
    return sz_;
  }
  PstSegTree(const std::vector<int>& a) : n_(a.size()) {
    root_.emplace_back(newNode());
    std::function<void(int, int, int)> build = [&](int l, int r, int p) {
      if (r - l == 1) {
        tree_[p].val = a[l];
      } else {
        int m = (l + r) / 2;
        build(l, m, tree_[p].l = newNode());
        build(m, r, tree_[p].r = newNode());
        pushUp(p);
      }
    };
    build(0, n_, root_.back());
  }
  // single point update, p is current version, q is new version
  void update(int pos, int val, int p, int q) {
    update(pos, val, 0, n_, p, q);
  }
  // segment query, p is current version, q is new version
  LL query(int L, int R, int p) {
    return query(L, R, 0, n_, p);
  }
};
// https://www.luogu.com.cn/problem/P3919

// Bitree inside a Persistent Segment Tree to get dynamic k-th smallest number(online)
class BitPstSegTree {
  static inline constexpr int N = 1e9 + 2;
  struct Node {
    int l, r, val;
  };
  std::vector<Node> tree_;
  void pushUp(int p) {
    tree_[p].val = tree_[tree_[p].l].val + tree_[tree_[p].r].val;
  }
  void add(int x, int val, int l, int r, int p) {
    if (r - l == 1) {
      tree_[p].val += val;
    } else {
      int m = (l + r) / 2;
      if (x < m) {
        if (tree_[p].l == 0) tree_[p].l = newNode();
        add(x, val, l, m, tree_[p].l);
      } else {
        if (tree_[p].r == 0) tree_[p].r = newNode();
        add(x, val, m, r, tree_[p].r);
      }
      pushUp(p);
    }
  }
  int query(int k, int l, int r, std::vector<int>& p, std::vector<int>& q) {
    if (r - l == 1) return l;
    int m = (l + r) / 2, now = 0;
    for (auto x : q) now += tree_[tree_[x].l].val;
    for (auto x : p) now -= tree_[tree_[x].l].val;
    if (now >= k) {
      for (auto& x : q) x = tree_[x].l;
      for (auto& x : p) x = tree_[x].l;
      return query(k, l, m, p, q);
    }
    for (auto& x : q) x = tree_[x].r;
    for (auto& x : p) x = tree_[x].r;
    return query(k - now, m, r, p, q);
  }
  void add(int x, int val, int p) { add(x, val, 0, N, p); }
  int n_;
  std::vector<int> a_{0};
  int newNode() {
    int sz_ = tree_.size();
    tree_.emplace_back(Node());
    return sz_;
  }
  inline int lowbit(int x) { return x & -x; }
 public:
  BitPstSegTree(const std::vector<int>& x) : n_(x.size()) {
    // 0 is son of all leave nodes
    for (int i = 0; i <= n_; ++i) newNode();
    a_.insert(a_.end(), x.begin(), x.end());
    for (int i = 1; i <= n_; ++i) {
      for (int j = i; j <= n_; j += lowbit(j)) add(a_[i], 1, j);
    }
  }
  void modify(int x, int y) {
    for (int i = x; i <= n_; i += lowbit(i)) {
      add(a_[x], -1, i);
      add(y, 1, i);
    }
    a_[x] = y;
  }
  int query(int k, int l, int r) {
    std::vector<int> p, q;
    for (int i = --l; i; i -= lowbit(i)) p.emplace_back(i);
    for (int i = r; i; i -= lowbit(i)) q.emplace_back(i);
    return query(k, 0, N, p, q);
  }
};
// https://www.luogu.com.cn/problem/P2617
