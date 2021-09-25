#pragma once
#include <bits/stdc++.h>
using LL = long long;

// a will becomes next lexicographical order of a, satisfies $-1 < a_0 < a_1 < \cdots, a_{n - 1} < mx$
bool nextBinom(std::vector<int>& a, int mx) {
  int n = a.size(), i = 1;
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

// Returns the original value corresponding to the array value after discretization
template <typename T>
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
  for (int i = 1, n_ = a.size(); i < n_; ++i) {
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

// Disjoint Set Union
class DSU {
  std::vector<int> p_;
 public:
  DSU(int n) : p_(n) { iota(p_.begin(), p_.end(), 0); }
  int find(int x) {
    return x == p_[x] ? x : p_[x] = find(p_[x]);
  }
  bool merge(int x, int y) {
    int px = find(x), py = find(y);
    if (px == py) return false;
    // do something, small to big;
    return true;
  }
};


// Bit Tree Mininal version
template<typename T>
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

template<typename T>
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

template<typename T>
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

// There are two versions: sum and min/max
// min/max version slightly hard: you should record min/max, and second min/max value: https://codeforces.com/gym/102992/problem/J
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
    n_ = a.size();
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

// length of longest increasing subsquence
int LIS(std::vector<int>& a) {
  std::vector<int> b;
  for (auto x : a) {
    if (auto it = std::lower_bound(b.begin(), b.end(), x); it == b.end()) {
      b.emplace_back(x);
    } else
      *it = x;
  }
  return b.size();
}
// length of longest non-decreasing subsquence
int LNDS(std::vector<int>& a) {
  std::vector<int> b;
  for (auto x : a) {
    if (auto it = std::upper_bound(b.begin(), b.end(), x); it == b.end()) {
      b.emplace_back(x);
    } else
      *it = x;
  }
  return b.size();
}
// longest increasing subsquence
auto LISP(std::vector<int>& a) {
  std::vector<int> b, pb, pa(a.size());
  std::iota(pa.begin(), pa.end(), 0);
  for (int i = 0, na = a.size(); i < na; ++i) {
    if (auto it = std::upper_bound(b.begin(), b.end(), a[i]);
      it == b.end()) {
      if (!pb.empty()) pa[i] = pb.back();
      b.emplace_back(a[i]);
      pb.emplace_back(i);
    } else {
      *it = a[i];
      int t = it - b.begin();
      pb[t] = i;
      if (t > 0) pa[i] = pb[t - 1];
    }
  }
  std::stack<int> c;
  int now = pb.back();
  c.push(a[now]);
  while (now != pa[now]) {
    now = pa[now];
    c.push(a[now]);
  }
  return c;
}

// monicDeque: index of every max element of SubInterval of length m
std::vector<int> monicDequeMax(std::vector<int>& a, int m) {
  std::vector<int> r;
  std::deque<int> Q;
  for (int i = 0, na = a.size(); i < na; ++i) {
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
  int n = a.size();
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

// CDQ divided and conquer for partial order of dimension 3
struct cdqNode {
  int x, y, z, id, w;
  bool operator<(const cdqNode& A) const {
    if (x == A.x) return y == A.y ? z < A.z : y < A.y;
    return x < A.x;
  }
};
// ans[i] is the number of element less or equal to a[i]
std::vector<int> cdq(std::vector<cdqNode>& a, int k) {
  // sort by y
  std::vector<int> ans(a.size());
  std::sort(a.begin(), a.end());
  int last = 0;
  for (int i = 1, na = a.size(); i < na; ++i) {
    if (a[i].x != a[i - 1].x || a[i].y != a[i - 1].y ||
      a[i].z != a[i - 1].z) {
      int t = i - last - 1;
      for (int j = last; j < i; ++j) {
        ans[a[j].id] = t;
        a[j].w = 0;
      }
      a[i - 1].w = i - last;
      last = i;
    }
  }
  int t = a.size() - last - 1;
  for (int i = last, na = a.size(); i < na; ++i) {
    ans[a[i].id] = t;
    a[i].w = 0;
  }
  a.back().w = a.size() - last;
  Bitree<LL> A(k);
  auto cmpy = [](const cdqNode& lhs, const cdqNode& rhs) {
    return lhs.y < rhs.y;
  };
  std::function<void(int, int)> divide = [&](int l, int r) {
    if (r - l <= 1) return;
    int m = (l + r) / 2;
    divide(l, m);
    divide(m, r);
    std::sort(a.begin() + l, a.begin() + m, cmpy);
    std::sort(a.begin() + m, a.begin() + r, cmpy);
    int t = l;
    for (int i = m; i < r; ++i) {
      while (t < m && a[t].y <= a[i].y) {
        A.add(a[t].z, a[t].w);
        ++t;
      }
      ans[a[i].id] += A.sum(a[i].z);
    }
    for (int i = l; i < t; ++i) A.add(a[i].z, -a[i].w);
  };
  divide(0, a.size());
  return ans;
}
// https://www.luogu.com.cn/problem/P3810

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
