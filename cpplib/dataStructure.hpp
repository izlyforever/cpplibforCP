#pragma once
#include <bits/stdc++.h>
using LL = long long;

// a will becomes next lexicographical order of a, satisfies $-1 < a_0 < a_1 < \cdots, a_{n - 1} < mx$
bool nextBinom(std::vector<int> &a, int mx) {
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
    for (auto x : a) std::cout << x << " ";
    std::cout << "\n";
  } while (nextBinom(a, mx));
}

// Error Correction Code: O(n m + k^k n)
class ECC {
  std::vector<std::vector<int>> a;  // origin data: n rows, m cols.
  int k;                // Maximum number of differences allowed
  std::vector<std::vector<int>> bad;  // difference with current answer
  int n, m, mxId;
  void updateMxId(int i) {
    if (bad[i].size() > bad[mxId].size()) mxId = i;
  }
  bool dfs(int c) {  // remain time that current answer can change
    auto bd = bad[mxId];
    if ((int)bd.size() <= k) return true;
    if ((int)bd.size() - k > c) return false;
    // Note that bd is O(k) instead of O(m)
    std::vector<int> f(bd.size() - k);
    iota(f.begin(), f.end(), 0);
    int tMxId = mxId;
    do {
      mxId = tMxId;
      std::queue<int> tmp;
      for (auto x : f) {
        tmp.push(r[bd[x]]);
        for (int i = 0; i < n; ++i) {
          if (a[i][bd[x]] == r[bd[x]]) {
            bad[i].emplace_back(bd[x]);
          }
          if (a[i][bd[x]] == a[mxId][bd[x]]) {
            bad[i].erase(
              std::find(bad[i].begin(), bad[i].end(), bd[x]));
          }
        }
        r[bd[x]] = a[mxId][bd[x]];
      }
      for (int i = 0; i < n; ++i) updateMxId(i);
      if (dfs(c - f.size())) return true;
      for (auto x : f) {
        for (int i = 0; i < n; ++i) {
          if (a[i][bd[x]] == r[bd[x]]) {
            bad[i].emplace_back(bd[x]);
          }
          if (a[i][bd[x]] == tmp.front()) {
            bad[i].erase(
              std::find(bad[i].begin(), bad[i].end(), bd[x]));
          }
        }
        r[bd[x]] = tmp.front();
        tmp.pop();
      }
    } while (nextBinom(f, bd.size()));
    return false;
  }
public:
  std::vector<int> r;   // m cols vector, current answer
  ECC(std::vector<std::vector<int>> _a) : a(_a), r(a[0]) {
    n = a.size(); m = r.size();
    bad.resize(n);
    mxId = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j)
        if (a[i][j] != r[j]) {
          bad[i].emplace_back(j);
        }
      updateMxId(i);
    }
  }
  void setK(int _k) { k = _k; }
  bool solve() { return dfs(k); }
};

// Returns the original value corresponding to the array value after discretization
template <typename T>
std::vector<T> discrete(std::vector<T> &a) {
  auto b = a;
  std::sort(b.begin(), b.end());
  b.erase(std::unique(b.begin(), b.end()), b.end());
  for (auto &x : a) {
    x = std::lower_bound(b.begin(), b.end(), x) - b.begin();
  }
  return b;
}

// Disjoint Set Union
class DSU {
  std::vector<int> p;
public:
  DSU(int n) : p(n) { iota(p.begin(), p.end(), 0); }
  int find(int x) {
    return x == p[x] ? x : p[x] = find(p[x]);
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
  std::vector<T> s;
  BitreeMin() {}
  BitreeMin(int n) : s(n + 1, std::numeric_limits<T>::max()) {}
  int lowbit(int n) { return n & (-n); }
  void modify(int id, T p) {
    int ns = s.size();
    while (id < ns) {
      s[id] = std::min(s[id], p);
      id += lowbit(id);
    }
  }
  // cal minial value in [1, id]
  T min(int id) {
    T r = std::numeric_limits<T>::max()
    while (id) {
      r = std::min(r, s[id]);
      id -= lowbit(id);
    }
    return r;
  }
};

template<typename T>
struct Bitree {
  std::vector<T> s;
  Bitree() {}
  Bitree(int n) : s(n + 1) {}
  int lowbit(int n) { return n & (-n); }
  void add(int id, T p) {
    int ns = s.size();
    while (id < ns) {
      s[id] += p;
      id += lowbit(id);
    }
  }
  T sum(int id) {
    T r = 0;
    while (id) {
      r += s[id];
      id -= lowbit(id);
    }
    return r;
  }
  // find minimal index s.t. sum(id) >= x, sum must be increased
  int search(T val) {
    T sum = 0;
    int id = 0;
    for (int i = std::__lg(s.size()); ~i; --i) {
      if (id + (1 << i) < (int)s.size() && sum + s[id + (1 << i)] < val) {
        id += (1 << i);
        sum += s[id];
      }
    }
    return ++id;
  }
};

template<typename T>
class BitreePlus {
  int n;
  // c[i] = a[i] - a[i - 1], b_i = (i - 1) * c_i
  Bitree<T> B, C;
  void add(int id, T p) {
    C.add(id, p);
    B.add(id, (id - 1) * p);
  }
public:
  BitreePlus() {}
  BitreePlus(int _n) : n(_n), B(n), C(n) {}
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
    for (int i = std::__lg(n); ~i; --i)
      if (int idi = id + (1 << i); idi <= n) {
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
  int n;
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
    tag.resize(4 * n);
    sm.resize(4 * n);
  }
public:
  SegmentTree(int _n) : n(_n) { resize(); }
  SegmentTree(const std::vector<LL> &a) {
    n = a.size();
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
    build(0, n, 1);
  }
  void add(int L, int R, LL v) { rangeAdd(--L, R, v, 0, n, 1); }
  LL query(int L, int R) { return query(--L, R, 0, n, 1); }
};
// https://www.luogu.com.cn/problem/P3372

// Persistent Segment Tree, reference: https://zhuanlan.zhihu.com/p/250565583
class PstSegTree {
  struct Node {
    int l, r;
    LL val;
  };
  void pushUp(int p) {
    tree[p].val = tree[tree[p].l].val + tree[tree[p].r].val;
  }
public:
  int n;
  std::vector<int> root;  // version number
  std::vector<Node> tree;
  int newNode() {
    int sz = tree.size();
    tree.emplace_back(Node());
    return sz;
  }
  PstSegTree(const std::vector<int> &a) : n(a.size()) {
    root.emplace_back(newNode());
    std::function<void(int, int, int)> build = [&](int l, int r, int p) {
      if (r - l == 1) {
        tree[p].val = a[l];
      } else {
        int m = (l + r) / 2;
        build(l, m, tree[p].l = newNode());
        build(m, r, tree[p].r = newNode());
        pushUp(p);
      }
    };
    build(0, n, root.back());
  }
  // single point update, p is current version, q is new version
  void update(int pos, int val, int l, int r, int p, int q) {
    tree[q] = tree[p];
    if (r - l == 1) {
      tree[q].val = val;
    } else {
      int m = (l + r) / 2;
      if (pos < m)
        update(pos, val, l, m, tree[p].l, tree[q].l = newNode());
      else
        update(pos, val, m, r, tree[p].r, tree[q].r = newNode());
      pushUp(q);
    }
  }
  // segment query, p is current version, q is new version
  LL query(int L, int R, int l, int r, int p) {
    if (L <= l && R >= r) return tree[p].val;
    int m = (l + r) / 2;
    LL ans = 0;
    if (L < m) ans += query(L, R, l, m, tree[p].l);
    if (R > m) ans += query(L, R, m, r, tree[p].r);
    return ans;
  }
};
// https://www.luogu.com.cn/problem/P3919

// Bitree inside a Persistent Segment Tree to get dynamic k-th smallest number(online)
class BitPstSegTree {
  static inline constexpr int N = 1e9 + 2;
  struct Node {
    int l, r, val;
  };
  std::vector<Node> tree;
  void pushUp(int p) {
    tree[p].val = tree[tree[p].l].val + tree[tree[p].r].val;
  }
  void add(int x, int val, int l, int r, int p) {
    if (r - l == 1) {
      tree[p].val += val;
    } else {
      int m = (l + r) / 2;
      if (x < m) {
        if (tree[p].l == 0) tree[p].l = newNode();
        add(x, val, l, m, tree[p].l);
      } else {
        if (tree[p].r == 0) tree[p].r = newNode();
        add(x, val, m, r, tree[p].r);
      }
      pushUp(p);
    }
  }
  int query(int k, int l, int r, std::vector<int> &p, std::vector<int> &q) {
    if (r - l == 1) return l;
    int m = (l + r) / 2, now = 0;
    for (auto x : q) now += tree[tree[x].l].val;
    for (auto x : p) now -= tree[tree[x].l].val;
    if (now >= k) {
      for (auto &x : q) x = tree[x].l;
      for (auto &x : p) x = tree[x].l;
      return query(k, l, m, p, q);
    }
    for (auto &x : q) x = tree[x].r;
    for (auto &x : p) x = tree[x].r;
    return query(k - now, m, r, p, q);
  }
  void add(int x, int val, int p) { add(x, val, 0, N, p); }
  int n;
  std::vector<int> a{0};
  int newNode() {
    int sz = tree.size();
    tree.emplace_back(Node());
    return sz;
  }
  inline int lowbit(int x) { return x & -x; }
public:
  BitPstSegTree(const std::vector<int> &x) : n(x.size()) {
    // 0 is son of all leave nodes
    for (int i = 0; i <= n; ++i) newNode();
    a.insert(a.end(), x.begin(), x.end());
    for (int i = 1; i <= n; ++i) {
      for (int j = i; j <= n; j += lowbit(j)) add(a[i], 1, j);
    }
  }
  void modify(int x, int y) {
    for (int i = x; i <= n; i += lowbit(i)) {
      add(a[x], -1, i);
      add(y, 1, i);
    }
    a[x] = y;
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
int LIS(std::vector<int> &a) {
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
int LNDS(std::vector<int> &a) {
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
auto LISP(std::vector<int> &a) {
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
std::vector<int> monicDequeMax(std::vector<int> &a, int m) {
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
std::vector<int> monicStack(const std::vector<int> &a) {
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
int cartesian_build(std::vector<cNode> &tree, int n) {
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
  bool operator<(const cdqNode &A) const {
    if (x == A.x) return y == A.y ? z < A.z : y < A.y;
    return x < A.x;
  }
};
// ans[i] is the number of element less or equal to a[i]
std::vector<int> cdq(std::vector<cdqNode> &a, int k) {
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
  auto cmpy = [](const cdqNode &lhs, const cdqNode &rhs) {
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
  int l, r; // fa, sz \in [l, r]
  int f, d; // x \in [l, r] has real value f x - d, where f = 1, -1
  std::vector<int> fa;
  int find(int x) {
    return x == fa[x] ? x : fa[x] = find(fa[x]);
  }
  void merge(int x, int y) { // merge x to y
    fa[find(x)] = find(y);
  }
public:
  BlockAbs(int mx) : l(0), r(mx), f(1), d(0), fa(mx + 1) {
    std::iota(fa.begin(), fa.end(), 0);
  }
  void add(int x) { // |fi - d - x| = | i - (x + d) f|
    x = (x + d) * f;
    if ((l + r) < 2 * x) {
      f = -1;
      d = -x;
      if (x < r) {
        for (int i = r; i > x; --i) merge(i, 2 * x - i);
        r = x;
      }
    } else {
      f = 1;
      d = x;
      if (x > l) {
        for (int i = l; i < x; ++i) merge(i, 2 * x - i);
        l = x;
      }
    }
  }
  int query(int x) {
    return find(x) * f - d;
  }
};
// https://codeforces.com/gym/103104/problem/K

// Second Block minus version(online, space optim can be done if offline)
class BlockMinus {
  std::vector<int> fa, sz, a;
  int l, delta, mx; // real value x - delta
  int find(int x) {
    return x == fa[x] ? x : fa[x] = find(fa[x]);
  }
  void merge(int x, int y) { // merge x to y
    x = find(x); y = find(y);
    if (x == y) return;
    fa[x] = y;
    sz[y] += sz[x];
    sz[x] = 0;
  }
  void modifyPart(int ql, int qr, int x) {
    for (int i = ql; i < qr; ++i) {
      a[i] = find(a[i]);
      if (a[i] - delta > x) {
        --sz[a[i]];
        a[i] = find(a[i] - x);
        ++sz[a[i]];
      }
    }
  }
  void modifyAll(int x) {
    if (x < mx - delta - x) {
      for (int i = delta + 1; i <= x + delta; ++i) merge(i, i + x);
      delta += x;
    } else {
      for (int i = mx; i > x + delta; --i) merge(i, i - x);
      mx = x + delta;
    }
  }
  int queryPart(int ql, int qr, int x) {
    int ans = 0;
    for (int i = ql; i < qr; ++i) {
      if (find(a[i]) - delta == x) ++ans;
    }
    return ans;
  }
  int queryAll(int x) {
    x += delta;
    if (x > mx || find(x) != x) return 0;
    return sz[x];
  }
public:
  void init(const std::vector<int> &_a, int _l, int _r) {
    l = _l, delta = 0;
    a = {_a.begin() + _l, _a.begin() + _r};
    mx = *std::max_element(a.begin(), a.end());
    fa.resize(mx + 1); std::iota(fa.begin(), fa.end(), 0);
    sz.resize(mx + 1); for (auto x : a) ++sz[x];
  }
  void modify(int ql, int qr, int x) {
    if (x >= mx - delta) return;
    if (qr - ql == (int)a.size()) modifyAll(x);
    else modifyPart(ql - l, qr - l, x);
  }
  int query(int ql, int qr, int x) {
    if (qr - ql == (int)a.size()) return queryAll(x);
    return queryPart(ql - l, qr - l, x);
  }
};
// https://codeforces.com/contest/896/problem/E
