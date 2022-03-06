#pragma once
#include <bits/stdc++.h>

namespace Geomerty {
using Point = std::pair<double, double>;
double cross(const Point& op, const Point& sp, const Point& ep) {
  return (sp.first - op.first) * (ep.second - op.second)
  - (sp.second - op.second) * (ep.first - op.first);
}
bool crossLeft(const Point& op, const Point& sp, const Point& ep) {
  return (sp.first - op.first) * (ep.second - op.second)
  <= (sp.second - op.second) * (ep.first - op.first);
}
double dist2(const Point& p, const Point& q) {
  double x = q.first - p.first, y = q.second - p.second;
  return x * x + y * y;
};
double dist(const Point& p, const Point& q) {
  double x = q.first - p.first, y = q.second - p.second;
  return std::sqrt(x * x + y * y);
};

std::vector<Point> convexHull(std::vector<Point> p) {
  // note: parameter passing for argument of type 'std::pair<double, double>' 
  // when C++17 is enabled changed to match C++14 in GCC 10.1
  std::sort(p.begin(), p.end()); // compare with double stl_heap.h will change
  p.erase(std::unique(p.begin(), p.end()), p.end());
  int n = (int)p.size();
  std::vector<Point> q(n + 1);
  int top = 0;
  for (int i = 0; i < n; ++i) {
    while (top > 1 && crossLeft(q[top - 1], p[i], q[top - 2])) --top;
    q[top++] = p[i];
  }
  int len = top;
  for (int i = n - 2; i >= 0; --i) {
    while (top > len && crossLeft(q[top - 1], p[i], q[top - 2])) --top;
    q[top++] = p[i];
  }
  top -= n > 1;
  q.resize(top);
  return q;
}

double diameter(std::vector<Point> p) {
  auto q = convexHull(p);
  if (q.size() < 2) return 0;
  if (q.size() == 2) return dist2(q[0], q[1]);
  int n = (int)q.size();
  q.emplace_back(q[0]);
  double ans = 0;
  for (int i = 0, j = 2; i < n; ++i) {
    while (cross(q[i], q[i + 1], q[j]) < cross(q[i], q[i + 1], q[j + 1])) j = (j + 1) % n;
    ans = std::max({ans, dist2(q[i], q[j]), dist2(q[i + 1], q[j])});
  }
  return std::sqrt(ans);
} // float version: https://www.luogu.com.cn/problem/P6247
// Int version: https://www.luogu.com.cn/problem/P1452

double minDist(std::vector<Point> a) {
  double d = DBL_MAX;
  int n = (int)a.size();
  if (n <= 1) return d;
  std::sort(a.begin(), a.end());
  std::function<void(int, int)> merge = [&](int l, int r) {
    if (r - l <= 1) return;
    if (r - l == 2) {
      d = std::min(d, dist(a[l], a[l + 1]));
      return;
    }
    int m = (l + r) / 2;
    merge(l, m);
    merge(m + 1, r);
    std::vector<Point> p;
    for (int i = m - 1; i >= l && a[m].first - a[i].first < d; --i) {
      p.emplace_back(a[i].second, a[i].first);
    }
    for (int i = m; i < r && a[i].first - a[m].first < d; ++i) {
      p.emplace_back(a[i].second, a[i].first);
    }
    std::sort(p.begin(), p.end());
    for (int i = 0, np = (int)p.size(); i < np; ++i) {
      for (int j = i + 1; j < np && p[j].first - p[i].first < d; ++j) {
        d = std::min(d, dist(p[i], p[j]));
      }
    }
  };
  merge(0, n);
  return d;
} // https://www.luogu.com.cn/problem/P1429
} // namespace Geomerty

// a is k * n matrix: has n k-dimension vectors, return r[i]: number of vector less than i
std::vector<int> partialOrder(std::vector<std::vector<int>>& a) {
  int k = (int)a.size(), n = a[0].size();
  using Node = std::vector<std::pair<int, int>>;
  std::vector<Node> f(k, Node(n));
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < n; ++j) f[i][j] = {a[i][j], j};
    std::sort(f[i].begin(), f[i].end());
  }
  int sn = std::sqrt(n);
  constexpr int N = 4e4 + 2;
  using Data = std::vector<std::bitset<N>>;
  std::vector<Data> bs(k, Data(n / sn + 1));;
  for (int i = 0; i < k; ++i) {
    std::bitset<N> now;
    for (int j = 0; j < n; ++j) {
      if (j % sn == 0) bs[i][j / sn] = now;
      now.set(f[i][j].second);
    }
    if (n % sn == 0) bs[i][n / sn] = now;
  }
  auto getbst = [&](int i, int val) -> std::bitset<N> {
    int j = std::lower_bound(f[i].begin(), f[i].end(), std::make_pair(val, INT_MIN)) - f[i].begin();
    std::bitset<N> r = bs[i][j / sn];
    for (int t = j / sn * sn; t < j; ++t) r.set(f[i][t].second);
    return r;
  };
  std::vector<int> r(n);
  for (int j = 0; j < n; ++j) {
    std::bitset<N> now; now.set();
    for (int i = 0; i < k; ++i) {
      now &= getbst(i, a[i][j]);
    }
    r[j] = now.count();
  }
  return r;
} // http://cogs.pro:8081/cogs/problem/problem.php?pid=vSJzQVejP

namespace rectangleUnion {
class SegTree {
  static inline const int INF = 1e9 + 2;
  struct Node {
    int l = 0, r = 0, val = 0, mn = 0;
  };
  std::vector<Node> tree;
  int newNode() {
    int sz = tree.size();
    tree.emplace_back(Node());
    return sz;
  }
  void pushUp(int p) {
    tree[p].mn = std::min(tree[tree[p].l].mn, tree[tree[p].r].mn);
  }
  void updateNode(int p, int x) {
    tree[p].val += x;
    tree[p].mn += x;
  }
  void pushDown(int p) {
    if (tree[p].l) updateNode(tree[p].l, tree[p].val);
    if (tree[p].r) updateNode(tree[p].r, tree[p].val);
    tree[p].val = 0;
  }
  void add(int L, int R, int val, int l, int r, int p) {
    if (L <= l && r <= R) {
      updateNode(p, val);
      return;
    }
    auto m = (l + r) / 2;
    if (m > L) {
      if (tree[p].l == 0) tree[p].l = newNode();
      add(L, R, val, l, m, tree[p].l);
    }
    if (m < R) {
      if (tree[p].r == 0) tree[p].r = newNode();
      add(L, R, val, m, r, tree[p].r);
    }
  }
  int query(int l, int r, int p) {
    if (tree[p].mn) return r - l;
    auto m = (l + r) / 2;
    int ans = 0;
    if (tree[p].l) ans += query(l, m, tree[p].l);
    if (tree[p].r) ans += query(m, r, tree[p].r);
    return ans;
  }
 public:
  SegTree() {
    newNode();
  }
  void add(int L, int R, int val) {
    add(L, R, val, -INF, INF, 0);
  }
  int query() {
    return query(-INF, INF, 0);
  }
};

struct Edge {
  int x, l, r, val;
  bool operator<(const Edge& A) const {
    return x < A.x;
  }
};
LL rectangleUnion(const std::vector<std::tuple<int, int, int, int>>& rectangle) {
  std::vector<Edge> a;
  a.reserve(2 * rectangle.size());
  // make sure x1 < x2, y1 < y2
  for (auto [x1, y1, x2, y2] : rectangle) {
    a.push_back({x1, y1, y2, 1});
    a.push_back({x2, y1, y2, -1});
  }
  std::sort(a.begin(), a.end());
  SegTree A;
  A.add(a[0].l, a[0].r, a[0].val);
  LL ans = 0;
  for (int i = 1, n = a.size(); i < n; ++i) {
    if (a[i].x != a[i - 1].x) {
      ans += 1LL * (a[i].x - a[i - 1].x) * A.query();
    }
    A.add(a[i].l, a[i].r, a[i].val);
  }
  return ans;
} // https://www.luogu.com.cn/problem/T110664
} // namespace rectangleUnion
