#pragma once
#include <bits/stdc++.h>
#include "fenwick.hpp"
using LL = long long;

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
  for (int i = 1, na = (int)a.size(); i < na; ++i) {
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
  int t = (int)a.size() - last - 1;
  for (int i = last, na = (int)a.size(); i < na; ++i) {
    ans[a[i].id] = t;
    a[i].w = 0;
  }
  a.back().w = (int)a.size() - last;
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
