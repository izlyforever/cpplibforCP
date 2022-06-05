#pragma once
#include <bits/stdc++.h>

// Disjoint Set Union
class DSU {
  std::vector<int> p_;
  // std::vector<int> sz_;
 public:
  DSU(int n) : p_(n) { iota(p_.begin(), p_.end(), 0); }
  int find(int x) {
    return x == p_[x] ? x : p_[x] = find(p_[x]);
    // other impl not recomand
    // while (x != p_[x]) x = p_[x] = p_[p_[x]];
    // return x;
  }
  bool merge(int x, int y) {
    int px = find(x), py = find(y);
    if (px == py) return false;
    // do something, small to big;
    return true;
  }
};
