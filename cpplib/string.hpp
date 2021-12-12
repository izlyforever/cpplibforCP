#pragma once
#include <bits/stdc++.h>
#include "template.hpp"
using LL = long long;

// Trie build by lowercase characters(change charToInt for other charset)
class Trie {
  using Node = std::array<int, 26>;
  std::vector<Node> nxt_;
  // 0: no point, 1: has point unvisited, 2: has point visited
  std::vector<int> tag_;
  void addNode(int fa, int c) {
    nxt_[fa][c] = nxt_.size();
    nxt_.emplace_back(Node());
    tag_.emplace_back(0);
  }
  int charToInt(char x) { return x - 'a';}
 public:
  Trie() : nxt_(1), tag_(1) {}
  void insert(std::string s) {
    int p = 0;
    for (auto x : s) {
      int c = x - 'a';
      if (nxt_[p][c] == 0) addNode(p, c);
      p = nxt_[p][c];
    }
    tag_[p] = 1;
  }
  int find(std::string s) {
    int p = 0;
    for (auto x : s) {
      int c = charToInt(x);
      if (nxt_[p][c] == 0) return 0;
      p = nxt_[p][c];
    }
    if (tag_[p] != 1) return tag_[p];
    tag_[p] = 2;
    return 1;
  }
};
// https://www.luogu.com.cn/problem/P2580

// 01 Trie find maximal xor sum
class Trie01 {
  using Node = std::array<int, 2>;
  std::vector<Node> ch_;
  void addNode(int fa, int c) {
    ch_[fa][c] = ch_.size();
    ch_.emplace_back(Node());
  }
 public:
  Trie01() : ch_(1) {}
  void insert(int x) {
    for (int i = 30, p = 0; i >= 0; --i) {
      int c = (x >> i) & 1;
      if (ch_[p][c] == 0) addNode(p, c);
      p = ch_[p][c];
    }
  }
  int getMax(int x) {
    int r = 0;
    for (int i = 30, p = 0; i >= 0; --i) {
      int c = (x >> i) & 1;
      if (ch_[p][c ^ 1]) {
        p = ch_[p][c ^ 1];
        r |= (1 << i);
      } else {
        p = ch_[p][c];
      }
    }
    return r;
  }
  int getAns(const std::vector<int>& a) {
    int r = 0;
    for (auto x : a) {
      insert(x);
      r = std::max(r, getMax(x));
    }
    return r;
  }
};
// https://www.luogu.com.cn/problem/P4551

// 01-Trie (Fusion Tree) xor sum(support modifty, add 1 to all)
class FusionTree {
  using Node = std::array<int, 4>;
  std::vector<Node> ch_;
  #define lsonFT ch_[p][0]
  #define rsonFT ch_[p][1]
  // ch_[p][2] = size of subtree rooted by p
  // ch_[p][3] = xor sum of subtree rooted by p
  void addNode(int p, int c) {
    ch_[p][c] = ch_.size();
    ch_.emplace_back(Node());
  }
  void pushUp(int p) {
    ch_[p][3] = 0;
    if (lsonFT) ch_[p][3] ^= (ch_[lsonFT][3] << 1);
    if (rsonFT) ch_[p][3] ^= (ch_[rsonFT][3] << 1) | (ch_[rsonFT][2] & 1);
  }
  // note ch_[lson][2] = ch_[p][2] - ch_[rsonFT] is a lazy tag_
  void insert(int p, int x) {
    ++ch_[p][2];
    if (!x) return;
    if (!ch_[p][x & 1]) addNode(p, x & 1);
    insert(ch_[p][x & 1], x >> 1);
    pushUp(p);
  }
  void erase(int p, int x) {
    --ch_[p][2];
    if (!x) return;
    erase(ch_[p][x & 1], x >> 1);
    pushUp(p);
  }
  void addAll(int p) {
    if (!ch_[p][2]) return;
    if (rsonFT) addAll(rsonFT);
    if (!lsonFT) addNode(p, 0);
    ch_[lsonFT][2] = ch_[p][2] - (rsonFT ? ch_[rsonFT][2] : 0);
    std::swap(lsonFT, rsonFT);
    pushUp(p);
  }
 public:
  FusionTree() : ch_(1) {}
  void insert(int x) {
    insert(0, x);
  }
  void erase(int x) {
    erase(0, x);
  }
  void addAll() {
    addAll(0);
  }
  int getVal() {
    return ch_[0][3];
  }
};
// https://www.luogu.com.cn/problem/P6018

std::vector<int> prefixFunction(std::string s) {
  int n = (int)s.size();
  std::vector<int> p(n);
  for (int i = 1; i < n; ++i) {
    int j = p[i - 1];
    while (j > 0 && s[i] != s[j]) j = p[j - 1];
    if (s[i] == s[j]) ++j;
    p[i] = j;
  }
  return p;
}

// KMP based on prefixFunction. return all match postion in t
std::vector<int> kmp(std::string s, std::string t) {
  std::vector<int> ans;
  int n = (int)s.size(), m = (int)t.size();
  if (n > m) return ans;
  auto p = prefixFunction(s);
  for (int i = 0, j = 0; i < m; ++i) {
    while (j > 0 && s[j] != t[i]) j = p[j - 1];
    if (s[j] == t[i] && ++j == n) ans.emplace_back(i - n + 1);
  }
  return ans;
}

// ans[i] = count(s[0, ..., i] occur in s)
std::vector<int> countPrefix(std::string s) {
  auto p = prefixFunction(s);
  int n = (int)s.size();
  std::vector<int> ans(n + 1);
  for (auto x : p) ++ans[x];
  for (int i = n - 1; i > 0; --i) ans[p[i - 1]] += ans[i];
  for (int i = 0; i <= n; ++i) ++ans[i];
  return ans;
}
// / ans[i] = count(s[0, ..., i] occur in t)
std::vector<int> countPrefix(std::string s, std::string t) {
  auto p = prefixFunction(s);
  int n = (int)s.size(), m = (int)t.size();
  std::vector<int> ans(n + 1);
  for (int i = 0, j = 0; i < m; ++i) {
    while (j > 0 && t[i] != s[j]) j = p[j - 1];
    if (t[i] == s[j]) ++j;
    ++ans[j];
  }
  ++ans[0];
  for (int i = n; i > 0; --i) ans[p[i - 1]] += ans[i];
  return ans;
}
// https://codeforces.com/problemset/problem/432/D

// z[i] = longest common prefix of s and s[i,...,n - 1], z[0] = 0.
std::vector<int> zFunction(std::string s) {
  int n = (int)s.size();
  std::vector<int> z(n);
  for (int i = 1, l = 0, r = 0; i < n; ++i) {
    if (i <= r && z[i - l] < r - i + 1) {
      z[i] = z[i - l];
    } else {
      int j = std::max(0, r - i + 1);
      while (i + j < n && s[j] == s[i + j]) ++j;
      z[i] = j;
      if (i + j - 1 > r) {
        l = i;
        r = i + j - 1;
      }
    }
  }
  return z;
}

// KMP based on zFunction. return all match postion in t
std::vector<int> kmpZ(std::string s, std::string t) {
  std::vector<int> ans;
  int n = (int)s.size(), m = (int)t.size();
  if (n > m) return ans;
  auto z = zFunction(s);
  for (int i = 0, l = 0, r = -1; i < m; ++i) {
    if (i > r || z[i - l] >= r - i + 1) {
      int j = std::max(0, r - i + 1);
      while (j < n && i + j < m && s[j] == t[i + j]) ++j;
      if (j == n) ans.emplace_back(i);
      if (i + j - 1 > r) {
        l = i;
        r = i + j - 1;
      }
    }
  }
  return ans;
}

class Automaton {
  static inline constexpr int CHAR = 26;
  using Node = std::array<int, CHAR>;
  std::vector<Node> nxt_;
  std::vector<int> cnt_, fail_, last_;
  int charToInt(char x) { return x - 'a';}
  void addNode(int fa, int c) {
    nxt_[fa][c] = nxt_.size();
    nxt_.emplace_back(Node());
    cnt_.emplace_back(0);
    fail_.emplace_back(0);
    last_.emplace_back(0);
  }
 public:
  Automaton() : nxt_(1), cnt_(1), fail_(1), last_(1) {}
  void insert(std::string s) {
    int p = 0;
    for (auto x : s) {
      int c = charToInt(x);
      if (nxt_[p][c] == 0) addNode(p, c);
      p = nxt_[p][c];
    }
    ++cnt_[p];
  }
  void build() {
    std::queue<int> Q;
    for (int c = 0; c < CHAR; ++c) {
      if (nxt_[0][c]) Q.push(nxt_[0][c]);
    }
    while (!Q.empty()) {
      int p = Q.front(); Q.pop();
      for (int c = 0; c < CHAR; ++c) {
        if (int& q = nxt_[p][c]; q != 0) {
          fail_[q] = nxt_[fail_[p]][c];
          Q.push(q);
          // match count optim
          last_[q] = cnt_[fail_[q]] ? fail_[q] : last_[fail_[q]];
        } else {
          q = nxt_[fail_[p]][c];
        }
      }
    }
  }
  // case-by-case analysis
  int query(std::string s) {
    int p = 0, r = 0;
    auto add = [&](int & x) {
      r += x; x = 0;
    };
    for (auto x : s) {
      int c = charToInt(x);
      p = nxt_[p][c];
      if (cnt_[p]) add(cnt_[p]);
      int q = p;
      while (last_[q]) {
        q = last_[q];
        if (cnt_[q]) add(cnt_[q]);
      }
    }
    return r;
  }
};
// https://www.luogu.com.cn/problem/P3808

// suffix array: SA-IS in O(N)
// assume a.back() = 0 and other elements are small postiove integral
std::vector<int> SAIS(std::vector<int> a) {
  enum TYPE {L, S};
  int n = (int)a.size() - 1, mx = *std::max_element(a.begin(), a.end()) + 1;
  std::vector<int> SA(n + 1, -1);
  std::vector<int> bucket(mx), lbucket(mx), sbucket(mx);
  for (auto x : a) ++bucket[x];
  for (int i = 1; i < mx; ++i) {
    bucket[i] += bucket[i - 1];
    lbucket[i] = bucket[i - 1];
    sbucket[i] = bucket[i] - 1;
  }
  // Determine the type of L, S and the position of type *
  std::vector<TYPE> type(n + 1);
  type[n] = S;
  for (int i = n - 1; i >= 0; --i) {
    type[i] = (a[i] < a[i + 1] ? S : (a[i] > a[i + 1] ? L : type[i + 1]));
  }
  // induce sort(from type * induces type L, from type L induce type S)
  auto inducedSort = [&]() {
    for (int i = 0; i <= n; ++i) {
      if (SA[i] > 0 && type[SA[i] - 1] == L) {
        SA[lbucket[a[SA[i] - 1]]++] = SA[i] - 1;
      }
    }
    for (int i = 1; i < mx; ++i) {
      sbucket[i] = bucket[i] - 1;
    }
    for (int i = n; i >= 0; --i) {
      if (SA[i] > 0 && type[SA[i] - 1] == S) {
        SA[sbucket[a[SA[i] - 1]]--] = SA[i] - 1;
      }
    }
  };
  // sort LMS-prefix according to induce sort
  std::vector<int> pos;
  for (int i = 1; i <= n; ++i) {
    if (type[i] == S && type[i - 1] == L) {
      pos.emplace_back(i);
    }
  }
  for (auto x : pos) SA[sbucket[a[x]]--] = x;
  inducedSort();
  // Give a name of LMS-substring(thus is S1) according to the order of LMS-prefix
  auto isLMSchar = [&](int i) {
    return i > 0 && type[i] == S && type[i - 1] == L;
  };
  auto equalSubstring = [&](int x, int y) {
    do {
      if (a[x] != a[y]) return false;
      ++x; ++y;
    } while (!isLMSchar(x) && !isLMSchar(y));
    return a[x] == a[y];
  };
  // Note: two LMS-substrings are the same only if they are adjacent since we have sort LMS-prefix
  std::vector<int> name(n + 1, -1);
  int lx = -1, cnt = 0;
  bool flag = true; // means no same LMS-substring
  for (auto x : SA) if (isLMSchar(x)) {
    if (lx >= 0 && !equalSubstring(lx, x)) ++cnt;
    if (lx >= 0 && cnt == name[lx]) flag = false;
    name[x] = cnt;
    lx = x;
  }
  std::vector<int> S1;
  for (auto x : name) if (x != -1) S1.emplace_back(x);
  auto getSA1 = [&]() {
    int n1 = S1.size();
    std::vector<int> SA1(n1);
    for (int i = 0; i < n1; ++i) SA1[S1[i]] = i;
    return SA1;
  };
  auto SA1 = flag ? getSA1() : SAIS(S1);
  // induce sort SA again according to the order of S1
  lbucket[0] = sbucket[0] = 0;
  for (int i = 1; i < mx; ++i) {
    lbucket[i] = bucket[i - 1];
    sbucket[i] = bucket[i] - 1;
  }
  std::fill(SA.begin(), SA.end(), -1);
  // scan SA1 in reverse order sicne S type is reorder in SA
  for (int i = SA1.size() - 1; i >= 0; --i) {
    SA[sbucket[a[pos[SA1[i]]]]--] = pos[SA1[i]];
  }
  inducedSort();
  return SA;
}

// the lexicographical smallest suffix is s[ans[0],...]
std::vector<int> SAIS(const std::string& s) {
  // If charset of s is lowercase letter, using th following f
  // auto f = [](char x) -> int { return int(x - 'a') + 1;};
  auto f = [](char x) -> int { return int(x) + 1;};
  std::vector<int> a;
  for (auto c : s) a.emplace_back(f(c));
  a.emplace_back(0);
  auto sa = SAIS(a);
  return std::vector<int>(sa.begin() + 1, sa.end());
}
// https://www.luogu.com.cn/problem/P3809

std::vector<int> getHeight(const std::string& s) {
  int n = (int)s.size();
  auto sa = SAIS(s);
  std::vector<int> rk(n);
  for (int i = 0; i < n; ++i) rk[sa[i]] = i;
  std::vector<int> ht(n);
  for (int i = 0, k = 0; i < n; ++i) {
    if (k) --k;
    if (rk[i] == 0) k = 0;
    else {
      // won't out of range
      while (s[i + k] == s[sa[rk[i] - 1] + k]) ++k;
      ht[rk[i]] = k;
    }
  }
  return ht;
}

LL diffSubstringCount(const std::string& s) {
  int n = (int)s.size();
  auto ht = getHeight(s);
  return 1LL * (n + 1) * n / 2 - std::accumulate(ht.begin(), ht.end(), 0LL);
}
// https://www.luogu.com.cn/problem/P2408

template<typename T, typename enable = IntegerT<T>>
int minPresent(std::vector<T>& a) {
  int k = 0, i = 0, j = 1, n = (int)a.size();
  while (k < n && i < n && j < n) {
    if (a[(i + k) % n] == a[(j + k) % n]) {
      ++k;
    } else {
      a[(i + k) % n] > a[(j + k) % n] ? i += k + 1 : j += k + 1;
      if (i == j) ++i;
      k = 0;
    }
  }
  return std::min(i, j);
}

// Lyndon decomposition using Duval algorithm
std::vector<std::string> duval(const std::string& s) {
  std::vector<std::string> r;
  int n = (int)s.size(), i = 0;
  while (i < n) {
    int j = i + 1, k = i;
    while (j < n && s[k] <= s[j]) {
      if (s[k] < s[j]) k = i;
      else ++k;
      ++j;
    }
    while (i <= k) {
      r.emplace_back(s.substr(i, j - k));
      i += j - k;
    }
  }
  return r;
}
// https://www.luogu.com.cn/problem/P6114

// s(i - ans[i], i + ans[i])  is longest Palindrome
std::vector<int> Manacher(std::string s) {
  int n = (int)s.size(); // assume n is odd
  std::vector<int> d(n);
  for (int i = 0, l = 0, r = -1; i < n; ++i) {
    int k = i > r ? 1 : std::min(d[l + r - i], r - i);
    while (k <= i && i + k < n && s[i - k] == s[i + k]) ++k;
    d[i] = k--;
    if (i + k > r) {
      l = i - k;
      r = i + k;
    }
  }
  return d;
}
// add a special char to deal with the case of even length
