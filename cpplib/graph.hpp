#pragma once
#include <bits/stdc++.h>
using LL = long long;

using edge = std::vector<std::pair<int, int>>;
using Edge = std::tuple<int, int, int>;

// dfs order of rooted tree
class DfsTour {
  int n_, cnt_;
  std::vector<int> l_, r_;
  std::vector<std::vector<int>> e_;
 public:
  DfsTour(int n) : n_(n), cnt_(0), l_(n), r_(n), e_(n) {}
  void addEdge(int u, int v) {
    if (u == v) return;
    e_[u].emplace_back(v);
    e_[v].emplace_back(u);
  }
  void dfs(int u, int fa) {
    l_[u] = ++cnt_;
    for (auto v : e_[u]) if (v != fa) {
      dfs(v, u);
    }
    r_[u] = cnt_;
  }
};

// Euler Tour: rooted by rt(length 2n - 1)
std::vector<int> EulerTour(std::vector<std::vector<int>>& e, int rt) {
  std::vector<int> r;
  std::function<void(int, int)> dfs = [&](int u, int fa) -> void {
    r.emplace_back(u);
    for (auto v : e[u]) if (v != fa) {
      dfs(v, u);
      r.emplace_back(u);
    }
  };
  dfs(rt, rt);
  return r;
}

class LCA {
  int n_;
  std::vector<int> fa_, dep_, sz_, son_, top_;
 public:
  LCA(std::vector<std::vector<int>>& e, int rt = 1) : n_(e.size()) {
    fa_.resize(n_);
    dep_.resize(n_);
    sz_.resize(n_);
    son_.resize(n_);
    fa_[rt] = rt;
    dep_[rt] = 0;
    std::function<int(int)> pdfs = [&](int u) -> int {
      sz_[u] = 1;
      for (auto v : e[u]) if (v != fa_[u]) {
        dep_[v] = dep_[u] + 1;
        fa_[v] = u;
        sz_[u] += pdfs(v);
        if (sz_[v] > sz_[son_[u]]) son_[u] = v;
      }
      return sz_[u];
    };
    top_.resize(n_);
    std::function<void(int, int)> dfs = [&](int u, int t) -> void {
      top_[u] = t;
      if (son_[u] == 0) return;
      dfs(son_[u], t);
      for (auto v : e[u]) if (v != fa_[u] && v != son_[u]) dfs(v, v);
    };
    pdfs(rt);
    dfs(rt, rt);
  }
  int lca(int u, int v) {
    while (top_[u] != top_[v]) {
      if (dep_[top_[u]] > dep_[top_[v]]) {
        u = fa_[top_[u]];
      } else {
        v = fa_[top_[v]];
      }
    }
    return dep_[u] < dep_[v] ? u : v;
  }
};

// Minimum Spanning Tree
LL Prim(const std::vector<edge>& e) {
  LL r = 0;
  int n = (int)e.size(), cnt = 0;
  std::priority_queue<std::pair<int, int>> Q;
  std::vector<int> vis(n);
  Q.push({0, 0});
  while (!Q.empty()) {
    auto [w, u] = Q.top();
    Q.pop();
    if (vis[u]) continue;
    ++cnt;
    r -= w;
    vis[u] = 1;
    for (auto [v, c] : e[u]) if (!vis[v]) {
      Q.push({-c, v});
    }
  }
  return cnt == n ? r : INT64_MAX;
}

// LiuZhu: Minimum tree diagram
LL LiuZhu(std::vector<Edge> e, int n, int rt) { // e has no self-loop
  LL ans = 0;
  while (1) {
    std::vector<int> in(n, INT_MAX), pre(n, -1);
    for (auto [u, v, w] : e) if (u != v && in[v] > w) {
      in[v] = w;
      pre[v] = u;
    }
    for (int i = 0; i < n; ++i) {
      if (i != rt && pre[i] == -1) return -1;
    }
    // judge if there is a loop
    int cnt = 0;
    std::vector<int> vis(n, -1), id(n, -1);
    for (int i = 0; i < n; ++i) if (i != rt) {
      ans += in[i];
      int v = i;
      // note there may be a path_ of form '6'
      while (vis[v] != i && id[v] == -1 && v != rt) {
        vis[v] = i;
        v = pre[v];
      }
      if (id[v] == -1 && v != rt) {
        int u = v;
        do {
          id[u] = cnt;
          u = pre[u];
        } while (u != v);
        ++cnt;
      }
    }
    if (cnt == 0) break;
    // update nodes and edges
    for (int i = 0; i < n; ++i) if (id[i] == -1) id[i] = cnt++;
    for (auto& [u, v, w] : e) {
      if (id[u] != id[v]) w -= in[v];
      u = id[u];
      v = id[v];
    }
    rt = id[rt];
    n = cnt;
  }
  return ans;
}


// Topological sorting: Kahn algorithm for DAG. return lexicographical smallest one
std::vector<int> TopSort(std::vector<std::set<int>>& e) {
  std::vector<int> r;
  std::priority_queue<int> Q;
  int n = (int)e.size();
  std::vector<int> in(n);
  for (auto& x : e) for (auto i : x) ++in[i];
  for (int i = 0; i < n; ++i) if (in[i] == 0) Q.push(-i);
  while (!Q.empty()) {
    int u = -Q.top();
    r.emplace_back(u);
    Q.pop();
    for (auto v : e[u]) {
      if (--in[v] == 0) {
        Q.push(-v);
      }
    }
  }
  return r;
}
// https://www.luogu.com.cn/problem/U107394

// Euler Path: return lexicographical smallest one or empty(use set if there is no multi-edge)
std::stack<int> EulerPathS(std::vector<std::multiset<int>> e) {
  int cnt = std::count_if(e.begin(), e.end(), [](auto x) {
    return x.size() & 1;
  });
  if (cnt > 2) return std::stack<int>();
  std::stack<int> ans;
  std::function<void(int)> Hierholzer = [&](int u) {
    while (!e[u].empty()) {
      int v = *e[u].begin();
      e[u].erase(e[u].begin());
      e[v].erase(e[v].find(u));
      Hierholzer(v);
    }
    ans.push(u);
  };
  for (int i = 0, ne = (int)e.size(); i < ne; ++i) {
    if (!e[i].empty() && ((e[i].size() & 1) || (cnt == 0))) {
      Hierholzer(i);
      break;
    }
  }
  return ans;
}
// Euler Path start form rt: return lexicographical smallest one(assume existence and no multi-edge)
std::stack<int> EulerPath(std::vector<std::set<int>> e, int rt) {
  std::stack<int> ans;
  std::function<void(int)> Hierholzer = [&](int u) {
    while (!e[u].empty()) {
      int v = *e[u].begin();
      e[u].erase(e[u].begin());
      e[v].erase(e[v].find(u));
      Hierholzer(v);
    }
    ans.push(u);
  };
  Hierholzer(rt);
  return ans;
}
// https://www.luogu.com.cn/problem/P2731

namespace Floyd {
constexpr int N = 1003;
LL dp[N][N], path[N][N];
void Floyd(int n) {
  memset(path, -1, sizeof(path));
  for(int k = 0; k != n; ++k)
    for(int i = 0; i != n; ++i)
      for(int j = 0; j != n; ++j) if (dp[i][j] > dp[i][k] + dp[k][j]) {
        path[i][j] = k;
      }
}
std::vector<int> getPath(int x, int y) {
  if (path[x][y] == -1) {
    if (x == y) return std::vector<int>{x};
    return std::vector<int>{x, y};
  }
  auto left = getPath(x, path[x][y]);
  auto now = getPath(path[x][y], y);
  left.insert(left.end(), now.begin(), now.end());
  return left;
}
} // namespace Floyd

std::vector<LL> Dijkstra(const std::vector<edge>& e, int s) {
  std::priority_queue<std::pair<LL, int>> Q;
  std::vector<LL> d(e.size(), INT64_MAX);
  d[s] = 0;
  Q.push({0, s});
  while (!Q.empty()) {
    auto [du, u] = Q.top();
    Q.pop();
    if (d[u] != -du) continue;
    for (auto [v, w] : e[u]) if (d[v] > d[u] + w) {
      d[v] = d[u] + w;
      Q.emplace(-d[v], v);
    }
  }
  return d;
}

bool BellmanFord(std::vector<Edge>& e, int n, int s = 0) {
  std::vector<int> dist(n + 1, INT_MAX);
  dist[s] = 0;
  for (int i = 0; i <= n; ++i) {
    bool judge = false;
    for (auto [u, v, w] : e) if (dist[u] != INT_MAX) {
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        judge = true;
      }
    }
    if (!judge) return true;
  }
  return false;
}

bool spfa(std::vector<edge>& e, int s = 0) {
  int n = (int)e.size();
  std::queue<int> Q;
  std::vector<int> dist(n, INT_MAX), cnt(n), inQ(n);
  Q.push(s);
  inQ[s] = 1;
  dist[s] = 0;
  ++cnt[s];
  while (!Q.empty()) {
    int u = Q.front();
    Q.pop();
    inQ[u] = 0;
    for (auto [v, w]: e[u]) {
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        if (!inQ[v]) {
          Q.push(v);
          inQ[v] = 1;
          if (++cnt[v] == n) return false;
        }
      }
    }
  }
  return true;
}

// Kosaraju: Strongly Connected Components
class Scc {
  int n_, nScc_;
  std::vector<int> vis_, color_, order_;
  std::vector<std::vector<int>> e_, e2_;
  void dfs(int u) {
    vis_[u] = true;
    for (auto v : e_[u]) if (!vis_[v]) dfs(v);
    order_.emplace_back(u);
  }
  void dfs2(int u) {
    color_[u] = nScc_;
    for (auto v : e2_[u]) if (!color_[v]) dfs2(v);
  }
  void Kosaraju() {
    for (int i = 0; i < n_; ++i) if (!vis_[i]) dfs(i);
    for (auto it = order_.rbegin(); it != order_.rend(); ++it) if (!color_[*it]) {
      ++nScc_;
      dfs2(*it);
    }
  }
 public:
  Scc(int n) : n_(n) {
    nScc_ = 0;
    e_.resize(n_);
    e2_.resize(n_);
    vis_.resize(n_);
    color_.resize(n_);
  }
  void addEdge(int u, int v) {
    nScc_ = 0;
    e_[u].emplace_back(v);
    e2_[v].emplace_back(u);
  }
  int solve() {
    if (nScc_ == 0) Kosaraju();
    return nScc_;
  }
};

// n_ / 2 pairs (2i, 2i + 1)
struct twoSAT {
  int n_, nScc_;
  std::vector<int> vis_, color_, order_;
  std::vector<std::vector<int>> e_, e2_;
  twoSAT(int n) : n_(n * 2) {
    nScc_ = 0;
    e_.resize(n_);
    e2_.resize(n_);
    vis_.resize(n_);
    color_.resize(n_);
  }
  void addEdge(int u, int v) {
    e_[u].emplace_back(v);
    e2_[v].emplace_back(u);
  }
  void dfs(int u) {
    vis_[u] = true;
    for (auto v : e_[u]) if (!vis_[v]) dfs(v);
    order_.emplace_back(u);
  }
  void dfs2(int u) {
    color_[u] = nScc_;
    for (auto v : e2_[u]) if (!color_[v]) dfs2(v);
  }
  void Kosaraju() {
    for (int i = 0; i < n_; ++i) if (!vis_[i]) dfs(i);
    for (auto it = order_.rbegin(); it != order_.rend(); ++it) if (!color_[*it]) {
      ++nScc_;
      dfs2(*it);
    }
  }
  std::vector<int> solve() {
    Kosaraju();
    // choose component with max color_ number
    std::vector<int> choose(nScc_ + 1);
    for (int i = 0; i < n_; i += 2) {
      int c1 = color_[i], c2 = color_[i + 1];
      if (c1 == c2) return std::vector<int>();
      if (choose[c1] || choose[c2]) continue;
      choose[std::max(c1, c2)] = 1;
    }
    std::vector<int> r(n_ / 2);
    for (int i = 0; i * 2 < n_; ++i) r[i] = (choose[color_[i * 2]] ? 1 : -1);
    return r;
  }
};

std::vector<int> cutVertex(std::vector<std::vector<int>>& e) {
  int n = (int)e.size(), cnt = 0;
  std::vector<int> dfs(n), low(n), flag(n), r;
  std::function<void(int, int)> Tarjan = [&](int u, int fa) -> void {
    low[u] = dfs[u] = ++cnt;
    int ch = 0;
    for (auto v : e[u]) {
      if (dfs[v] == 0) {
        ++ch;
        Tarjan(v, u);
        low[u] = std::min(low[u], low[v]);
        if (u != fa && low[v] >= dfs[u]) flag[u] = 1;
      } else if (v != fa) {
        low[u] = std::min(low[u], dfs[v]);
      }
    }
    if (u == fa && ch > 1) flag[u] = 1;
  };
  for (int i = 0; i < n; ++i) if (dfs[i] == 0) Tarjan(i, i);
  for (int i = 0; i < n; ++i) if (flag[i]) r.emplace_back(i);
  return r;
}
// https://www.luogu.com.cn/problem/P3388

class CutEdge {
  int n_, cnt_;
  std::vector<std::vector<int>> g_;
  std::vector<int> e, flag, dfs, low;
  void Tarjan(int u, int inEdgeNum) {
    low[u] = dfs[u] = ++cnt_;
    for (auto i : g_[u]) {
      int v = e[i];
      if (dfs[v] == 0) {
        Tarjan(v, i);
        low[u] = std::min(low[u], low[v]);
        if (low[v] > dfs[u]) flag[i] = flag[i ^ 1] = 1;
      } else if ((i ^ 1) != inEdgeNum) {
        low[u] = std::min(low[u], dfs[v]);
      }
    }
  }
 public:
  CutEdge(int n) : n_(n), cnt_(0), g_(n), dfs(n), low(n) {}
  void addEdge(int u, int v) {
    if (u == v) return;
    g_[u].emplace_back(e.size());
    e.emplace_back(v);
    flag.emplace_back(0);
    g_[v].emplace_back(e.size());
    e.emplace_back(u);
    flag.emplace_back(0);
  }
  int solve() {
    for (int i = 0; i < n_; ++i) if (dfs[i] == 0) Tarjan(i, -1);
    int r = 0;
    for (auto x : flag) r += x;
    return r / 2;
  }
};
// https://www.luogu.com.cn/problem/T103481

// S-T max-Flow: Dinic $O(n_^2 m)$
class Dinic {
  int n_;
  // e[i] = {endPoint, conpacity} and e[i ^ 1] is opposite edge of e[i]
  // g_[u] = {edges start from u}
  std::vector<std::pair<int, int>> e;
  std::vector<std::vector<int>> g_;
  std::vector<int> cur_, h_;
  // h_[i] = dist(s, i), so h_[t] != -1 means there is a path from s to t
  bool bfs(int s, int t) {
    h_.assign(n_, -1);
    std::queue<int> Q;
    h_[s] = 0;
    Q.push(s);
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      for (auto i : g_[u]) {
        auto [v, c] = e[i];
        if (c > 0 && h_[v] == -1) {
          h_[v] = h_[u] + 1;
          Q.push(v);
        }
      }
    }
    return h_[t] != -1;
  }
  LL dfs(int u, int t, LL f) {
    if (u == t || f == 0) return f;
    LL r = f;
    for (int& i = cur_[u], ng = g_[u].size(); i < ng; ++i) {
      int j = g_[u][i];
      auto [v, c] = e[j];
      if (c > 0 && h_[v] == h_[u] + 1) {
        int a = dfs(v, t, std::min(r, LL(c)));
        e[j].second -= a;
        e[j ^ 1].second += a;
        r -= a;
        if (r == 0) return f;
      }
    }
    return f - r;
  }
 public:
  Dinic(int n) : n_(n), g_(n) {}
  void addEdge(int u, int v, int c) {
    if (u == v) return;
    g_[u].emplace_back(e.size());
    e.emplace_back(v, c);
    g_[v].emplace_back(e.size());
    e.emplace_back(u, 0);
  }
  LL maxFlow(int s, int t) {
    LL r = 0;
    while (bfs(s, t)) {
      cur_.assign(n_, 0);
      r += dfs(s, t, INT64_MAX);
    }
    return r;
  }
};

// S-T max-Flow: HLPP $O(n_^2 \sqrt{m})$
class HLPP {
  int n_;
  // e[i] = {endPoint, conpacity} and e[i ^ 1] is opposite edge of e[i]
  // g_[u] = {edges start from u}
  std::vector<std::pair<int, int>> e;
  std::vector<std::vector<int>> g_;
  std::vector<int> h_;
  std::vector<LL> ex_;
  void addFlow(int i, int a) {
    ex_[e[i ^ 1].first] -= a;
    ex_[e[i].first] += a;
    e[i].second -= a;
    e[i ^ 1].second += a;
  };
  // d[u] = dist(u, t)
  bool init(int s, int t) {
    std::queue<int> Q;
    Q.push(t);
    h_[t] = 0;
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      for (auto i : g_[u]) {
        int v = e[i].first;
        if (e[i ^ 1].second > 0 && h_[v] == n_) {
          h_[v] = h_[u] + 1;
          Q.push(v);
        }
      }
    }
    return h_[s] == n_;
  }
 public:
  HLPP(int n) : n_(n), g_(n), h_(n, n), ex_(n) {}
  void addEdge(int u, int v, int c) {
    if (u == v) return;
    g_[u].emplace_back(e.size());
    e.emplace_back(v, c);
    g_[v].emplace_back(e.size());
    e.emplace_back(u, 0);
  }
  LL maxFlow(int s, int t) {
    if (init(s, t)) return 0;
    std::vector<int> gap(n_ + 1, 0), vis(n_);
    for (auto x : h_) ++gap[x];
    std::priority_queue<std::pair<int, int>> pq;
    // overload if ex_[u] > 0 after push, so lift height is needed.
    auto push = [&](int u) -> bool {
      if (ex_[u] == 0 || h_[u] == n_) return false;
      for (auto i : g_[u]) {
        auto [v, c] = e[i];
        if (c == 0 || (h_[u] != h_[v] + 1 && u != s)) continue;
        int a = std::min(ex_[u], LL(c));
        addFlow(i, a);
        if (!vis[v]) {
          pq.push({h_[v], v});
          vis[v] = 1;
        }
        if (ex_[u] == 0) return false;
      }
      return true;
    };
    ex_[s] = INT64_MAX;
    push(s);
    h_[s] = n_;
    vis[s] = vis[t] = 1;
    while (!pq.empty()) {
      int u = pq.top().second;
      pq.pop();
      vis[u] = 0;
      while (push(u)) {
        if (--gap[h_[u]] == 0) {
          for (int i = 0; i < n_; ++i) if (h_[i] > h_[u]) h_[i] = n_;
        }
        h_[u] = n_ - 1;
        for (auto i : g_[u]) {
          auto [v, c] = e[i];
          if (c > 0 && h_[u] > h_[v]) h_[u] = h_[v];
        }
        ++gap[++h_[u]];
      }
    }
    return ex_[t];
  }
};
// https://vjudge.net/problem/LibreOJ-127

// Global minimum cut of undirected graph: Stoer-Wagner $O(n_^3)$ implement
class StoerWagner {
  int n_;
  std::vector<std::vector<int>> g_;
  std::vector<int> del;
  void merge(int s, int t) {
    del[s] = 1;
    for (int i = 0; i < n_; ++i) {
      g_[i][t] = (g_[t][i] += g_[s][i]);
    }
  }
 public:
  StoerWagner(int n) : n_(n), g_(n, std::vector<int>(n)), del(n) {}
  void addEdge(int u, int v, int c) {
    if (u == v) return;
    g_[u][v] += c;
    g_[v][u] += c;
  }
  // the graph will be destory after minCut
  int minCut() {
    auto f = [&](int cnt, int& s, int& t) -> int {
      std::vector<int> vis(n_), d(n_);
      auto push = [&](int x){
        vis[x] = 1;
        d[x] = 0;
        for (int i = 0; i < n_; ++i) if (!del[i] && !vis[i]) d[i] += g_[x][i];
      };
      for (int i = 0; i < cnt; ++i) {
        push(t);
        s = t;
        t = std::max_element(d.begin(), d.end()) - d.begin();
      }
      return d[t];
    };
    int s = 0, t = 0, r = INT_MAX;
    for (int i = n_ - 1; i > 0; --i) {
      r = std::min(r, f(i, s, t));
      merge(s, t);
    }
    return r == INT_MAX ? 0 : r;
  }
};
// https://www.luogu.com.cn/problem/P5632

// Minimum cost maximum flow

class Flow {
  static inline constexpr LL INF = INT64_MAX >> 1;
  int n_;
  // e_[i] = {endPoint, conpacity} and e_[i ^ 1] is opposite edge of e_[i]
  // g_[u] = {edges start from u}
  std::vector<std::tuple<int, int, int>> e_;
  std::vector<std::vector<int>> g_;
  std::vector<int> path_;
  std::vector<LL> h_;
  // h_[i] = dist(s, i), h_[t] != -1 means there is a path_ from s to t, and h_[t] will be the potential.
  // path_[v]: short path_ form s to v, path_[v] is the previous node of v
  bool Dijkstra(int s, int t) {
    std::priority_queue<std::pair<LL, int>> Q;
    std::fill(path_.begin(), path_.end(), -1);
    std::vector<LL> d(n_, INF);
    d[s] = 0;
    Q.push({0, s});
    while (!Q.empty()) {
      auto [du, u] = Q.top();
      Q.pop();
      if (d[u] != -du) continue;
      for (auto i : g_[u]) {
        auto [v, w, c] = e_[i];
        c += h_[u] - h_[v];
        if (w > 0 && d[v] > d[u] + c) {
          d[v] = d[u] + c;
          path_[v] = i;
          Q.push({-d[v], v});
        }
      }
    }
    for (int i = 0; i < n_; ++i) {
      if ((h_[i] += d[i]) > INF) h_[i] = INF;
    }
    return h_[t] != INF;
  }
 public:
  Flow(int n) : n_(n), g_(n), h_(n), path_(n) {}
  void addEdge(int u, int v, int w, int c) {
    if (u == v) return;
    g_[u].emplace_back(e_.size());
    e_.emplace_back(v, w, c);
    g_[v].emplace_back(e_.size());
    e_.emplace_back(u, 0, -c);
  }
  std::pair<LL, LL> maxFlow(int s, int t) {
    LL flow = 0, cost = 0;
    while (Dijkstra(s, t)) {
      int f = INT_MAX, now = t;
      std::vector<int> r;
      while (now != s) {
        r.emplace_back(path_[now]);
        f = std::min(f, std::get<1>(e_[path_[now]]));
        now = std::get<0>(e_[path_[now] ^ 1]);
      }
      for (auto i : r) {
        std::get<1>(e_[i]) -= f;
        std::get<1>(e_[i ^ 1]) += f;
      }
      flow += f;
      cost += LL(f) * h_[t];
    }
    return {flow, cost};
  }
};
// https://www.luogu.com.cn/problem/P3381

// $O(m \sqrt{m})$, we will get TLE if the answer greater than INT_MAX
int circle3count(const std::vector<std::pair<int, int>>& edge, int n) {
  std::vector<int> d(n), vis(n, -1);
  for (auto [u, v] : edge) ++d[u], ++d[v];
  std::vector<std::vector<int>> e(n);
  // Giving Orienting to Edge
  for (auto [u, v] : edge) {
    if (d[u] < d[v] || (d[u] == d[v] && u < v)) {
      e[u].emplace_back(v);
    } else {
      e[v].emplace_back(u);
    }
  }
  int ans = 0;
  for (int i = 0; i < n; ++i) {
    for (auto u : e[i]) vis[u] = i;
    for (auto u : e[i]) {
      for (auto v : e[u]) if (vis[v] == i) ++ans;
    }
  }
  return ans;
}
// https://www.luogu.com.cn/problem/P1989


// $O(m \sqrt{m})$
LL circle4count(const std::vector<std::pair<int, int>>& edge, int n) {
  std::vector<int> d(n), c(n, -1), id(n);
  for (auto [u, v] : edge) ++d[u], ++d[v];
  std::iota(id.begin(), id.end(), 0);
  std::sort(id.begin(), id.end(), [&](int i, int j) {
    return d[i] < d[j] || (d[i] == d[j] && i < j);
  });
  std::vector<std::vector<int>> e(n);
  for (auto [u, v] : edge) {
    e[u].emplace_back(v);
    e[v].emplace_back(u);
  }
  // x -> y -> z and x -> w -> z
  LL ans = 0;
  for (int i = 0; i < n; ++i) {
    for (auto u : e[i]) if (id[i] < id[u]) {
      for (auto v : e[u]) if (id[i] < id[v]) ans += c[v]++;
    }
    for (auto u : e[i]) if (id[i] < id[u]) {
      for (auto v : e[u]) if (id[i] < id[v]) c[v] = 0;
    }
  }
  return ans;
}
// https://www.luogu.com.cn/blog/221955/san-yuan-huan-si-yuan-huan-ji-shuo
