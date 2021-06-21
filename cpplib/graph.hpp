#pragma once
#include <bits/stdc++.h>

using edge = std::vector<std::pair<int, int>>;
using Edge = std::tuple<int, int, int>;

// dfs order of rooted tree
class DfsTour {
	int n, cnt;
	std::vector<int> l, r;
	std::vector<std::vector<int>> e;
public:
	DfsTour(int _n) : n(_n), cnt(0), l(n), r(n), e(n) {}
	void addEdge(int u, int v) {
		if (u == v) return;
		e[u].emplace_back(v);
		e[v].emplace_back(u);
	}
	void dfs(int u, int fa) {
		l[u] = ++cnt;
		for (auto v : e[u]) if (v != fa) {
			dfs(v, u);
		}
		r[u] = cnt;
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
	int n;
	std::vector<int> fa, dep, sz, son, top;
public:
	LCA(std::vector<std::vector<int>> &e, int rt = 1) : n(e.size()) {
		fa.resize(n);
		dep.resize(n);
		sz.resize(n);
		son.resize(n);
		fa[rt] = rt;
		dep[rt] = 0;
		std::function<int(int)> pdfs = [&](int u) -> int {
			sz[u] = 1;
			for (auto v : e[u]) if (v != fa[u]) {
				dep[v] = dep[u] + 1;
				fa[v] = u;
				sz[u] += pdfs(v);
				if (sz[v] > sz[son[u]]) son[u] = v;
			}
			return sz[u];
		};
		top.resize(n);
		std::function<void(int, int)> dfs = [&](int u, int t) -> void {
			top[u] = t;
			if (son[u] == 0) return;
			dfs(son[u], t);
			for (auto v : e[u]) if (v != fa[u] && v != son[u]) dfs(v, v);
		};
		pdfs(rt);
		dfs(rt, rt);
	}
	int lca(int u, int v) {
		while (top[u] != top[v]) {
			if (dep[top[u]] > dep[top[v]]) {
				u = fa[top[u]];
			} else {
				v = fa[top[v]];
			}
		}
		return dep[u] < dep[v] ? u : v;
	}
};

// Minimum Spanning Tree
LL Prim(const std::vector<edge> &e) {
	LL r = 0;
	int n = e.size(), cnt = 0;
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
			// note there may be a path of form '6'
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
		for (auto &[u, v, w] : e) {
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
std::vector<int> TopSort(std::vector<std::set<int>> &e) {
	std::vector<int> r;
	std::priority_queue<int> Q;
	int n = e.size();
	std::vector<int> in(n);
	for (auto &x : e) for (auto i : x) ++in[i];
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
	for (int i = 0, ne = e.size(); i < ne; ++i) {
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

std::vector<LL> Dijkstra(int s, const std::vector<edge> &e) {
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

bool BellmanFord(std::vector<Edge> &e, int n, int x = 0) {
	std::vector<int> dist(n + 1, INT_MAX);
	dist[x] = 0;
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

bool spfa(std::vector<edge> &e, int x = 0) {
	int n = e.size();
	std::queue<int> Q;
	std::vector<int> dist(n, INT_MAX), cnt(n), inQ(n);
	Q.push(x);
	inQ[x] = 1;
	dist[x] = 0;
	++cnt[x];
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
struct Scc {
	int n, nScc;
	std::vector<int> vis, color, order;
	std::vector<std::vector<int>> e, e2;
	Scc(int _n) : n(_n * 2) {
		nScc = 0;
		e.resize(n);
		e2.resize(n);
		vis.resize(n);
		color.resize(n);
	}
	void addEdge(int u, int v) {
		e[u].emplace_back(v);
		e2[v].emplace_back(u);
	}
	void dfs(int u) {
		vis[u] = true;
		for (auto v : e[u]) if (!vis[v]) dfs(v);
		order.emplace_back(u);
	}
	void dfs2(int u) {
		color[u] = nScc;
		for (auto v : e2[u]) if (!color[v]) dfs2(v);
	}
	void Kosaraju() {
		for (int i = 0; i < n; ++i) if (!vis[i]) dfs(i);
		for (auto it = order.rbegin(); it != order.rend(); ++it) if (!color[*it]) {
			++nScc;
			dfs2(*it);
		}
	}
};

// n / 2 pairs (2i, 2i + 1)
struct twoSAT {
	int n, nScc;
	std::vector<int> vis, color, order;
	std::vector<std::vector<int>> e, e2;
	twoSAT(int _n) : n(_n * 2) {
		nScc = 0;
		e.resize(n);
		e2.resize(n);
		vis.resize(n);
		color.resize(n);
	}
	void addEdge(int u, int v) {
		e[u].emplace_back(v);
		e2[v].emplace_back(u);
	}
	void dfs(int u) {
		vis[u] = true;
		for (auto v : e[u]) if (!vis[v]) dfs(v);
		order.emplace_back(u);
	}
	void dfs2(int u) {
		color[u] = nScc;
		for (auto v : e2[u]) if (!color[v]) dfs2(v);
	}
	void Kosaraju() {
		for (int i = 0; i < n; ++i) if (!vis[i]) dfs(i);
		for (auto it = order.rbegin(); it != order.rend(); ++it) if (!color[*it]) {
			++nScc;
			dfs2(*it);
		}
	}
	std::vector<int> solve() {
		Kosaraju();
		// choose component with max color number
		std::vector<int> choose(nScc + 1);
		for (int i = 0; i < n; i += 2) {
			int c1 = color[i], c2 = color[i + 1];
			if (c1 == c2) return std::vector<int>();
			if (choose[c1] || choose[c2]) continue;
			choose[std::max(c1, c2)] = 1;
		}
		std::vector<int> r(n / 2);
		for (int i = 0; i * 2 < n; ++i) r[i] = (choose[color[i * 2]] ? 1 : -1); 
		return r;
	}
};

std::vector<int> cutVertex(std::vector<std::vector<int>>& e) {
	int n = e.size(), cnt = 0;
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
	int n, cnt;
	std::vector<std::vector<int>> g;
	std::vector<int> e, flag, dfs, low;
	void Tarjan(int u, int inEdgeNum) {
		low[u] = dfs[u] = ++cnt;
		for (auto i : g[u]) {
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
	CutEdge(int _n) : n(_n), cnt(0), g(n), dfs(n), low(n) {}
	void addEdge(int u, int v) {
		if (u == v) return;
		g[u].emplace_back(e.size());
		e.emplace_back(v);
		flag.emplace_back(0);
		g[v].emplace_back(e.size());
		e.emplace_back(u);
		flag.emplace_back(0);
	}
	int solve() {
		for (int i = 0; i < n; ++i) if (dfs[i] == 0) Tarjan(i, -1);
		int r = 0;
		for (auto x : flag) r += x;
		return r / 2;
	}
};
// https://www.luogu.com.cn/problem/T103481

// S-T max-Flow: Dinic $O(n^2 m)$
class Dinic {
	int n;
	// e[i] = {endPoint, conpacity} and e[i ^ 1] is opposite edge of e[i]
	// g[u] = {edges start from u}
	std::vector<std::pair<int, int>> e;
	std::vector<std::vector<int>> g;
	std::vector<int> cur, h;
	// h[i] = dist(s, i), so h[t] != -1 means there is a path from s to t
	bool bfs(int s, int t) {
		h.assign(n, -1);
		std::queue<int> Q;
		h[s] = 0;
		Q.push(s);
		while (!Q.empty()) {
			int u = Q.front();
			Q.pop();
			for (auto i : g[u]) {
				auto [v, c] = e[i];
				if (c > 0 && h[v] == -1) {
					h[v] = h[u] + 1;
					Q.push(v);
				}
			}
		}
		return h[t] != -1;
	}
	LL dfs(int u, int t, LL f) {
		if (u == t || f == 0) return f;
		LL r = f;
		for (int &i = cur[u], ng = g[u].size(); i < ng; ++i) {
			int j = g[u][i];
			auto [v, c] = e[j];
			if (c > 0 && h[v] == h[u] + 1) {
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
	Dinic(int _n) : n(_n), g(n) {}
	void addEdge(int u, int v, int c) {
		if (u == v) return;
		g[u].emplace_back(e.size());
		e.emplace_back(v, c);
		g[v].emplace_back(e.size());
		e.emplace_back(u, 0);
	}
	LL maxFlow(int s, int t) {
		LL r = 0;
		while (bfs(s, t)) {
			cur.assign(n, 0);
			r += dfs(s, t, INT64_MAX);
		}
		return r;
	}
};

// S-T max-Flow: HLPP $O(n^2 \sqrt{m})$
class HLPP {
	int n;
	// e[i] = {endPoint, conpacity} and e[i ^ 1] is opposite edge of e[i]
	// g[u] = {edges start from u}
	std::vector<std::pair<int, int>> e;
	std::vector<std::vector<int>> g;
	std::vector<int> h;
	std::vector<LL> ex;
	void addFlow(int i, int a) {
		ex[e[i ^ 1].first] -= a;
		ex[e[i].first] += a;
		e[i].second -= a;
		e[i ^ 1].second += a;
	};
	// d[u] = dist(u, t)
	bool init(int s, int t) {
		std::queue<int> Q;
		Q.push(t);
		h[t] = 0;
		while (!Q.empty()) {
			int u = Q.front();
			Q.pop();
			for (auto i : g[u]) {
				int v = e[i].first;
				if (e[i ^ 1].second > 0 && h[v] == n) {
					h[v] = h[u] + 1;
					Q.push(v);
				}
			}
		}
		return h[s] == n;
	}
public:
	HLPP(int _n) : n(_n), g(n), h(n, n), ex(n) {}
	void addEdge(int u, int v, int c) {
		if (u == v) return;
		g[u].emplace_back(e.size());
		e.emplace_back(v, c);
		g[v].emplace_back(e.size());
		e.emplace_back(u, 0);
	}
	LL maxFlow(int s, int t) {
		if (init(s, t)) return 0;
		std::vector<int> gap(n + 1, 0), vis(n);
		for (auto x : h) ++gap[x];
		std::priority_queue<std::pair<int, int>> pq;
		// overload if ex[u] > 0 after push, so lift height is needed.
		auto push = [&](int u) -> bool {
			if (ex[u] == 0 || h[u] == n) return false;
			for (auto i : g[u]) {
				auto [v, c] = e[i];
				if (c == 0 || (h[u] != h[v] + 1 && u != s)) continue;
				int a = std::min(ex[u], LL(c));
				addFlow(i, a);
				if (!vis[v]) {
					pq.push({h[v], v});
					vis[v] = 1;
				}
				if (ex[u] == 0) return false;
			}
			return true;
		};
		ex[s] = INT64_MAX;
		push(s);
		h[s] = n;
		vis[s] = vis[t] = 1;
		while (!pq.empty()) {
			int u = pq.top().second;
			pq.pop();
			vis[u] = 0;
			while (push(u)) {
				if (--gap[h[u]] == 0) {
					for (int i = 0; i < n; ++i) if (h[i] > h[u]) h[i] = n;
				}
				h[u] = n - 1;
				for (auto i : g[u]) {
					auto [v, c] = e[i];
					if (c > 0 && h[u] > h[v]) h[u] = h[v];
				}
				++gap[++h[u]];
			}
		}
		return ex[t];
	}
};
// https://vjudge.net/problem/LibreOJ-127

// Global minimum cut of undirected graph: Stoer-Wagner $O(n^3)$ implement
class StoerWagner {
	int n;
	std::vector<std::vector<int>> g;
	std::vector<int> del;
	void merge(int s, int t) {
		del[s] = 1;
		for (int i = 0; i < n; ++i) {
			g[i][t] = (g[t][i] += g[s][i]);
		}
	}
public:
	StoerWagner(int _n) : n(_n), g(n, std::vector<int>(n)), del(n) {}
	void addEdge(int u, int v, int c) {
		if (u == v) return;
		g[u][v] += c;
		g[v][u] += c;
	}
	// the graph will be destory after minCut
	int minCut() {
		auto f = [&](int cnt, int &s, int &t) -> int {
			std::vector<int> vis(n), d(n);
			auto push = [&](int x){
				vis[x] = 1;
				d[x] = 0;
				for (int i = 0; i < n; ++i) if (!del[i] && !vis[i]) d[i] += g[x][i];
			};
			for (int i = 0; i < cnt; ++i) {
				push(t);
				s = t;
				t = std::max_element(d.begin(), d.end()) - d.begin();
			}
			return d[t];
		};
		int s = 0, t = 0, r = INT_MAX;
		for (int i = n - 1; i > 0; --i) {
			r = std::min(r, f(i, s, t));
			merge(s, t);
		}
		return r == INT_MAX ? 0 : r;
	}
};
// https://www.luogu.com.cn/problem/P5632

// Minimum cost maximum flow
class Flow {
	static inline const int INF = 1e9;
	int n;
	// e[i] = {endPoint, conpacity} and e[i ^ 1] is opposite edge of e[i]
	// g[u] = {edges start from u}
	std::vector<std::tuple<int, int, int>> e;
	std::vector<std::vector<int>> g;
	std::vector<int> h, path;
	// h[i] = dist(s, i), h[t] != -1 means there is a path from s to t, and h[t] will be the potential.
	// path[v]: short path form s to v, path[v] is the previous node of v
	bool Dijkstra(int s, int t) {
		std::priority_queue<std::pair<int, int>> Q;
		std::fill(path.begin(), path.end(), -1);
		std::vector<int> d(n, INF);
		d[s] = 0;
		Q.push({0, s});
		while (!Q.empty()) {
			auto [du, u] = Q.top();
			Q.pop();
			if (d[u] != -du) continue;
			for (auto i : g[u]) {
				auto [v, w, c] = e[i];
				c += h[u] - h[v];
				if (w > 0 && d[v] > d[u] + c) {
					d[v] = d[u] + c;
					path[v] = i;
					Q.push({-d[v], v});
				}
			}
		}
		for (int i = 0; i < n; ++i) {
			if ((h[i] += d[i]) > INF) h[i] = INF;
		}
		return h[t] != INF;
	}
public:
	Flow(int _n) : n(_n), g(n), h(n), path(n) {}
	void addEdge(int u, int v, int w, int c) {
		if (u == v) return;
		g[u].emplace_back(e.size());
		e.emplace_back(v, w, c);
		g[v].emplace_back(e.size());
		e.emplace_back(u, 0, -c);
	}
	std::pair<LL, LL> maxFlow(int s, int t) {
		LL flow = 0, cost = 0;
		while (Dijkstra(s, t)) {
			int f = INT_MAX, now = t;
			std::vector<int> r;
			while (now != s) {
				r.emplace_back(path[now]);
				f = std::min(f, std::get<1>(e[path[now]]));
				now = std::get<0>(e[path[now] ^ 1]);
			}
			for (auto i : r) {
				std::get<1>(e[i]) -= f;
				std::get<1>(e[i ^ 1]) += f;
			}
			flow += f;
			cost += LL(f) * h[t];
		}
		return {flow, cost};
	}
};
// https://www.luogu.com.cn/problem/P3381

// $O(m \sqrt{m})$, we will get TLE if the answer greater than INT_MAX
int circle3count(const std::vector<std::pair<int, int>> &edge, int n) {
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
LL circle4count(const std::vector<std::pair<int, int>> &edge, int n) {
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