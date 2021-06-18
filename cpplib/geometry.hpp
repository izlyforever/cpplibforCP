#pragma once
#include <bits/stdc++.h>

namespace Geomerty {
using Point = std::pair<double, double>;
bool crossLeft(const Point &op, const Point &sp, const Point &ep) {
	return (sp.first - op.first) * (ep.second - op.second) 
	<= (sp.second - op.second) * (ep.first - op.first);
}
double cross(const Point &op, const Point &sp, const Point &ep) {
	return (sp.first - op.first) * (ep.second - op.second) 
	- (sp.second - op.second) * (ep.first - op.first);	
}
double dist2(const Point &p, const Point &q) {
	double x = q.first - p.first, y = q.second - p.second;
	return x * x + y * y;
};

std::vector<Point> convexHull(std::vector<Point> p) {
	std::sort(p.begin(), p.end());
	p.erase(std::unique(p.begin(), p.end()), p.end());
	int n = p.size();
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
	int n = q.size();
	q.emplace_back(q[0]);
	double ans = 0;
	for (int i = 0, j = 2; i < n; ++i) {
		while (cross(q[i], q[i + 1], q[j]) < cross(q[i], q[i + 1], q[j + 1])) j = (j + 1) % n;
		ans = std::max({ans, dist2(q[i], q[j]), dist2(q[i + 1], q[j])});
	}
	return std::sqrt(ans);
} // float version: https://www.luogu.com.cn/problem/P6247
// Int version: https://www.luogu.com.cn/problem/P1452

double dist (const Point& p, const Point &q) {
	double x = q.first - p.first, y = q.second - p.second;
	return std::sqrt(x * x + y * y);
};

double minDist(std::vector<Point> a) {
	double d = DBL_MAX;
	int n = a.size();
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
		for (int i = 0, np = p.size(); i < np; ++i) {
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
std::vector<int> partialOrder(std::vector<std::vector<int>> &a) {
	int k = a.size(), n = a[0].size();
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
