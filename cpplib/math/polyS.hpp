#pragma once
#include <bits/stdc++.h>
using LL = long long;

class PolyS : public std::vector<int> {
	static inline std::vector<int> rev, roots{0, 1};
	static int powMod(int x, int n) {
		int r(1);
		while (n) {
			if (n&1) r = 1LL * r * x % M;
			n >>= 1; x = 1LL * x * x % M;
		}
		return r;
	}
	void dft() {
		int n = size();
		if ((int)rev.size() != n) {
			int k = __builtin_ctz(n) - 1;
			rev.resize(n);
			for (int i = 0; i < n; ++i) {
				rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
			}
		}
		if ((int)roots.size() < n) {
			int k = __builtin_ctz(roots.size());
			roots.resize(n);
			while ((1 << k) < n) {
				int e = powMod(g, (M - 1) >> (k + 1));
				for (int i = 1 << (k - 1); i < (1 << k); ++i) {
					roots[2 * i] = roots[i];
					roots[2 * i + 1] = 1LL * roots[i] * e % M;
				}
				++k;
			}
		}
		for (int i = 0; i < n; ++i) if (rev[i] < i) {
			std::swap((*this)[i], (*this)[rev[i]]);
		}
		for (int k = 1; k < n; k *= 2) {
			for (int i = 0; i < n; i += 2 * k) {
				for (int j = 0; j < k; ++j) {
					int u = (*this)[i + j];
					int v = 1LL * (*this)[i + j + k] * roots[k + j] % M;
					int x = u + v, y = u - v;
					if (x >= M) x -= M;
					if (y < 0) y += M;
					(*this)[i + j] = x;
					(*this)[i + j + k] = y;
				}
			}
		}
	}
	void idft() {
		int n = size();
		std::reverse(begin() + 1, end());
		dft();
		int inv = powMod(n, M - 2);
		for (int i = 0; i < n; ++i) {
			(*this)[i] = 1LL * (*this)[i] * inv % M;
		}
	}
	void standard() {
		while (!empty() && !back()) pop_back();
	}
	void reverse() {
		std::reverse(begin(), end());
		standard();
	}
public:
	static inline const int M = 998244353, g = 3;
	PolyS() {}
	PolyS(const std::vector<int> &a) : std::vector<int>{a} { standard();}
	PolyS(const int &x) : std::vector<int>{x} { standard();}
	int at(int id) const {
		if (id < 0 || id >= (int)size()) return 0;
		return (*this)[id];
	}
	PolyS operator-() const {
		auto A = *this;
		for (auto &x : A) x = (x == 0 ? 0 : M - x);
		return A;
	}	
	PolyS mulXn(int n) const {
		auto A = *this;
		if (!A.empty()) A.insert(A.begin(), n, 0);
		return A;
	}
	PolyS modXn(int n) const {
		if (n > (int)size()) return *this;
		return PolyS({begin(), begin() + n});
	}
	PolyS divXn(int n) const {
		if ((int)size() <= n) return PolyS();
		return PolyS({begin() + n, end()});
	}
	PolyS &operator+=(const PolyS &rhs) {
		if ((int)size() < (int)rhs.size()) resize(rhs.size());
		for (int i = 0, rs = rhs.size(); i < rs; ++i) {
			if (((*this)[i] += rhs[i]) >= M) (*this)[i] -= M;
		}
		standard();
		return *this;
	}
	PolyS &operator-=(const PolyS &rhs) {
		if (size() < rhs.size()) resize(rhs.size());
		for (int i = 0, rs = rhs.size(); i < rs; ++i) {
			if (((*this)[i] -= rhs[i]) < 0) (*this)[i] += M;
		}
		standard();
		return *this;
	}
	PolyS &operator*=(PolyS rhs) {
		int n = size(), m = rhs.size(), tot = std::max(1, n + m - 1);
		int sz = 1 << std::__lg(tot * 2 - 1);
		resize(sz);
		rhs.resize(sz);
		dft();
		rhs.dft();
		for (int i = 0; i < sz; ++i) {
			(*this)[i] = 1LL * (*this)[i] * rhs[i] % M;
		}
		idft();
		standard();
		return *this;
	}
	PolyS &operator/=(PolyS rhs) {
		int n = size(), m = rhs.size();
		if (n < m) return (*this) = PolyS();
		reverse();
		rhs.reverse();
		(*this) *= rhs.inv(n - m + 1);
		resize(n - m + 1);
		reverse();
		return *this;
	}
	PolyS &operator%=(PolyS rhs) {
		return (*this) -= (*this) / rhs * rhs; 
	}
	PolyS operator+(const PolyS &rhs) const {
		return PolyS(*this) += rhs;
	}
	PolyS operator-(const PolyS &rhs) const {
		return PolyS(*this) -= rhs;
	}
	PolyS operator*(PolyS rhs) const {
		return PolyS(*this) *= rhs;
	}
	PolyS operator/(PolyS rhs) const {
		return PolyS(*this) /= rhs;
	}
	PolyS operator%(PolyS rhs) const {
		return PolyS(*this) %= rhs;
	}
	PolyS powModPoly(int n, PolyS p) {
		PolyS r(1), x(*this);
		while (n) {
			if (n&1) (r *= x) %= p;
			n >>= 1; (x *= x) %= p;
		}
		return r;
	}
	int inner(const PolyS &rhs) {
		int r = 0, n = std::min(size(), rhs.size());
		for (int i = 0; i < n; ++i) {
			r = (r + 1LL * (*this)[i] * rhs[i]) % M;
		}
		return r;
	}
	PolyS derivation() const {
		if (empty()) return PolyS();
		int n = size();
		std::vector<int> r(n - 1);
		for (int i = 1; i < n; ++i) r[i - 1] =  1LL * (*this)[i] * i % M;
		return PolyS(r);
	}
	PolyS integral() const {
		if (empty()) return PolyS();
		int n = size();
		std::vector<int> r(n + 1), inv(n + 1, 1);
		for (int i = 2; i <= n; ++i) inv[i] = 1LL * (M - M / i) * inv[M % i] % M;
		for (int i = 0; i < n; ++i) r[i + 1] = 1LL * (*this)[i] * inv[i + 1] % M;
		return PolyS(r);
	}
	PolyS inv(int n) const {
		// assert((*this)[0] != 0);
		PolyS x(powMod((*this)[0], M - 2));
		int k = 1;
		while (k < n) {
			k *= 2;
			x *= (PolyS(2) - modXn(k) * x).modXn(k);
		}
		return x.modXn(n);
	}
	// assume a[0] = 1
	PolyS log(int n) const {
		return (derivation() * inv(n)).integral().modXn(n);
	}
	// assume a[0] = 0
	PolyS exp(int n) const {
		PolyS x(1);
		int k = 1;
		while (k < n) {
			k *= 2;
			x = (x * (PolyS(1) - x.log(k) + modXn(k))).modXn(k);
		}
		return x.modXn(n);
	}
	// assume a[0] = 1;
	PolyS sqrt(int n) const {
		const static int inv2 = (M + 1) / 2;
		PolyS x(1);
		int k = 1;
		while (k < n) {
			k *= 2;
			x += modXn(k) * x.inv(k);
			x = x.modXn(k) * inv2;
		}
		return x.modXn(n);
	}
	// transpose convolution {\rm MULT}(F(x),G(x))=\sum_{i\ge0}(\sum_{j\ge 0}f_{i+j}g_j)x^i
	PolyS mulT(PolyS rhs) const {
		if (rhs.size() == 0) return PolyS();
		int n = rhs.size();
		std::reverse(rhs.begin(), rhs.end());
		return ((*this) * rhs).divXn(n - 1);
	}
	int eval(int x) {
		int r = 0, t = 1;
		for (int i = 0, n = size(); i < n; ++i) {
			r = (r + 1LL * (*this)[i] * t) % M;
			t = 1LL * t * x % M;
		}
		return r;
	}
	// multi-evaluation(new tech)
	std::vector<int> evals(std::vector<int> x) const {
		if (size() == 0) return std::vector<int>(x.size());
		int n = x.size();
		std::vector<int> ans(n);
		std::vector<PolyS> g(4 * n);
		std::function<void(int, int, int)> build = [&](int l, int r, int p) {
			if (r - l == 1) {
				// g[p] = std::vector<int>{1, x[i] ? M - x[l] : 0};
				g[p] = PolyS({1, x[l] ? M - x[l] : 0});
			} else {
				int m = (l + r) / 2;
				build(l, m, 2 * p);
				build(m, r, 2 * p + 1);
				g[p] = g[2 * p] * g[2 * p + 1];
			}
		};
		build(0, n, 1);
		std::function<void(int, int, int, const PolyS &)> solve = [&](int l, int r, int p, const PolyS &f) {
			if (r - l == 1) {
				ans[l] = f.at(0);
			} else {
				int m = (l + r) / 2;
				solve(l, m, 2 * p, f.mulT(g[2 * p + 1]).modXn(m - l));
				solve(m, r, 2 * p + 1, f.mulT(g[2 * p]).modXn(r - m));
			}
		};
		solve(0, n, 1, mulT(g[1].inv(size())).modXn(n));
		return ans;
	} // https://www.luogu.com.cn/problem/P5050
}; // https://www.luogu.com.cn/training/3015#information