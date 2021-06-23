// Main reference: https://www.luogu.com.cn/blog/command-block/sheng-cheng-han-shuo-za-tan
#pragma once
#include <bits/stdc++.h>
#include "basic.hpp"
using LL = long long;

// using valT = decltype(T::a)::value_type;
template<typename T, typename valT>
class Poly : public T {
	// many function will fail for the case n > mod
	static inline const valT j = pow(valT(3), (valT::mod() - 1) / 4);
	static inline const valT inv2 = (valT::mod() + 1) / 2;
	static inline constexpr int maxN = 1e6 + 2;  // assume size(a) < maxN
	static inline const auto &Binom = BinomModp<valT>::Instance(maxN);
	// static inline constexpr Binom
public:
	using T::T;
	Poly (const T &x) : T(x) {}
	Poly mulXn(int n) const {
		auto b = this->a;
		b.insert(b.begin(), n, 0);
		return Poly(b);
	}
	Poly modXn(int n) const {
		if (n > (int)this->size()) return *this;
		return Poly({this->a.begin(), this->a.begin() + n});
	}
	Poly divXn(int n) const {
		if ((int)this->size() <= n) return Poly();
		return Poly({this->a.begin() + n, this->a.end()});
	}
	Poly operator-() const {
		auto A = *this;
		for (auto &x : A.a) x = -x;
		return A;
	}
	Poly &operator+=(const Poly &rhs) {
		if (this->size() < rhs.size()) this->a.resize(rhs.size());
		for (int i = 0, nr = rhs.size(); i < nr; ++i) this->a[i] += rhs.a[i];
		this->standard();
		return *this;
	}
	Poly &operator-=(const Poly &rhs) {
		if (this->size() < rhs.size()) this->a.resize(rhs.size());
		for (int i = 0, nr = rhs.size(); i < nr; ++i) this->a[i] -= rhs.a[i];
		this->standard();
		return *this;
	}
	Poly operator+(const Poly &rhs) const {
		return Poly(*this) += rhs;
	}
	Poly operator-(const Poly &rhs) const {
		return Poly(*this) -= rhs;
	}
	Poly operator*(const Poly &rhs) const {
		return this->mul(rhs);
	}
	Poly &operator*=(const Poly &rhs) {
		return (*this) = (*this) * rhs;
	}
	// assume a[0] \neq 0
	Poly inv(int n) const {
		// assert(this->a[0] != 0);
		Poly x(this->a[0].inv());
		int k = 1;
		while (k < n) {
			k *= 2;
			x *= (Poly(2) - this->modXn(k) * x).modXn(k);
		}
		return x.modXn(n);
	}
	Poly &operator/=(Poly rhs) {
		int n = this->size(), m = rhs.size();
		if (n < m) return (*this) = Poly();
		this->reverse();
		rhs.reverse();
		(*this) *= rhs.inv(n - m + 1);
		this->a.resize(n - m + 1);
		this->reverse();
		return *this;
	}
	Poly operator/(const Poly &rhs) const {
		return Poly(*this) /= rhs;
	}
	Poly &operator%=(const Poly &rhs) {
		return *this -= (*this) / rhs * rhs; 
	}
	Poly operator%(const Poly &rhs) const {
		return Poly(*this) %= rhs;
	}
	Poly powModPoly(LL n, const Poly &p) const {
		Poly r(1), x(*this);
		while (n) {
			if (n&1) r = r * x % p;
			n >>= 1; x = x * x % p;
		}
		return r;
	}
	valT inner(const Poly &rhs) const {
		valT r(0);
		int n = std::min(this->size(), rhs.size());
		for (int i = 0; i < n; ++i) r += this->a[i] * rhs.a[i];
		return r;
	}
	Poly derivation() const {
		if (this->a.empty()) return Poly();
		int n = this->size();
		std::vector<valT> r(n - 1);
		for (int i = 1; i < n; ++i) r[i - 1] = this->a[i] * valT(i);
		return Poly(r);
	}
	Poly integral() const {
		if (this->a.empty()) return Poly();
		int n = this->size();
		std::vector<valT> r(n + 1), inv(n + 1, 1);
		for (int i = 2; i <= n; ++i) inv[i] = valT(valT::mod() - valT::mod() / i) * inv[valT::mod() % i];
		for (int i = 0; i < n; ++i) r[i + 1] = this->a[i] * inv[i + 1];
		return Poly(r);
	}
	// assume a[0] = 1
	Poly log(int n) const {
		return (derivation() * inv(n)).integral().modXn(n);
	}
	// assume a[0] = 0
	Poly exp(int n) const {
		Poly x(1);
		int k = 1;
		while (k < n) {
			k *= 2;
			x = (x * (Poly(1) - x.log(k) + this->modXn(k))).modXn(k);
		}
		return x.modXn(n);
	}
	Poly sin(int n) const {
		auto A = *this;
		for (auto &x : A.a) x *= j;
		A = A.exp(n);
		A -= A.inv(n);
		auto m = -j * inv2;
		for (auto &x : A.a) x *= m;
		return A;
	}
	Poly cos(int n) const {
		auto A = *this;
		for (auto &x : A.a) x *= j;
		A = A.exp(n);
		A += A.inv(n);
		for (auto &x : A.a) x *= inv2;
		return A;
	}
	// assume a[0] = 0
	Poly asin(int n) const { // unkown for acos 
		auto D = this -> derivation();
		auto A = (Poly(1) - (*this) * (*this)).modXn(n - 1);
		D = D * A.sqrt(n - 1).inv(n - 1);
		D = D.modXn(n - 1);
		return D.integral();
	}
	// assume a[0] = 0
	Poly atan(int n) const {
		auto D = this -> derivation();
		auto A = (Poly(1) + (*this) * (*this)).modXn(n - 1);
		D = D * A.inv(n - 1);
		D = D.modXn(n - 1);
		return D.integral();
	}
	// assum a[0] = 1
	Poly sqrt(int n) const {
		Poly x(1);
		int k = 1;
		while (k < n) {
			k *= 2;
			x += this->modXn(k) * x.inv(k);
			x = x.modXn(k) * inv2;
		}
		return x.modXn(n);
	}
	// transpose {\rm MULT}(F(x),G(x))=\sum_{i\ge0}(\sum_{j\ge 0}f_{i+j}g_j)x^i
	Poly mulT(Poly rhs) const {
		if (rhs.size() == 0) return Poly();
		int n = rhs.size();
		std::reverse(rhs.a.begin(), rhs.a.end());
		return ((*this) * rhs).divXn(n - 1);
	}
	// compose poly common algorithm F(A(x)) in $O(n^2)$, however Brent-Kung algorithm with $(n \log n)^{1.5}$ may be slower.
	Poly compose(Poly A, int n) const {
		A = A.modXn(n);
		int sn = std::sqrt(n);
		std::vector<Poly> G(sn);
		G[0] = {1};
		for (int i = 1; i < sn; ++i) G[i] = (G[i - 1] * A).modXn(n);
		auto B = (G.back() * A).modXn(n);
		Poly now{1}, ans;
		for (int i = 0; i < n; i += sn) {
			Poly sm;
			for (int j = 0; j < sn && i + j < n; ++j) {
				auto m = this->at(i + j);
				auto tmp = G[j];
				for (auto &x : tmp.a) x *= m;
				sm += tmp;
			}
			ans += (now * sm).modXn(n);
			now = (now * B).modXn(n);
		}
		return ans;
	}
	// compose inverse (assmue a[0] = 0, a[1] \neq 0) based on Langrange: $[x^n]G(x) = \frac{1}{n}[x^{n - 1}](\frac{x}{F(x)})^n$
	Poly composeInv(int n) const {
		auto A = this->divXn(1).inv(n - 1);
		int sn = std::sqrt(n);
		std::vector<Poly> G(sn + 1);
		G[0] = {1};
		for (int i = 1; i <= sn; ++i) G[i] = (G[i - 1] * A).modXn(n - 1);
		std::vector<valT> ans(n), inv(n);
		auto M = valT::mod();
		inv[1] = 1;
		for (int i = 2; i < n; ++i) inv[i] = inv[M % i] * valT(M - M / i);
		Poly now{1};
		for (int i = 0; i < n; i += sn) {
			for (int j = 1; j <= sn && i + j < n; ++j) {
				valT tmp;
				auto &sg = G[j];
				for (int k = 0, sk = i + j - 1; k <= sk; ++k) {
					tmp += now.at(k) * sg.at(sk - k);
				}
				ans[i + j] = tmp * inv[i + j];
			}
			now = (now * G.back()).modXn(n - 1);
		}
		return Poly(ans);
	}
	Poly toFallingPowForm() {
		int n = this->size();
		std::vector<valT> x(n);
		for (int i = 0; i < n; ++i) x[i] = i;
		auto y = this->evals(x);
		auto tmp = Binom.ifac;
		for (int i = 1; i < n; i += 2) tmp[i] = -tmp[i];
		Poly A = Poly(y) * Poly(tmp);
		return A.modXn(n);
	}
	Poly fromFallingPowForm() {
		int n = this->size();
		Poly A = ((*this) * Poly(Binom.ifac)).modXn(n);
		std::vector<valT> x(n), y = A.a;
		for (int i = 0; i < n; ++i) x[i] = i;
		y.resize(n);
		for (int i = 0; i < n; ++i) y[i] *= Binom.fac[i];
		return Lagrange(x, y);
	}
	valT eval(valT x) const {
		valT r(0), t(1);
		for (int i = 0, n = this->size(); i < n; ++i) {
			r += this->a[i] * t;
			t *= x;
		}
		return r;
	}
	// multi-evaluation(new tech)
	std::vector<valT> evals(std::vector<valT> x) const {
		if (this->size() == 0) return std::vector<valT>(x.size());
		int n = x.size();
		std::vector<valT> ans(n);
		std::vector<Poly> g(4 * n);
		std::function<void(int, int, int)> build = [&](int l, int r, int p) {
			if (r - l == 1) {
				g[p] = Poly({1, -x[l]});
			} else {
				int m = (l + r) / 2;
				build(l, m, 2 * p);
				build(m, r, 2 * p + 1);
				g[p] = g[2 * p] * g[2 * p + 1];
			}
		};
		build(0, n, 1);
		std::function<void(int, int, int, const Poly &)> solve = [&](int l, int r, int p, const Poly &f) {
			if (r - l == 1) {
				ans[l] = f.at(0);
			} else {
				int m = (l + r) / 2;
				solve(l, m, 2 * p, f.mulT(g[2 * p + 1]).modXn(m - l));
				solve(m, r, 2 * p + 1, f.mulT(g[2 * p]).modXn(r - m));
			}
		};
		solve(0, n, 1, mulT(g[1].inv(this->size())).modXn(n));
		return ans;
	} // https://www.luogu.com.cn/problem/P5050

	// \sum_{i = 0}^{n - 1} a_i / (1 - b_i x)
	static std::vector<valT> sumFraction(std::vector<valT> a, std::vector<valT> b, int N) {
		std::function<std::pair<Poly, Poly>(int, int)> solve = [&](int l, int r) -> std::pair<Poly, Poly> {
			if (r - l == 1) return {Poly(a[l]), Poly({1, - b[l]})};
			int m = (l + r) / 2;
			auto [pl, ql] = solve(l, m);
			auto [pr, qr] = solve(m, r);
			return {pl * qr + pr * ql, ql * qr};
		};
		auto [p, q] = solve(0, a.size());
		p *= q.inv(N);
		auto ans = p.a;
		ans.resize(N);
		return ans;
	} // https://codeforces.com/gym/102978/problem/D
	
	// compute $h(m), \cdots, h(m + cnt - 1)$ accroding to $h(0), h(1), \cdots, h(d)$
	static std::vector<valT> valToVal(std::vector<valT> h, valT m, int cnt) { // m > h.size()
		int d = h.size() - 1;
		for (int i = 0; i <= d; ++i) {
			h[i] *= Binom.ifac[i] * Binom.ifac[d - i];
			if ((d - i) & 1) h[i] = -h[i];
		}
		std::vector<valT> f(d + cnt);
		auto now = m - valT(d);
		for (int i = 0; i < d + cnt; ++i) {
			f[i] = now.inv();
			++now;
		}
		h = (Poly(f) * Poly(h)).a;
		h.resize(d + cnt);
		h = std::vector<valT>(h.begin() + d, h.end());
		now = 1;
		for (int i = 0; i <= d; ++i) now *= m - valT::raw(i);
		h[0] *= now;
		for (int i = 1; i < cnt; ++i) {
			now *= m + valT::raw(i);
			now *= (m + valT(i - d - 1)).inv();
			h[i] *= now;
		}
		return h;
	}; // https://www.luogu.com.cn/problem/P5667
	static Poly Lagrange(std::vector<valT> x, std::vector<valT> y) {
		std::function<Poly(int l, int r)> mulP = [&](int l, int r) -> Poly {
			if (r - l == 1) return Poly({-x[l], 1});
			int m = (l + r) / 2;
			return mulP(l, m) * mulP(m, r);
		};
		int n = x.size();
		auto A = mulP(0, n).derivation();
		auto z = A.evals(x);
		for (int i = 0; i < n; ++i) y[i] /= z[i];
		std::function<std::pair<Poly, Poly>(int, int)> solve = [&](int l, int r) -> std::pair<Poly, Poly> {
			if (r - l == 1) {
				return {Poly(y[l]), Poly({-x[l], 1})};
			}
			int m = (l + r) / 2;
			auto [pl, ql] = solve(l, m);
			auto [pr, qr] = solve(m, r);
			return {pl * qr + pr * ql, ql * qr};
		};
		auto [p, q] = solve(0, x.size());
		return p;
	}
	
	// $a_n = \sum_{i = 1}^{k} f_i a_{n - i}$: https://oi-wiki.org/math/linear-recurrence/
	// find n-th term of The recursive formula for the constant coefficient of order k in $O(k \log k \log n)$
	static valT linearRecursion(std::vector<valT> a, std::vector<valT> f, LL n) {
		if (n < (int)a.size()) return a[n];
		int m = f.size();
		std::reverse(f.begin(), f.end());
		std::vector<valT> g(m);
		g.emplace_back(1);
		Poly A = Poly({0, 1}), p = Poly(g) - Poly(f);
		Poly R = A.powModPoly(n, p);
		return R.inner(a);
	} // https://www.luogu.com.cn/problem/P4723

	// ans[i] = 1^i + 2^i + ... + (n - 1)^i, 0 < i < k
	static std::vector<valT> prefixPowSum(int n, int k) {
		// Poly Numerator = Poly({0, n}).exp(k + 1).divXn(1);
		// Poly denominator  = Poly({0, 1}).exp(k + 1).divXn(1);
		std::vector<valT> a(k), b(k);
		for (int i = 0; i < k; ++i) a[i] = b[i] = Binom.ifac[i + 1];
		valT cur = 1;
		for (int i = 0; i < k; ++i) a[i] *= (cur *= valT::raw(n));
		auto Numerator = Poly(a), denominator = Poly(b);
		
		auto f = (Numerator * denominator.inv(k)).modXn(k) - Poly(1);
		auto ans = f.a;
		ans.resize(k);
		valT now(1);
		for (int i = 2; i < k; ++i) {
			now *= valT(i);
			ans[i] *= now;
		}
		return ans;
	}
	// $\prod_{i = 0}^{n - 1} (x + i)$ in $O(n \log n)$
	static Poly prod(int n) {
		std::function<Poly(int)> solve = [&](int n) -> Poly {
			if (n == 1) return Poly({0, 1});
			int k = n / 2;
			auto A = solve(k);
			std::vector<valT> tmp(k + 1);
			valT now{1};
			for (int i = 0; i <= k; ++i) {
				tmp[i] = now * Binom.ifac[i];
				now *= k;
			}
			auto B = A;
			for (int i = 0; i < B.size(); ++i) B[i] *= Binom.fac[i];
			B = B.mulT(tmp).modXn(k + 1);
			for (int i = 0; i < B.size(); ++i) B[i] *= Binom.ifac[i];
			A *= B;
			if (2 * k != n) {
				B = A;
				for (auto &x : B.a) x *= valT::raw(n - 1);
				A = A.mulXn(1) + B;
			}
			return A;
		};
		return solve(n);
	}
	static Poly prodS(int n) { // $O(n \log^2 n)$
		std::function<Poly(int l, int r)> solve = [&](int l, int r) -> Poly {
			if (r - l == 1) return Poly({l, 1});
			int m = (l + r) / 2;
			return solve(l, m) * solve(m, r);
		};
		return solve(0, n);
	}
	// Stirling number
	static std::vector<valT> stirling1row(int n) {
		auto B = prod(n).a;
		B.resize(n + 1);
		return B;
	}
	std::vector<valT> stirling1col(int n, int k) {
		if (k > n)  return std::vector<valT>(n + 1);
		auto B = Poly({1, -1}).log(n + 2 - k).divXn(1);
		B = (-B).log(n + 1 - k);
		for (auto &x : B.a) x *= valT::raw(k);
		auto ans = B.exp(n + 1 - k).mulXn(k).a;
		ans.resize(n + 1);
		auto ifacK = Binom.ifac[k];
		for (int i = 0; i <= n; ++i) ans[i] *= Binom.fac[i] * ifacK;
		return ans;
	}
	std::vector<valT> stirling2row(int n) {
		auto tmp = Binom.ifac, a = Binom.ifac;
		for (int i = 1; i <= n; i += 2) tmp[i] = -tmp[i];
		for (int i = 0; i <= n; ++i) a[i] *= pow(valT::raw(i), n);
		auto ans = (Poly(a) * Poly(tmp)).a;
		ans.resize(n + 1);
		return ans;
	}
	std::vector<valT> stirling2col(int n, int k) {
		if (k > n)  return std::vector<valT>(n + 1);
		auto A = Poly(Binom.ifac);
		A = A.divXn(1).modXn(n + 1 - k);
		A = A.log(n + 1 - k);
		for (auto &x : A.a) x *= valT::raw(k);
		A = A.exp(n + 1 - k).mulXn(k);
		auto ans = A.a;
		ans.resize(n + 1);
		auto ifacK = Binom.ifac[k];
		for (int i = 0; i <= n; ++i) ans[i] *= Binom.fac[i] * ifacK;
		return ans;
	}
}; // https://www.luogu.com.cn/training/3015#information


template<typename valT>
class PolyBase {
protected:
	void standard() {
		while (!a.empty() && !a.back()) a.pop_back();
	}
	void reverse() {
		std::reverse(a.begin(), a.end());
		standard();
	}
public:
	std::vector<valT> a;
	PolyBase() {}
	PolyBase(valT x) { if (x != valT(0)) a = {x};}
	PolyBase(const std::vector<valT> &_a) : a(_a) { standard();}
	int size() const { return a.size();}
	valT &operator[](int id) { return a[id];}
	valT at(int id) const {
		if (id < 0 || id >= (int)a.size()) return 0;
		return a[id];
	}
};







// There will be `PolyNTT`, `PolyFFT`, `PolyFFTDynamic`, `PolyMFT` provided to suit for different module $M$.

// - PolyMFT: $M > \text{INT_MAX}$
// - PolyFFTDynamic: else if $M$ is uncertain.
// - PolyNTT: else if $M$ is fixed NTT-friendly, such as $M = 998244353$
// - PolyFFT: else
// - PolyOrigin for testing
