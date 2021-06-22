#pragma once
#include <bits/stdc++.h>
#include "basic.hpp"
using LL = long long;

// note that p[1] = 2(p[0] is meanless)
class Prime {
	static inline constexpr int N = 5e6 + 2;
	bool isp[N]{};
	std::vector<int> p{0, 2}, pi;
	// $O(N \log \log N)$ but faster when N < 1e9
	std::vector<int> initPrime() {
		isp[2] = true;
		for (int i = 3; i < N; i += 2) isp[i] = true;
		int sq = int(std::sqrt(N + 0.1)) | 1;
		for (int i = 3; i <= sq; i += 2) if (isp[i]) {
			p.emplace_back(i);
			for (int j = i * i; j < N; j += i << 1) isp[j] = false;
		}
		for (int i = sq + 2; i < N; i += 2) if (isp[i]) p.emplace_back(i);
		return p;
	}
	// $O(N)$ but slower when N < 1e9
	std::vector<int> initPrimeS() {
		isp[2] = true;
		for (int i = 3; i < N; i += 2) isp[i] = true;
		for (int i = 3; i < N; i += 2) {
			if (isp[i]) p.emplace_back(i);
			// use t to avoid overflow
			for (int j = 2, t = (N - 1) / i + 1, np = p.size(); j < np && p[j] < t; ++j) {
				isp[i * p[j]] = false; // isp[x] 
				// It will be O(nloglogn) if we remove following code
				if (i % p[j] == 0) break; // % is time-consuming
			}
		}
		return p;
	}
	static inline constexpr int M = 6;
	std::vector<int> phi[M + 1];
	void initPi() {
		pi[2] = 1;
		for (int i = 3; i < N; ++i) {
			pi[i] = pi[i - 1];
			if (isp[i]) ++pi[i];
		}
		std::vector<int> sz(M + 1, 1);
		for (int i = 1; i <= M; ++i) sz[i] = p[i] * sz[i - 1];
		phi[0] = {1}; // phi[0] is meanless
		// optim since phi[j][i] = phi[j][i - 1] - phi[j / p[i]][i - 1]
		for (int i = 1; i <= M; ++i) {
			phi[i].resize(sz[i]);
			for (int j = 0; j < p[i]; ++j) {
				for (int k = 0, jsz = j * sz[i - 1]; k < sz[i - 1]; ++k) {
					phi[i][jsz + k] = j * phi[i - 1].back() + phi[i - 1][k];
				}
			}
			for (int k = 0; k < sz[i - 1]; ++k) {
				for (int j = 0, kp = k * p[i]; j < p[i]; ++j) {
					phi[i][kp + j] -= phi[i - 1][k];
				}
			}
		}
	}
	// See if x \in (s - n, s] is prime assume that p.back() * p.back() >= s
	std::vector<int> seive(LL s, int n) { // O(N log s)
		std::vector<int> isP(n, 1); // isP[i] = 1 means s - i is prime
		for (int i = 1; 1LL * p[i] * p[i] <= s; ++i) {
			for (int j = s % p[i]; j < n; j += p[i]) isP[j] = 0;
		}
		return isP;
	}
	Prime() : pi(N) {
		initPrime();
		initPi();
	}
	bool isPrimeS(LL n) {
		if (n == 2) return true;
		if (n == 1 || n % 2 == 0) return false;
		for (LL i = 3; i * i <= n; i += 2) if (n % i == 0) return false;
		return true;
	}
public:
	int operator[](int i) { return p[i];}
	LL primephi(LL x, int s) {
		if (s <= M) return (x / phi[s].size()) * phi[s].back() + phi[s][x % phi[s].size()];
		if (x / p[s] <= p[s]) return primePi(x) - s + 1;
		if (x / p[s] / p[s] <= p[s] && x < N) {
			int s2x = pi[(int)(std::sqrt(x + 0.2))];
			LL ans = pi[x] - (s2x + s - 2) * (s2x - s + 1) / 2;
			for (int i = s + 1; i <= s2x; ++i) {
				ans += pi[x / p[i]];
			}
			return ans;
		}
		return primephi(x, s - 1) - primephi(x / p[s], s - 1);
	}
	LL primePi(LL x) {
		if (x < N) return pi[x];
		int ps2x = primePi(int(std::sqrt(x + 0.2)));
		int ps3x = primePi(int(std::cbrt(x + 0.2)));
		LL ans = primephi(x, ps3x) + 1LL * (ps2x + ps3x - 2) * (ps2x - ps3x + 1) / 2;
		for (int i = ps3x + 1, ed = ps2x; i <= ed; ++i) {
			ans -= primePi(x / p[i]);
		}
		return ans;
	}
	bool isPrime(LL n) {
		if (n < N) return isp[n];
		if (1LL * p.back() * p.back() > n) return isPrimeS(n);
		for (int i = 1; p[i] * p[i] <= n; ++i) if (n % p[i] == 0) return false;
		return true;
	}
	// DynamicProgramming version O(\frac{n}{\log n}) with n < 10^12
	LL primePiS(LL n) {
		int rn = std::sqrt(n + 0.2);
		std::vector<LL> R(rn + 1);
		for (int i = 1; i <= rn; ++i) R[i] = n / i - 1;
		int ln = n / (rn + 1) + 1;
		std::vector<LL> L(ln + 1);
		for (int i = 1; i <= ln; ++i) L[i] = i - 1;
		for (int p = 2; p <= rn; ++p) {
			if (L[p] == L[p - 1]) continue;
			for (int i = 1, tn = std::min(n / p / p, LL(rn)); i <= tn; ++i) {
				R[i] -= (i <= rn / p ? R[i * p] : L[n / i / p]) - L[p - 1];
			}
			for (int i = ln; i / p >= p; --i) {
				L[i] -= L[i / p] - L[p - 1];
			}
		}
		return R[1];
	}
	LL nthPrime(LL n) { // Newton method
		if (n < (int)p.size()) return p[n];
		LL ans = n * log(n), err = log(n) / log(10);
		LL m = primePi(ans);
		while (m < n || m > n + err) {
			ans += (n - m) / (log(m) - 1) * log(m) * log(m);
			m = primePi(ans);
		}
		int sn = std::sqrt(N);
		while (1) {
			auto isP = seive(ans, sn);
			for (int i = 0; i < sn; ++i) if (isP[i]) {
				if (m-- == n) return ans - i;
			}
			ans -= sn;
		}
	}
	Prime(const Prime&) = delete;
	static Prime& Instance() {
		static Prime instance;
		return instance;
	}
};

class Euler {
	static inline constexpr int N = 5e6 + 2;
	std::vector<int> phi, p{0, 2};
	std::unordered_map<int, LL> mpPhi; 
	std::vector<LL> sumPhi;
	void initPhi() { // $O(N)$
		for (int i = 1; i < N; i += 2) phi[i] = i;
		for (int i = 2; i < N; i += 2) phi[i] = i >> 1;
		for (int i = 3; i < N; i += 2) {
			if (phi[i] == i) p.emplace_back(i), --phi[i];
			for (int j = 2, t = (N - 1) / i + 1, np = p.size(); j < np && p[j] < t; ++j) {
				if (i % p[j] == 0) {
					phi[i * p[j]] = phi[i] * p[j];
					break;
				}
				phi[i * p[j]] = phi[i] * (p[j] - 1);
			}
		}
		for (int i = 2; i < N; i += 4) phi[i] = phi[i >> 1];
		for (int i = 4; i < N; i += 4) phi[i] = phi[i >> 1] << 1;
	}
	LL getPhiS(LL n) {
		if (n % 2 == 0) n /= 2;
		LL r = n;
		while (n % 2 == 0) n /= 2;
		for (LL i = 3; i * i <= n; i += 2) if (n % i  == 0) {
			r = r / i * (i - 1);
			while (n % i == 0) n /= i;
		}
		if (n > 1) r = r / n * (n - 1);
		return r;
	}
	Euler() : phi(N), sumPhi(N) {
		initPhi();
		for (int i = 1; i < N; ++i) sumPhi[i] = sumPhi[i - 1] + phi[i];
	}
public:
	int operator[](int i) { return phi[i];}
	LL getPhi(LL n) {
		if (n < phi.size()) return phi[n];
		if (1LL * p.back() * p.back() > n) return getPhiS(n);
		LL r = n;
		for (int i = 1; 1LL * p[i] * p[i] <= n; ++i) if (n % p[i] == 0) {
			r = r / p[i] * (p[i] - 1);
			while (n % p[i] == 0) n /= p[i];
		}
		if (n > 1) r = r / n * (n - 1);
		return r;
	}
	// min_25 $O(n^{\frac{2}{3}})$
	LL getSumPhi(int n) {
		if (n < N) return sumPhi[n];
		if (mpPhi.count(n)) return mpPhi[n];
		LL r = 1LL * (n + 1) * n / 2;
		for (int i = 2, j; i <= n; i = j + 1) {
			j = n / (n / i);
			r -= (j - i + 1) * getSumPhi(n / i);
		}
		return mpPhi[n] = r;
	}
	Euler(const Euler&) = delete;
	static Euler& Instance() {
		static Euler instance;
		return instance;
	}
};
// https://www.luogu.com.cn/problem/P4213

class Mobious {
	static inline constexpr int N = 5e6 + 2;
	std::vector<int> mu, sumMu, p{0, 2};
	std::unordered_map<int, int> mpMu;
	int getMuS(LL n){
		if (n % 4 == 0) return 0;
		int r = (n % 2 ? 1 : -1);
		if (n % 2 == 0) n /= 2;
		for (LL i = 3; i * i <= n; i += 2) if (n % i == 0) {
			n /= i;
			if (n % i == 0) return 0;
			r = -r;
		}
		return n > 1 ? -r : r;
	}
	void initMuS() { // $O(n log n)$
		mu[1] = 1;
		for (int i = 1; i < N; ++i) {
			for (int j = i * 2; j < N; j += i) {
				mu[j] -= mu[i];
			}
		}
	}
	void initMu() {
		for (int i = 1; i < N; i += 2) mu[i] = i;
		for (int i = 3; i < N; i += 2) {
			if (mu[i] == i) mu[i] = -1, p.emplace_back(i);
			for (int j = 2, t = (N - 1) / i + 1, np = p.size(); j < np && p[j] < t; ++j) {
				if (i % p[j] == 0) {
					mu[i * p[j]] = 0;
					break;
				}
				mu[i * p[j]] = -mu[i];
			}
		}
		for (int i = 2; i < N; i += 4) mu[i] = -mu[i >> 1];
	}
	Mobious() : mu(N), sumMu(N) {
		initMu();
		for (int i = 1; i < N; ++i) sumMu[i] = sumMu[i - 1] + mu[i];
	}
public:
	int operator[](int i) { return mu[i];}
	int getMu(LL n) {
		if (n < mu.size()) return mu[n];
		if (1LL * p.back() * p.back() > n) return getMuS(n);
		int r = 1;
		for (int i = 1; 1LL * p[i] * p[i] <= n; ++i) if (n % p[i] == 0) {
			n /= p[i];
			if (n % p[i] == 0) return 0;
			r = -r;
		}
		return n > 1 ? -r : r;
	}
	LL getAbsSum(LL n) { // Q(n) = Q(n-1) + |mu(n)|
		LL r = 0;
		for (LL i = 1; i * i < n; ++i) {
			r += mu[i] * (n / i / i);
		}
		return r;
	}
	// min_25 $O(n^{\frac{2}{3}})$
	int getSumMu(int n) {
		if (n < N) return sumMu[n];
		if (mpMu.count(n)) return mpMu[n];
		int r = 1;
		for (int i = 2, j; i <= n; i = j + 1) {
			j = n / (n / i);
			r -= 1LL * (j - i + 1) * getSumMu(n / i);	
		}
		return mpMu[n] = r;
	}
	Mobious(const Mobious&) = delete;
	static Mobious& Instance() {
		static Mobious instance;
		return instance;
	}
};
// https://www.luogu.com.cn/problem/P4213

// init numbers of (multi) prime factors less than N in $O(N)$
std::pair<std::vector<int>, std::vector<int>> npf(int N) {
	std::vector<int> np(N, 1), nps(N, 1), p{0, 2};
	nps[0] = nps[1] = 0;
	np[0] = np[1] = 0;
	for (int i = 3; i < N; i += 2) {
		if (nps[i] == 1) p.emplace_back(i);
		for (int j = 2, t, pSize = p.size(); j < pSize && (t = i * p[j]) < N; ++j) {
			nps[t] = nps[i] + 1;
			np[t] = np[i];
			if (i % p[j] == 0) break;
			++np[t];
		}
	}
	for (int i = 2; i < N; i += 4) np[i] = np[i >> 1] + 1;
	for (int i = 4; i < N; i += 4) np[i] = np[i >> 1];
	for (int i = 2; i < N; i += 2) nps[i] = nps[i >> 1] + 1;
	return {np, nps};
}

// list of different prime factors of n
std::vector<int> factor(int n, const std::vector<int> &sp) {
	std::vector<int> ans;
	while (n > 1) {
		int pn = sp[n];
		ans.emplace_back(pn);
		while (n % pn == 0) n /= pn;
	}
	return ans;
}

// list of prime factors of n
std::vector<std::pair<int, int>> Factor(int n, const std::vector<int> &sp) {
	std::vector<std::pair<int, int>> ans;
	while (n > 1) {
		int pn = sp[n], cnt = 0;
		while (n % pn == 0) n /= pn, ++cnt;
		ans.emplace_back(pn, cnt);
	}
	return ans;
}

// smallest primitive root or 0
int primitiveRoot(int n, const std::vector<int> &sp) {
	if (n < 2) return 0;
	if (n == 2 || n == 4) return n - 1;
	if (n % 4 == 0) return 0;
	int n2 = n % 2 == 0 ? n / 2 : n;
	int pn = sp[n2];
	while (n2 % pn == 0) n2 /= pn;
	if (n2 != 1) return 0;
	auto fp = factor(pn - 1, sp);
	auto check = [&](int i) {
		for (auto x : fp) {
			if (powMod(i, (pn - 1) / x, pn) == 1) return false;
		}
		return true;
	};
	int ans = 2;
	while (!check(ans)) ++ans;
	n2 = n % 2 == 0 ? n / 2 : n;
	if (n2 != pn) {
		int m = n2 / pn * (pn - 1);
		auto fm = factor(m, sp);
		for (auto x : fp) {
			if (powMod(ans, m / x, m) == 1) {
				ans += pn;
				break;
			}
		}
	}
	if (n2 != n && (ans % 2 == 0)) ans += n2;
	return ans;
}

// list of all primitive roots or empty
std::vector<int> primitiveRootAllS(int n, const std::vector<int> &sp) {
	int g = primitiveRoot(n, sp);
	if (g == 0) return {};
	if (n == 2 || n == 4) return {n - 1};
	int n2 = n & 1 ? n : n / 2;
	int pn = sp[n2], m = n2 / pn * (pn - 1), now = g;
	std::vector<int> ans{g};
	for (int i = 2; i < m; ++i) {
		now = 1LL * now * g % n;
		if (std::gcd(i, m) == 1) ans.emplace_back(now);
	}
	std::sort(ans.begin(), ans.end());
	return ans;
}

// list of all primitive roots or empty
std::vector<int> primitiveRootAll(int n, const std::vector<int> &sp) {
	if (n < 2) return {};
	if (n == 2 || n == 4) return {n - 1};
	if (n % 4 == 0) return {};
	int n2 = n % 2 == 0 ? n / 2 : n, pn = sp[n2];
	while (n2 % pn == 0) n2 /= pn;
	if (n2 != 1) return {};
	int m = (n & 1 ? n : n / 2) / pn * (pn - 1);
	std::vector<int> vis(n, -1), ans;
	for (int i = 2; i < n; ++i) if (vis[i] == -1 && std::gcd(i, n) == 1) {
		bool flag = true;
		int now = 1;
		for (int j = 1; j < m; ++j) {
			now = 1LL * now * i % n;
			if (now == 1) {
				flag = false;
				break;
			}
			if (std::gcd(j, m) == 1) vis[now] = i;
			else vis[now] = 0;
		}
		if (flag) {
			for (int j = 0; j < n; ++j) if (vis[j] == i) {
				ans.emplace_back(j);
			}
			return ans;
		}
	}
	return {};
} 
// https://www.luogu.com.cn/problem/P6091

// Probabilistic Method: Miller-Rabin prime test and PollardRho big number Decomposition
namespace PollardRho {
std::mt19937_64 rnd64(std::chrono::steady_clock::now().time_since_epoch().count());
LL powModll(LL x, LL n, LL p) {
	LL r = 1;
	while (n) {
		if (n&1) r = __int128(r) * x % p;
		n >>= 1; x = __int128(x) * x % p;
	}
	return r;
}

// m - 1 = m * 2 ^ t, return false if test is invaild
bool witness(LL a, LL n, LL m, int t) {
	LL x = powModll(a, m, n);
	if (x == 1 || x == n - 1) return false;
	while (t--) {
		x = __int128(x) * x %  n;
		if (x == n - 1) return false;
	}
	return true;
}
constexpr int TIMES = 52;
bool rabin(LL n) {
	if (n < 2) return false;
	if (n == 2) return true;
	if (n % 2 == 0) return false;
	LL m = n - 1;
	int t = __builtin_ctzll(m);
	m >>= t;
	for (int cnt = 0; cnt < TIMES; ++cnt) {
		LL a = rnd64() % (n - 1) + 1;
		if (witness(a, n, m, t)) return false;
	}
	return true;
}
LL pollardrho(LL n) {
	LL x = 0, y = 0, z = 1, i = 1, k = 2, c = rnd64() % (n - 1) + 1;
	while (true) {
		x = (__int128(x) * x + c) % n;
		z = __int128(y - x + n) * z % n;
		// optim times of compute gcd
		if (++i == k) {
			LL d = std::gcd(z, n);
			if (d > 1) return d;
			y = x;
			if (k > n) return n;
			k <<= 1;
		}
	}
}
LL spf(LL n) {
	if (rabin(n) || n == 1) return n;
	LL d = n;
	while (d == n) d = pollardrho(n);
	return std::min(spf(d), spf(n / d));
}
LL gpf(LL n, LL mxf = 1) {
	if (rabin(n)) return n;
	if (n <= mxf) return 1;
	LL d = n;
	while (d == n) d = pollardrho(n);
	LL res = gpf(d, mxf);
	return std::max(res, gpf(n / d, std::max(res, mxf)));
}
} // namespace PollardRho


// find smallest non-negetive $x$ s.t. $a^x = b \mod p$, or $-1$(assume $0^0 = 1$)
int babyStepGiantStep(int a, int b, int p) {
	a %= p; b %= p;
	if (p == 1 || b == 1) return 0;
	int cnt = 0, t = 1;
	for (int g = std::gcd(a, p); g != 1; g = std::gcd(a, p)) {
		if (b % g) return -1;
		p /= g; ++cnt;
		b /= g; t = 1LL * t * (a / g) % p;
		if (b == t) return cnt;
	}
	std::map<int, int> mp;
	int m = ceil(std::sqrt(p));
	int base = b;
	for (int i = 0; i != m; ++i) {
		mp[base] = i;
		base = 1LL * base * a % p;
	}
	base = powMod(a, m, p);
	for (int i = 1; i <= m; ++i) {
		t = 1LL * t * base % p;
		if (mp.count(t)) return (1LL * i * m - mp[t]) % p + cnt;
	}
	return -1;
}
// https://www.luogu.com.cn/problem/P4195

// find $x$ s.t. $x^2 = a \mod p$, or $-1$ in $O(\log^2 p)$ Tonelli-Shanks
int sqrtModpS(int a, int p) { // 0 <= a < p < INT_MAX
	if (a == 0 || p == 2) return a;
	auto power = [p](int x, int n) {
		int r = 1;
		while (n) {
			if (n&1) r = 1LL * r * x % p;
			n >>= 1; x = 1LL * x * x % p;
		}
		return r;
	};
	int q = (p - 1) / 2;
	if (power(a, q) != 1) return -1;
	if (q & 1) return power(a, (q + 1) / 2);
	int b; // find a non-quadratic residue
	std::mt19937 rnd(std::chrono::steady_clock::now().time_since_epoch().count());
	while (power(b = rnd() % (p - 1) + 1, q) == 1);
	int c = __builtin_ctz(q);
	q >>= c; // p - 1 = q << (c + 1)
	b = power(b, q);
	int x = power(a, (q + 1) / 2), t = power(a, q);
	// Keep x^2 = a t, t^{2^c} = 1, b^{2^c} = -1
	while (t != 1) {
		// return smallest r s.t. u^{2^r} = -1
		int cc = [p](int u) {
			int r = 0;
			while ((u = 1LL * u * u % p) != 1) ++r;
			return r;
		}(t);
		int d = power(b, 1 << (c - cc - 1)); // d^{2^{cc + 1}} = -1
		// update reason: (xd)^2 = a t d^2, (t d^2)^{2^{cc}} = 1, (d^2)^{2^cc} = -1
		x = 1LL * x * d % p;
		b = 1LL * d * d % p;
		t = 1LL * t * b % p;
		c = cc;
	}
	return x;
}

struct pseudoComplex {
	int x, y;
	static inline int p, m;
	static void setMod(int _p, int _m) { p = _p, m = _m;}
	pseudoComplex(int _x = 0, int _y = 0) : x(_x), y(_y) {};
	pseudoComplex operator*(const pseudoComplex& A) const {
		return pseudoComplex((1LL * x * A.x + 1LL * y * A.y % p * m) % p, (1LL * x * A.y + 1LL * y * A.x) % p);
	}
};
// find $x$ s.t. $x^2 = a \mod p$, or $-1$ in $O(\log p)$ Cipolla
int sqrtModp(int a, int p) {
	if (a == 0 || p == 2) return a;
	auto power = [p](int x, int n) {
		int r = 1;
		while (n) {
			if (n&1) r = 1LL * r * x % p;
			n >>= 1; x = 1LL * x * x % p;
		}
		return r;
	};
	int q = (p - 1) / 2;
	if (power(a, q) != 1) return -1;
	if (q & 1) return power(a, (q + 1) / 2);
	std::mt19937 rnd(std::chrono::steady_clock::now().time_since_epoch().count());
	int b, m; // find a non-quadratic residue
	do {
		b = rnd() % p;
		m = (1LL * b * b - a) % p;
		if (m < 0) m += p;
	} while (power(m, q) == 1);
	int n = (p + 1) / 2;
	pseudoComplex::setMod(p, m);
	pseudoComplex R(1, 0), A(b, 1);
	while (n) {
		if (n & 1) R = R * A;
		n >>= 1;   A = A * A;
	}
	return R.x;
}
// https://www.luogu.com.cn/problem/P5491

// return all pair $(i, j, lcm(i, j)$  with lcm(i, j) <= n, $O(n \log^2 n)$
std::vector<std::tuple<int, int, int>> lcmPair(int n) {
	std::vector<std::pair<int, int>> ed;
	for (int i = 1; i <= n; ++i) {
		for (int j = 1; j < i && i * j <= n; ++j) if (std::__gcd(i, j) == 1) {
			ed.emplace_back(j, i);
		}
	}
	std::vector<int> deg(n + 1);
	std::vector<std::tuple<int, int, int>> edge;
	auto addedge = [&](int i, int j, int d) {
		++deg[i];
		++deg[j];
		edge.emplace_back(i, j, d);
	};
	for (auto [u, v] : ed) {
		int uv = u * v;
		for (int i = u, j = v, d = uv; d <= n; i += u, j += v, d += uv) {
			addedge(i, j, d);
		}
	}
	return edge;
}


template<typename T>
class DirichletProduct {
	static inline auto &prime = Prime::Instance();
	static inline int n = 1e5 + 2; // assume they have the same lenth
public:
	std::vector<T> a;
	static void setLen(int _n) { n = _n;}
	DirichletProduct() : a(n + 1) {}
	DirichletProduct(const std::vector<T> &_a) : a(_a) {}
	// $O(n \log n)$: $c[k] = \sum_{i j = k} a[i] b[j]$
	DirichletProduct operator*(const DirichletProduct& A) const {
		std::vector<T> c(n + 1);
		for (int i = 1; i <= n; ++i) {
			for (int j = 1; i * j <= n; ++j) {
				c[i * j] += a[i] * A.a[j];
			}
		}
		return DirichletProduct(c);
	}
	DirichletProduct& operator*=(const DirichletProduct& A) {
		return *this = (*this) * A;
	}
	// $O(n \log \log n)$
	void mobious() {   // new_a[n] = \sum_{d | n} old_a[d]
		for (int i = 1; prime[i] <= n; ++i) {
			for (int j = 1; j * prime[i] <= n; ++j) {
				a[j * prime[i]] += a[j];
			}
		}
	}
	// $O(n \log \log n)$
	void mobiousInv() { // old_a[n] = \sum_{d | n} new_a[d]
		for (int i = prime.primePi(n); i; --i) {
			for (int j = n / prime[i]; j; --j) {
				a[j * prime[i]] -= a[j];
			}
		}
	}
	// O(n \log n) $c[i] = \sum_{j} a[i j] b[j]$
	DirichletProduct mulT(const DirichletProduct& A) {
		std::vector<T> c(n + 1);
		for (int i = 1; i <= n; ++i) {
			for (int j = 1; i * j <= n; ++j) {
				c[i] += a[i * j] * A.a[j];
			}
		}
		return DirichletProduct(c);
	}
	// $O(n \log \log n)$
	void mobiousT() {   // new_a[d] = \sum_{d | n} old_a[n]
		for (int i = 1; prime[i] <= n; ++i) {
			for (int j = 1; j * prime[i] <= n; ++j) {
				a[j] += a[j * prime[i]];
			}
		}
	}
	// $O(n \log \log n)$
	void mobiousInvT() { // old_a[d] = \sum_{d | n} new_a[n]
		for (int i = 1; prime[i] <= n; ++i) {
			for (int j = n / prime[i]; j; --j) {
				a[j] -= a[j * prime[i]];
			}
		}
	}
	friend std::ostream& operator<<(std::ostream &out, const DirichletProduct &A) {
		for (auto x : A.a) out << x << ' ';
		// out << '\n';
		return out;
	}
};
// reference(many mistakes): https://www.cnblogs.com/ivorysi/p/8889154.html
