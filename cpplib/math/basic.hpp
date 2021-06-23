#pragma once
#include <bits/stdc++.h>
using LL = long long;

int powMod(int x, int n, int M) {
	int r = 1;
	while (n) {
		if (n&1) r = 1LL * r * x % M;
		n >>= 1; x = 1LL * x * x % M;
	}
	return r;
}

template<typename T>
T floor(T a, T n) { // n > 0
	return a < 0 ? (a - n + 1) / n : a / n;
}
template<typename T>
T ceil(T a, T n) { // n > 0
	return a < 0 ? a / n : (a + n - 1) / n;
}

// never mixed it with cin and cout
namespace int128 {
__int128 read(){
	__int128 x = 0;
	bool negative = false;
	char ch = getchar();
	while (ch < '0' || ch > '9'){
		if (ch == '-') negative = true;
		ch = getchar();
	}
	while (ch >= '0' && ch <= '9') {
		x = x * 10 + ch - '0';
		ch = getchar();
	}
	return negative ?  -x : x;
}

void printS(__int128 x){
	if (x > 9) printS(x / 10);
	putchar(x % 10 + '0');
}
void print(__int128 x){ 
	if (x < 0) {
		putchar('-');
		x = -x;
	}
	printS(x);
}
} // namespace int128


// slightly faster than std::gcd
LL gcd(LL a, LL b) {
	if (!a || !b) return a | b;
	unsigned shift = __builtin_ctzll(a | b);
	a >>= __builtin_ctzll(a);
	do {
		b >>= __builtin_ctzll(b);
		if (a > b) std::swap(a, b);
		b -= a;
	} while (b);
	return a << shift;
}
// https://cp-algorithms.com/algebra/euclid-algorithm.html

// ax + by = gcd(a,b)
template<typename T>
std::tuple<T, T, T> exGcd(T a, T b) {
	if (b == 0) return {a, 1, 0};
	auto [d, y, x] = exGcd(b, a % b);
	return {d, x, y - a / b * x};
}

// Chinese remainder theorem: x = ai mod mi, m_i > 0, 0 <= a_i < m_i
std::pair<LL, LL> crt2(LL a1, LL m1, LL a2, LL m2) {
	auto [d, t1, t2] = exGcd(m1, m2);
	assert((a2 - a1) % d == 0);
	LL m = m1 / d * m2;
	LL ans = (a1 + (a2 - a1) / d * t1 % m2 * m1) % m;
	return std::make_pair(ans < 0 ? ans + m: ans, m);
}

std::pair<LL, LL> crt(const std::vector<std::pair<LL, LL>> &A) {
	auto ans = A[0];
	for (int i = 1, na = A.size(); i < na; ++i) {
		ans = crt2(ans.first, ans.second, A[i].first, A[i].second);
	}
	return ans;
}
// https://www.luogu.com.cn/problem/P1495

// $O(N)$ smallest prime factor
std::vector<int> spf(int N) {
	std::vector<int> sp(N), p{0, 2};
	for (int i = 2; i < N; i += 2) sp[i] = 2;
	for (int i = 1; i < N; i += 2) sp[i] = i;
	for (int i = 3; i < N; i += 2) {
		if (sp[i] == i) p.emplace_back(i);
		for (int j = 2, np = p.size(); j < np && p[j] <= sp[i] && i * p[j] < N; ++j) {
			sp[i * p[j]] = p[j]; // Note that sp[x] is assigned only once foreach x
		}
	}
	return sp;
}

class Binom {
	static inline const int N = 65;
	LL C[N][N];
	Binom() {
		for (int i = 0; i < N; ++i) C[i][0] = C[i][i] = 1;
		for (int i = 1; i < N; ++i) {
			for (int j = 1; j < i; ++j) {
				C[i][j] = C[i - 1][j] + C[i - 1][j - 1];
			}
		}
	}
public: 
	Binom(const Binom&) = delete;
	// Binom& operator=(const Binom&) = delete; // can be erase
	static Binom& Instance() {
		static Binom instance;
		return instance;
	}
	LL operator()(const int &m, const int &n) {
		assert(n < N && m < N);
		return C[m][n];
	}
};

template<typename valT>
class BinomModp {
	int sz;
	BinomModp() : sz(2), fac({1, 1}), ifac({1, 1}), inv({0, 1}) {}
	void init(int N) {
		if (sz > N) return;
		const int M = valT::mod();
		assert(N <= M);
		fac.resize(N), ifac.resize(N), inv.resize(N);
		for (int i = sz; i < N; ++i) inv[i] = inv[M % i] * valT::raw(M - M / i);
		for (int i = sz; i < N; ++i) fac[i] = fac[i - 1] * valT::raw(i);
		for (int i = sz; i < N; ++i) ifac[i] = ifac[i - 1] * inv[i];
		// another way to compute ifac
		// ifac[n - 1] = fac[n - 1].inv();
		// for (int i = n - 1; i > 0; --i) ifac[i - 1] = ifac[i] * T(i);
		sz = N;
	}
public:
	std::vector<valT> fac, ifac, inv;
	BinomModp(const BinomModp&) = delete;
	static BinomModp& Instance(int N) {
		static BinomModp instance;
		instance.init(N);
		return instance;
	}
	valT binom(int n, int k) {
		if (n < 0 || n < k) return valT(0);
		assert(n < sz);
		return fac[n] * ifac[k] * ifac[n - k];
	}
	// M is a small prime number in this case
	valT lucas(int n, int k) {
		valT r(1);
		const int M = valT::mod();
		while (n && k) {
			int np = n % M, kp = k % M;
			if (np < kp) return valT(0);
			r *= binom(np, kp); 
			n /= M; k /= M;
		}
		return r;
	}
};

// Calculate f(m) where f is the Lagrange interpolation on $f(0), f(1), \cdots, f(n - 1)$ 
template<typename valT>
valT Lagrange(const std::vector<valT> &f, int m) {
	int n = f.size();
	if (m < n) return f[m];
	auto &B = BinomModp<valT>::Instance(n);
	std::vector<valT> AP(n), BP(n);
	AP[0] = BP[n - 1] = valT(1);
	for (int i = 1; i < n; ++i) AP[i] = AP[i - 1] * valT::raw(m + 1 - i);
	for (int i = n - 2; ~i; --i) BP[i] = BP[i + 1] * valT::raw(m - 1 - i);
	valT ans = 0;
	for (int i = 0; i < n; ++i) {
		valT x = f[i] * AP[i] * BP[i] * B.ifac[i] * B.ifac[n - 1 - i];
		ans += (n - 1 - i) & 1 ? -x : x;
	}
	return ans;
}
// Lagrange theorem $f(x) =  \sum_{i = 0}^{n - 1} f_i \prod_{j \neq i} \frac{x - j}{i - j}$
// Simplies $f(m) = \sum_{i = 0}^{n - 1} (-1)^{n - 1 - i} f_i \binom{m}{i} \binom{m - i - 1}{n - 1 - i}$

// Calculate powSum in $O(k)$ based on Lagrange interpolation
template<typename valT>
valT powSum(int n, int k, const std::vector<int> &sp) {
	if (k == 0) return valT(n);
	std::vector<valT> f(k + 2);
	f[1] = valT(1);
	for (int i = 2, nf = f.size(); i < nf; ++i) {
		if (sp[i] == i) f[i] = pow(valT(i), k);
		else f[i] = f[sp[i]] * f[i / sp[i]];
	}
	for (int i = 1, nf = f.size(); i < nf; ++i) f[i] += f[i - 1];
	return Lagrange(f, n);
}
// https://codeforces.com/problemset/problem/622/F


template<typename valT>
class Matrix {
	static inline constexpr int N = 1003;
	int n;
public:
	valT a[N][N];
	Matrix() {}
	Matrix(int _n, valT x = 0): n(_n) {
		all(0);
		for (int i = 0; i < n; ++i) {
			a[i][i] = x;
		}
	}
	void all(valT x) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				a[i][j] = x;
			}
		}
	}
	Matrix operator+(const Matrix &A) const {
		Matrix R(n);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				R.a[i][j] += a[i][j];
			}
		}
		return R;
	}
	Matrix operator*(const Matrix &A) const {
		Matrix R(n);
		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < n; ++k) {
				for (int j = 0; j < n; ++j) {
					R.a[i][j] += a[i][k] * A.a[k][j];
				}
			}
		}
		return R;
	}
	void print() {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				std::cout << a[i][j] << ' ';
			}
			std::cout << '\n';
		}
	}
	friend Matrix pow(Matrix A, int n) {
		Matrix R(A.n, valT(1));
		while (n) {
			if (n&1) R = R * A;
			n >>= 1; A = A * A;
		}
		return R;
	}
};

// You must change mid to be random one or you will be hack
template<typename T> // don't use it
void quickSort(std::vector<T> &a) {
	std::function<void(int, int)> qSort = [&](int l, int r) {
		int i = l, j = r;
		auto mid = a[(l + r) / 2];
		while (i <= j) {
			while (i <= j && a[i] < mid) ++i;
			while (j >= i && a[j] > mid) --j;
			if (i <= j) {
				std::swap(a[i], a[j]);
				++i;
				--j;
			}
		}
		if (i < r) qSort(i, r);
		if (l < j) qSort(l, j);
	};
	qSort(0, a.size() - 1);
}

class MEX {
	// B may need to be bigger
	static inline constexpr int B = 20;
	std::array<std::map<int, int>, B> mp;
	std::map<int, int> S;
public:
	void insert(int x) {
		if (S[x]++) return;
		int mask = 0;
		for (int i = B - 1; i >= 0; --i) {
			mask |= 1 << i;
			++mp[i][x & mask];
		}
	}
	void erase(int x) {
		if (--S[x]) return;
		S.erase(x);
		int mask = 0;
		for (int i = B - 1; i >= 0; --i) {
			mask |= 1 << i;
			--mp[i][x & mask];
		}
	}
	// find mex(a_i ^ x)
	int solve(int x = 0) {
		int mask = 0, r = 0;
		for (int i = B - 1; i >= 0; --i) {
			mask |= x & (1 << i);
			if (mp[i][mask] == (1 << i)) {
				mask ^= 1 << i;
				r |= 1 << i;
			}
		}
		return r;
	}
};


// transform vector<int> to vector<valT>
template<typename valT>
std::vector<valT> trans(const std::vector<int> &a) {
	int n = a.size();
	std::vector<valT> ans(n);
	for (int i = 0; i < n; ++i) ans[i] = valT(a[i]);
	return ans;
}

// Shortest recursive relational formula: https://cmwqf.github.io/2020/07/18/%E6%B5%85%E8%B0%88Berlekamp-Massey%E7%AE%97%E6%B3%95/
template<typename valT>
static std::vector<valT> BerlekampMassey(const std::vector<valT> &a) {
	std::vector<valT> ans, lst;
	valT delta = 0;
	for (int i = 0, w = -1, n = a.size(); i < n; ++i) {
		valT t = 0;
		for (int j = 0, na = ans.size(); j < na; ++j) {
			t += ans[j] * a[i - 1 - j];
		}
		if (t == a[i]) continue;
		// first time ans fail
		if (w == -1) {
			w = i; delta = a[i];
			ans.emplace_back(0);
			continue;
		}
		auto now = ans;
		auto mul = (a[i] - t) / delta;
		if (i - w + lst.size() > ans.size()) ans.resize(i - w + lst.size());
		ans[i - w - 1] += mul;
		for (int j = 0, lj = lst.size(); j < lj; ++j) ans[i - w + j] -= mul * lst[j];
		if ((int)now.size() - i < (int)lst.size() - w) {
			w = i; delta = a[i] - t; std::swap(now, lst);
		}
	}
	return ans;
}
// https://www.luogu.com.cn/problem/P5487