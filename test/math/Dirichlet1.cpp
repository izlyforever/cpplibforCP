
// docs/test/math/Dirichlet1.cpp
#include <bits/stdc++.h>
#define clog(x) std::clog << (#x) << " is " << (x) << '\n';
using LL = long long;
#include "../../cpplib/math/numberTheory.hpp"
using UL = unsigned long long;
std::mt19937_64 rnd64(std::chrono::steady_clock::now().time_since_epoch().count());

int main() {
	//freopen("in", "r", stdin);
	std::cin.tie(nullptr)->sync_with_stdio(false);
	int n = 1e5;
	std::vector<UL> a(n + 1), e(n + 1, 1), mu(n + 1);
	e[0] = 0;
	mu[1] = 1;
	for (int i = 1; i <= n; ++i) {
		for (int j = i * 2; j <= n; j += i) {
			mu[j] -= mu[i];
		}
	}
	for (int i = 1; i <= n; ++i) a[i] = rnd64();
	auto b = a;

	auto c = DirichletProduct(a, e, n);
	mobiousTran(a, n);
	for (int i = 0; i <= n; ++i) assert(a[i] == c[i]);

	c = DirichletProduct(c, mu, n);
	mobiousTranInv(a, n);
	for (int i = 0; i <= n; ++i) assert(a[i] == c[i]);
	for (int i = 0; i <= n; ++i) assert(a[i] == b[i]);


	c = DirichletRevProduct(a, e, n);
	mobiousRevTran(a, n);
	for (int i = 0; i <= n; ++i) assert(a[i] == c[i]);

	c = DirichletRevProduct(c, mu, n);
	mobiousRevTranInv(a, n);
	for (int i = 0; i <= n; ++i) assert(a[i] == c[i]);
	for (int i = 0; i <= n; ++i) assert(a[i] == b[i]);
	
	return 0;
}