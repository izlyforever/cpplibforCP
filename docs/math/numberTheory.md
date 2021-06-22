# numberTheory.hpp

> please use `g++ -o main main.cpp -std=c++17 -O2` to complier examples below.

## Prime

It is a singleton which has member `isP, p, pi`, and method: `primePi` and `nthPrime`

$$
\psi(x,s) = \sum_{n \leq x} |\gcd(n,m_s) == 1| = \sum_{d \mid m_s} u(d)\lfloor \frac{x}{d} \rfloor
$$

where $m_s = p_1 \cdots p_s$

The key point is: if $s \geq \pi(\sqrt{x})$, then 

$$
\psi(x,s) = \pi(x) - s + 1
$$

> See [cnblog](https://www.cnblogs.com/izlyforever/p/computationOfPiX.html) or [origin paper package](https://chachabai.github.io/computationOfPiX/countPrime.zip) for detail.


### Example

``` C++
#include <bits/stdc++.h>
#define clog(x) std::clog << (#x) << " is " << (x) << '\n';
using LL = long long;
#include "cpplib/math/numberTheory.hpp"

int main() {
	//freopen("in", "r", stdin);
	std::cin.tie(nullptr)->sync_with_stdio(false);
	auto start = std::clock();
	auto &prime = Prime::Instance();
	std::clog << "Init time used: " << (std::clock() - start) / 1000 << "ms" << std::endl;	
	
	LL n = prime.primePi(123456789012LL);
	std::cout << n << '\n';

	LL x = prime.nthPrime(n);
	std::cout << x << '\n';

	// It must the same as n
	std::cout << prime.primePi(x) << '\n';
	
	std::clog << "Total time used: " << (std::clock() - start) / 1000 << "ms" << std::endl;
	return 0;
}
```

## Euler and Mobius 

Euler's Totient function and Mobius function. The have many in common.

The key points are: 

$$
\sum_{i = 1}^n \text{sumPhi}(\lfloor \frac{n}{i} \rfloor) = \frac{n(n + 1)}{2}, \quad \sum_{i = 1}^n \text{sumMu}(\lfloor \frac{n}{i} \rfloor) = 1
$$

and then `numberTheoryBlock` tech is used to give a $O(n^{\frac{2}{3}})$ implement.

### Example

``` C++
#include <bits/stdc++.h>
#define clog(x) std::clog << (#x) << " is " << (x) << '\n';
using LL = long long;
#include "cpplib/math/numberTheory.hpp"

int main() {
	//freopen("in", "r", stdin);
	std::cin.tie(nullptr)->sync_with_stdio(false);
	auto start = std::clock();
	auto &euler = Euler::Instance();
	auto &Mobius = Mobius::Instance();
	std::clog << "Init time used: " << (std::clock() - start) / 1000 << "ms" << std::endl;	
	
	int n = 1e9 + 7;
	clog(euler.getPhi(n));
	clog(euler.getSumPhi(n));

	clog(Mobius.getMu(n));
	clog(Mobius.getSumMu(n));
	clog(Mobius.getAbsSum(n));

	std::clog << "Total time used: " << (std::clock() - start) / 1000 << "ms" << std::endl;
	return 0;
}
```

## others

- npf: init numbers of (multi) prime factors less than N in $O(N)$
- factor: list of different prime factors of n
- Factor: list of prime factors of n
- primitiveRoot: smallest primitive root or 0
- primitiveRootAll: list of all primitive roots or empty
- PollardRho: Probabilistic Method: Miller-Rabin prime test and PollardRho big number Decomposition
- babyStepGiantStep: find smallest non-negetive $x$ s.t. $a^x = b \mod p$, or $-1$
- sqrtModp: find $x$ s.t. $x^2 = a \mod p$, or $-1$ in $O(\log^2 p)$
- lcmPair: return all pair $(i, j, lcm(i, j)$  with lcm(i, j) <= n, $O(n \log^2 n)$
- DirichletProduct: Dirichlet Product and fast Mobius transform.


### DirichletProduct Test

``` C++
// docs/test/math/DirichletTest1.cpp
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
```