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

## Euler and Mobious

Euler's Totient function and Mobious function. The have many in common.

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
	auto &mobious = Mobious::Instance();
	std::clog << "Init time used: " << (std::clock() - start) / 1000 << "ms" << std::endl;	
	
	int n = 1e9 + 7;
	clog(euler.getPhi(n));
	clog(euler.getSumPhi(n));

	clog(mobious.getMu(n));
	clog(mobious.getSumMu(n));
	clog(mobious.getAbsSum(n));

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
- DirichletProduct: Dirichlet Product and fast Mobious transform.


### DirichletProduct Test

``` C++
#include <bits/stdc++.h>
#define clog(x) std::clog << (#x) << " is " << (x) << '\n';
using LL = long long;
#include "../cpplib/math/numberTheory.hpp"

int main() {
	//freopen("in", "r", stdin);
	std::cin.tie(nullptr)->sync_with_stdio(false);
	std::vector<int> a{0, 1, 10, 100, 1000, 10000, 100000};
	std::vector<int> b{0, 1, 1, 1, 1, 1, 1};
	std::vector<int> c{0, 1, -1, -1, 0, -1, 1}; // mu

	DirichletProduct<int>::setLen(int(a.size() - 1));
	DirichletProduct A(a), B(b), C(c);
	auto D = A;
	std::cout << "Fast mobious Test: \n";
	std::cout << A * B << '\n';
	A.mobious();
	std::cout << A << '\n';
	A.mobiousInv();

	std::cout << "\nFast mobious invserse Test: \n";
	std::cout << D * C << '\n';
	D.mobiousInv();
	std::cout << D << '\n';	

	std::cout << "\nFast transpose mobious Test: \n";
	D = A;
	std::cout << A.mulT(B) << '\n';
	A.mobiousT();
	std::cout << A << '\n';
	A.mobiousInvT();

	std::cout << "\nFast transpose mobious invserse Test: \n";
	std::cout << D.mulT(C) << '\n'; // overflow for int occur
	D.mobiousInvT();
	std::cout << D << '\n';	
	
	std::cout << "\nSee origin A: \n";
	std::cout << A;
	
	return 0;
}
```