# numberTheory.hpp

$p$ is assume to be a prime number, and $M$ for arbitrary module.

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

__Constraints__

- primePi(n): $n < N^2$ 
- nthPrime(n): $n < (\frac{N}{\ln n})^2$

__Complexity__

- primePi(n): $O(n^{\frac{2}{3}})$
- nthPrime(n): $O(n^{\frac{2}{3}} \log^2 n)$


### Example

``` cpp
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

``` cpp
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


## npf

``` cpp
std::pair<std::vector<int>, std::vector<int>> npf(int N)
// auto [a, b] = npf(N)
```

init numbers of (multi) prime factors less than $N$

__Complexity__

- $O(N)$

__Example__

- $a[4] = 1, b[4] = 2$
- $a[21] = 2, b[21] = 2$
- $a[72] = 2, b[72] = $5



## factor and Factor

``` cpp
std::vector<int> factor(int n, const std::vector<int> &sp)
std::vector<std::pair<int, int>> Factor(int n, const std::vector<int> &sp)
```

factor  $n$ into prime number

__Complexity__

- $O(\log n)$

__Example__

- $\text{factor}(60, sp) = \{2, 3, 5\}$
- $\text{Factor}(60, sp) = \{\{2, 2\}, \{3, 1\}, \{5, 1\}\}$



## primitiveRoot

``` cpp
int primitiveRoot(int n, const std::vector<int> &sp)
```

return smallest primitive root or 0 if not exist.

$n$ have primitive root if and only if $n = 2, 4, p^n, 2 p^n$ where $p > 2$ is prime number.

__Complexity__

- $O(\log n)$

### Example

- $\text{primitiveRoot}(998244353) = 3$



## primitiveRootAll

``` cpp
std::vector<int> primitiveRootAll(int n, const std::vector<int> &sp)
```

return list of all primitive roots or empty if not exist

__Complexity__

- $O(n)$

### Example

- $\text{primitiveRoot}(5) = \{2, 3\}$



## PollardRho

Pollard-rho algorithm is a probabilistic method for __big number decomposition__, which is based on __big prime test__ probabilistic method: Miller-Rabin

``` cpp
bool PollardRho::rabin(LL n)
LL PollardRho::spf(LL n)
LL PollardRho::gpf(LL n, LL mxf = 1)
```

__Complexity__

- $O(n^{\frac{1}{4}} \log n)$



## babyStepGiantStep

find smallest non-negative $x$ s.t. $a^x = b \mod p$, or $-1$

__Constraints__

- $p$ is prime 

__Complexity__

- $O(\sqrt{p} \log p)$

### 

## sqrtModp

``` cpp
int sqrtModp(int a, int p) 
```

__Constraints__

- $p$ is prime 

__Complexity__

- $O(\log p)$



## lcmPair

``` cpp
std::vector<std::tuple<int, int, int>> lcmPair(int n)
```

return all pair $(i, j, \text{lcm}(i, j)$  with $\text{lcm}(i, j) \leq n$ 

__Complexity__

- $O(n \log^2 n)$



## DirichletProduct

give two function, $f, g$, we define Dirichlet Product of $f, g$ as $f \star g$ defined as
$$
(f \star g)(n) \doteq \sum_{d | n} f(d) g(\frac{n}{d})
$$

- $g \equiv 1$, then $f \star g$  is call Mobius transform, 
- $g \equiv \mu$, then  $f \star g$  is call Mobius inverse transform,  where $\mu$ is mobius function

__Complexity__

- DirichletProduct: $O(n \log n)$
- Mobius transform: $O(n \log \log n)$
- Mobius inverse transform: $O(n \log \log n)$


### DirichletProduct Test

``` cpp
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

