# polyALL.hpp

> The size: $N$ should be less than $10^6$ or $2^{22} \doteq 4 \cdot 10^6$ at least.

`poly.hpp` support almost every algorithm involved polynomial and __the module number $M$ can be any prime number__.

## how to choose 

There are `PolyNTT`, `PolyFFT`, `PolyFFTDynamic`, `PolyMFT` provided to suit for different module $M$.

- PolyMFT: $M > \text{INT_MAX}$
- PolyFFTDynamic: else if $M$ is uncertain.
- PolyNTT: else if $M$ is fixed NTT-friendly, such as $M = 998244353$ or `PolyS` instead
- PolyFFT: else

> `PolyOrigin` for testing.

## polyS.hpp

`PolyS.hpp` is a simple and small version of Poly.  It contains basic operators: `+, -, *, /`, log, exp, sqrt, and multi-evaluation for fiexed mod = 998244353.

> You may change `NTTS::M = 998244353` to other NTT-friendly prime number(and primitive root `NTTS::g = 3`). 

## Arbitrary module

However $M$ should be bigger than the size $N$ since some function need to assmue $1, \cdots, N - 1$ invertible in $\mod M$

Two ways to support it.

- FFT based: you should check if the precision sufficient
- NTT based: use 3 or 4 or more NTT-friendly modules, are then use __Chinese remainder theorem__

we choose `M0 = 595591169, M1 = 645922817, M2 = 897581057, M3 = 998244353` in PolyMFT using following sageMath code:

``` Python
ans = []
for i in range(23, 28):
    for j in range(1, 1000, 2):
        if(j * 2 ^i + 1 < 2^30 and is_prime(j*2^i+1) and primitive_root(j*2^i+1) == 3):
            ans.append(j*2^i+1)

for i in sorted(ans):
    print(i, "\t = 1 + ", factor(i-1))

# output:
# 167772161        = 1 +  2^25 * 5
# 469762049        = 1 +  2^26 * 7
# 595591169        = 1 +  2^23 * 71
# 645922817        = 1 +  2^23 * 7 * 11
# 897581057        = 1 +  2^23 * 107
# 998244353        = 1 +  2^23 * 7 * 17
```


## Method

- elementary: `+, -, *, /, %, +=, -=, *=, /=, %=, -(negative)`, inv, mulXn, modXn, divXn.
- fundamental: powModPoly, inner, derivation, integral, log, exp, sqrt, mulT,  evals, Lagrange, linearRecursion, prod, stirling number(stirling1row, stirling1col, stirling2row, stirling2col)
- mixed: sin, cos, asin, atan, compose, composeInv, toFallingPowForm, fromFallingPowForm, valToVal
- prefixPowSum: $1^i + 2^i + \cdots + (n - 1)^i,  0 < i < k$
- sumFraction: $\sum_{i = 0}^{n - 1} a_i / (1 - b_i x)$

__As an application, we compute $n!$ $O(\sqrt{n} \log^2 n)$ and $O(\sqrt{n} \log n$(introduced by min_25) in `poly.hpp`__

## Example

``` cpp
#include <bits/stdc++.h>
#define clog(x) std::clog << (#x) << " is " << (x) << '\n';
using LL = long long;
#include "cpplib/math.hpp"

template<typename T>
void debug(std::vector<T> a){
	for (auto &i : a) std::cout << i << ' ';
	std::cout << std::endl; 
}

int main() {
	//freopen("in", "r", stdin);
	std::cin.tie(nullptr)->sync_with_stdio(false);	

	std::vector<int> a{1, 2, 3, 4};
	std::vector<int> b{1, 2, 3};
	PolyS A1(a), B1(b);
	auto c1 = (A1 * B1).a;
	debug(c1);

	using modM = MInt<998244353>;
	PolyNTT A2(trans<modM>(a)), B2(trans<modM>(b));
	auto c2 = (A2 * B2).a;
	debug(c2);

	// you must setMod before using it
	ModInt::setMod(998244353);
	PolyFFTDynamic A3(trans<ModInt>(a)), B3(trans<ModInt>(b));
	auto c3 = (A3 * B3).a;
	debug(c3);

	ModLL::setMod(998244353);
	PolyMFT A4(trans<ModLL>(a)), B4(trans<ModLL>(b));
	auto c4 = (A4 * B4).a;
	debug(c4);

	return 0;
}
```
