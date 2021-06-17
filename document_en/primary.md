<head>
	<script type="text/x-mathjax-config">
		MathJax.Hub.Config({
		  tex2jax: {
			skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'],
			inlineMath: [['$','$']],
			processEscapes: true
		  }
		});
	</script>
	<script type="text/javascript" async
	  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/latest.js?config=TeX-MML-AM_CHTML">
	</script>
</head>

# primary.hpp

`primary.hpp` is the foundation of `math.hpp`

## mod.hpp

There are three classes: `MInt`, `ModInt`, `ModLL`. Only `MInt` is a class template and `template<typename valT>`(occur in other place) will be one of them.

> You must init them a mod before use them.

### methods

- Elementary arithmetics: `+, -, *, /, +=, -=, *=, /=`
- C-style operator: `++, --, <<, <<=`(There is no `>>` since $((a + M) \text{>>} x) \neq (a \text{>>} x)$ in general)
- `raw` for constant-factor speedup.
- `pow`, `>>, <<` are friend methods. 
- `inv` is not based on `pow`, since $M$ is not assume to be prime number.


### Example

``` C++
#include <bits/stdc++.h>
#include "cpplib/math/mod.hpp"

constexpr int M = 998244353;
using mod = MInt<M>;

int main() {
	//freopen("in", "r", stdin);
	std::cin.tie(nullptr)->sync_with_stdio(false);
	int a, b;
	std::cin >> a >> b;

	mod a1(a), b1(b);
	std::cout << a1 + b1 << '\n';

	ModInt::setMod(M);
	ModInt a2(a), b2(b);
	std::cout << a2 - b2 << '\n';

	ModLL::setMod(M);
	ModInt a3(a), b3(b);
	std::cout << a3 * b3 << '\n';
	std::cout << a3 / b3 << '\n';
	return 0;
}
```

## fft.hpp

> fast Fourier transform

It contains `dft` and `idft` in namespace `FFT`.

> The size of $a$ is a pow of $2$ before you use `dft(a)` or `idft(a)`

## ntt.hpp

> number theory transform

It contains `dft` and `idft` in template classs `NTT`.

> The template $M$ should be NTT-friendly, and size of $a$ is a pow of $2$ less that $10^6$ before you use `dft(a)` or `idft(a)`


## fmt.hpp

> fast Mobious transform


## basic.hpp

- $\text{floor}(a, n) = \lfloor \frac{a}{n} \rfloor$, $\text{ceil}(a, n) = \lceil \frac{a}{n} \rceil$
- $\text{powMod}(a, n, p) = a^n \mod p$ 
- int128: input and output
- gcd, exGcd
- crt2, crt
- Binom, BinomModp
- Lagrange
- $\displaystyle \text{powSum}(n, k) = \sum_{i = 0}^n i^k$
- matrix multiplication(class template): you may use `Matrix<MInt<998244353>>` to defined a matrix
- quickSort (don't use it)
- $\text{MEX.solve}(x) = MEX_{a_i \in S} (a_i \oplus x)$
- spf: smallest prime factor(foundmental important)
- trans: transform `vector<int>` to `vector<valT>`
- BerlekampMassey: find shortest recursive relational formula