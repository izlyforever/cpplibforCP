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

# mixed.hpp

- GospersHack: n choose k, brute-force(you should implement it to meet for needs)
- $n$-th Fibonacci number
- floorSum: $\displaystyle \sum_{i = 0}^{n - 1} \lfloor \frac{a \cdot i + b}{m} \rfloor$
- sumNum: $\displaystyle \sum_{\sum c_i x_i = m} \frac{(\sum x_i)!}{\prod (x_i !)}$
- decInc: count min time: every time --n or ++m s.t. $n \mid m$
- FirstInRange: finds min x s.t. $L \leq (A x) \mod M \leq R$ (or -1 if it does not exist)
- Gauss: Gauss-Jordan Elimination $Ax = b$, float version
- GaussModp: Gauss-Jordan Elimination $Ax = b$, mod version
- simplex: Simplex algorithm for linear programming
- Karatsuba: Polynomial multiplication with arbitrary modulus $O(n^{\log_2 3})$ 
- KaratsubaParallel: parallel version of Karatsuba(you may need `-lpthread` to complier)
- quadrangleItvDp: Segment DP optim: $O(n^3)$ to $O(n^2)$
- quadrangleRollDp: roll DP optim $O(n^3)$ to $O(n^2)$


## Gauss

``` C++
#include <bits/stdc++.h>
#define clog(x) std::clog << (#x) << " is " << (x) << '\n';
using LL = long long;
#include "../cpplib/math/mixed.hpp"

int main() {
	//freopen("in", "r", stdin);
	std::cin.tie(nullptr)->sync_with_stdio(false);	

	std::vector<std::vector<double>> A{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
	std::vector<double> b{1.0, 2.0, 3.0};
	for (auto x : Gauss(A, b)) std::cout << x << '\n';
	
	return 0;
}
```