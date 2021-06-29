# mixed.hpp

## BigInt

See [BigInt](https://github.com/izlyforever/BigInt) for detail.

**Complexity**

- add, sub: $O(n)$
- mul, div: $O(n \log n)$

> where $n$ is the length of number in decimal.

## GospersHack

``` cpp
void GospersHack(int n, int k)
```

brute-force all case: $n$ choose $k$, 1 stand for choosen

> you should implement it to meet to feed your needs

**Complexity**

- $O(\binom{n}{k})$



## Fib

``` cpp
int Fib(int n, int M)
```

return $n$-th Fibonacci number mod $M$.

**Complexity**

- $O(\log n)$



## floorSum

``` cpp
LL floorSum(int n, int m, int a, int b)
```

$\displaystyle \text{floorSum}(n, m, a, b) = \sum_{i = 0}^{n - 1} \lfloor \frac{a \cdot i + b}{m} \rfloor$

**Complexity**

- $O(\log m)$



## sumNum

``` cpp
int sumNum(const std::vector<int> &c, int m, int M)
```

$\displaystyle \text{sumNum}(c, m, M) = \sum_{\sum c_i x_i = m} \frac{(\sum x_i)!}{\prod (x_i !)} \mod M$

**Complexity**

- $O(m)$



## decInc

``` cpp
int decInc(int n, int m)
```

count min time: every time `--n` or `++m` s.t. $n \mid m$

**Complexity**

- $O(\sqrt{m})$

**Example**

- $\text{decInc}(3, 3) = 0$
- $\text{decInc}(5, 3) = 2$
- $\text{decInc}(4, 9) = 1$



## FirstInRange

``` cpp
int FirstInRange(int a, int m, int l, int r)
```

finds min $x$ s.t. $l \leq a x \mod m \leq r$ (or -1 if it does not exist)

**Constraints**

- $0 \leq l \leq r < m$
- $0 \leq a < M$

**Complexity**

- $O(\log m)$



## Gauss

``` cpp
std::vector<double> Gauss(std::vector<std::vector<double>> A, std::vector<double> b)
```

Gauss-Jordan Elimination $Ax = b$, float version, Inspire by [spookywooky](https://codeforces.com/profile/spookywooky)

**Complexity**

- $O(n^3)$



## GaussModp

``` cpp
std::vector<valT> GaussModp(std::vector<std::vector<valT>> A, std::vector<valT> b)
```

Gauss-Jordan Elimination $Ax = b$, mod version

**Complexity**

- $O(n^3)$

## simplex

``` cpp
using VD = std::vector<double>
VD simplex(VD c, std::vector<VD> Aq, VD bq, std::vector<VD> Alq, VD blq)
```

Simplex algorithm for linear programming: compute max $cx$

with constraints: $Aq \cdot x = bq$ and $Alq \cdot x \leq blq$



## Karatsuba

``` cpp
using VL = std::vector<LL>
VL Karatsuba(VL a, VL b, LL p)
```

Polynomial multiplication with arbitrary modulus 

**Complexity**

- $O(n^{\log_2 3})$

> There is a parallel version of Karatsuba(you may need `-lpthread` to complier): `KaratsubaParallel`



## quadrangleItvDp

``` cpp
std::vector<std::vector<T>> quadrangleItvDp(std::vector<std::vector<T>> w, int n)
```

$f_{l, r} = \min_{l \leq k < r} f_{l, k} + f_{k + 1, r} + w(l, r) \qquad (1 \leq l < r \leq n)$

It is a common tech for Segment DP optim: $O(n^3)$ to $O(n^2)$

**Complexity**

- $O(n^2)$



## quadrangleRollDp

``` cpp
std::vector<std::vector<T>> quadrangleRollDp(std::vector<std::vector<T>> w, int n, int m)
```

$f_{i, j} = \min_{k < j} f_{i - 1, k} + w(k + 1, j) \quad (1 \leq i \leq n, 1 \leq j \leq m)$

It is a common tech for oll DP optim: $O(n^3)$ to $O(n^2)$

**Complexity**

- $O(n m)$



## PalindromeNumber

It is a singleton class. Recall that $n$  is a `Palindrome Number` if it's digits representation is a palindrome, for example `1, 121, 33, 23532` are `Palindrome Number` which `132, 112` are not.

``` cpp
LL nthPalindrome(int k)
```

return $k$-th Palindrome Number.

**Constraints**

- $k < 10^9$

**Complexity**

- $O(\log k)$

**Example**

- $\text{nthPalindrome}(1) = 1$
- $\text{nthPalindrome}(10) = 11$
- $\text{nthPalindrome}(20) = 111$



``` cpp
int Palindrome(LL n)
```

return numbers of Palindrome less that $n$

**Constraints**

- $n < 10^{18}$

**Complexity**

- $O(\log n)$

**Example**

- $\text{Palindrome}(1) = 0$
- $\text{Palindrome}(111) = 19$



``` cpp
LL solve(LL n, int k) { return nthPalindrome(k + Palindrome(n)); }
```

return $k$-th Palindrome Number greater than or equal to $n$

