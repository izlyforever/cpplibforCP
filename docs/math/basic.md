# baisc.hpp

## powMod

``` cpp
int powMod(int x, int n, int M)
```

$\text{powMod}(x, n, M) = x^n \mod p$

__Constraints__

- $0 \leq x < M$
- $n \geq 0$
- $M \geq 1$

__Complexity__

- $O(\log n)$



## floor  and ceil

``` cpp
template<typename T>
T floor(T a, T n)
template<typename T>
T ceil(T a, T n)
```

$\text{floor}(a, n) = \lfloor \frac{a}{n} \rfloor$

$\text{ceil}(a, n) = \lceil \frac{a}{n} \rceil$

__Constraints__

- $n > 0$
- $T$ can be (unsigned) `int`, `LL`

__Complexity__

- $O(1)$



## `__int128` : input and output

``` cpp
__int128 int128::read()
void int128::print()
```

__Constraints__

- Never use it with `cin` and `cout`

__Complexity__

- $O(\log n)$



## gcd 

``` cpp
// Binary GCD: slightly faster than std::gcd
LL gcd(LL a, LL b) 
```

__Complexity__

- $O(\log \text{lcm}(a, b))$

__Reference__

- [cp-algorithm](https://cp-algorithms.com/algebra/euclid-algorithm.html)



## exGcd

```cpp
template<typename T>
std::tuple<T, T, T> exGcd(T a, T b)
// auto [d, x, y] = exGcd(a, b)
```

where $d = \gcd(a, b) = ax + by$

__Constraints__

- $T$ can be (unsigned) `int`, `long long`

__Complexity__

- $O(\log \text{lcm}(a, b))$



## crt2

``` cpp
std::pair<LL, LL> crt2(LL a1, LL m1, LL a2, LL m2)
// auto [a, m] = crt2(a1, m1, a2, m2) 
```

The Chinese Remainder Theorem shows that
$$
x \equiv a \mod m
$$
equivalence to
$$
x \equiv a_i \mod m_i, \quad i = 1, 2
$$
__Constraints__

- $\text{lcm}(m_1, m_2)$ is in `LL`
- $0 \leq a_1 < m_1, 0 \leq a_2 < m_2$

__Complexity__

- $O(\log \text{lcm}(m_1, m_2))$



## crt

``` cpp
std::pair<LL, LL> crt(const std::vector<std::pair<LL, LL>> &A)
// n = A.size(), a[i] = A[i].first, m[i] = A[i].second, 
```

> use `crt2` above $n - 1$ times,  we have

$$
x \equiv a \mod m
$$

equivalence to
$$
x \equiv a_i \mod m_i, \quad i = 1, 2, \cdots, n
$$
__Constraints__

- $\text{lcm}(m_i)$ is in `LL`
- $0 \leq a_i < m_i$

__Complexity__

- $O(\log \text{lcm}(m_i))$

__Dependence__

- crt2



## spf (foundmental important)

``` cpp
std::vector<int> spf(int N)
// auto sp = spf(N)
```

where `sp[x]` is smallest prime factor of $x$ 

__Complexity__

- $O(\log N)$



## Binom

It is a const singleton class with static const int number $N = 65$, and $C[N][N]$.

$c[i][j] = \binom{i}{j}$

## exmaple

```cpp
Binom C;
std::cout << C(4, 2) << '\n'; // the answer is 6
```



## BinomModp

It is a singleton template class, with typename `valT`

__Members and Methods__

- `fac`, `ifac`, `inv`
- $\text{binom}(n, k) \doteq \binom{n}{k} = fac[n] \cdot ifac[k] \cdot ifac[n - k]$ if $0 \leq k \leq n$, 0 else.

__Constraints__

- `valT` should be `MInt`, `ModInt` or `ModLL` defined in [mod.hpp](mod.md)
- $valT::mod() \geq n$ 

__Complexity__

- $O(n)$



## Lagrange

``` cpp
template<typename valT>
valT Lagrange(const std::vector<valT> &f, int m)
```

Calculate $f(m)$ where $f$ is the Lagrange interpolation on $f(0), f(1), \cdots, f(n - 1)$

__Complexity__

- $O(n)$

__Dependence__

- BinomModp



## powSum

``` cpp
valT powSum(int n, int k, const std::vector<int> &sp)
```

$\displaystyle \text{powSum}(n, k) = \sum_{i = 0}^n i^k$,  where`sp[x]` is smallest prime factor of $x$

__Complexity__

- $O(k)$

__Dependence__

- Lagrange
- spf



## Matrix

It is a class template, contains matrix add, multiplication, and pow

__Complexity__

- `+`: $O(N^2)$
- `*`: $O(N^3)$ 
- `pow(A, n)`: $O(N^3 \log n)$



## MEX

It is a class contains static const int $B = 20$ and set $S$, you can Insert/erase element to $S$, and

$\text{MEX.solve}(x) = MEX_{a_i \in S} (a_i \oplus x)$

__Constraints__

- $\forall x \in S, x < 2^B$

__Complexity__

- $O(|S| \log |S|)$



## trans

 transform `vector<int>` to `vector<valT>`

__Complexity__

- $O(n)$



## BerlekampMassey(every useful)

``` cpp
static std::vector<valT> BerlekampMassey(const std::vector<valT> &a)
```

It return shortest recursive relational formula of $a$.

__Complexity__

- $O(n^2)$

__Example__

- `BerlekampMassey({1, 2, 4, 8, 16}) = {2}`
- `BerlekampMassey({1, 1, 2, 3, 5}) = {1, 1}`

