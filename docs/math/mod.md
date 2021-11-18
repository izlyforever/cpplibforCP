There are three classes: `MInt`, `ModInt`, `ModLL`. Only `MInt` is a class template.

**`template<typename valT>`(occur in other place) will be one of them.**

> You must assign them a mod before use them.

## methods

- Elementary arithmetics: `+, -, *, /, +=, -=, *=, /=`
- C-style operator: `++, --` and bit operator `<<, <<=`
- `raw` for constant-factor speedup.
- `pow`, `>>, <<` are friend methods.
- `inv` is not based on `pow`, since $M$ is not assume to be prime number.

> There is no bit operator `>>` since $((a + M) \text{>>} x) \neq (a \text{>>} x)$ in general

## Example

``` cpp
#include <bits/stdc++.h>
#include "cpplib/math/mod.hpp"

constexpr int M = 998244353;
using mod = MInt<M>;

int main() {
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