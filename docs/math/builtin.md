# builtin.hpp

## ctz

``` cpp
int ctz32(unsigned x);
int ctz64(unsigned long long x);
```

Slow than __builtin_ctz and __builtin_ctzll but funny.

## bitCount

we may use `__builtin_popcount`, or `__builtin_popcountll` in g++(but not in clang++), but it is slower than the following methods

``` cpp
int bitCount(unsigned n);  
int bitCountll(unsigned long long n); // The fastest so far for 64bit
int BitCountTable(unsigned n);    // The fastest so far for 32bit
int BitCountTablell(unsigned long long n);
```

`bitCountTable` use static Table of length 256, and `bitCount` and `bitCountll` use MIT HAKMEM, since $2^6 = 64 > 32$, and $2^8 = 256 > 64$


## parity

``` cpp
bool parity(unsigned n);
bool parityll(unsigned long long n);
bool parityTable(unsigned n);
bool parityTablell(unsigned long long n);
bool parityMIT(unsigned n);
bool parityMITll(unsigned long long n);
```

All above are sightly slow than __builtin_parity and __builtin_parityll

## acos, asin
