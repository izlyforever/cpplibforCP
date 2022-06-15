# BigInt

A high performance and easy to use class for CPP17


## Notes

- `vector` is used to restore the number. 
- `ngt`: `true` if the number is negetive, `false` otherwise.
- NTT is used to give a $O(n \log n)$ algorithm for multiplication.
- $\frac{a}{b}$ is calculated with the ideal of compute $\frac{a \frac{10^{2m}}{b}}{10^{2m}}$, where $m$ is the length of $b$.
- Other basic methods are trival.
