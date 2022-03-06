# primary.hpp

`primary.hpp` is the foundation of `math.hpp`.

## [basic.hpp](basic.md)

## [mod.hpp](mod.md)

## fft.hpp

> fast Fourier transform

It contains `dft` and `idft` in namespace `FFT`.

> The size of $a$ is a pow of $2$ before you use `dft(a)` or `idft(a)`

## ntt.hpp

> number theory transform

It contains `dft` and `idft` in template classs `NTT`.

> The template $M$ should be NTT-friendly, and size of $a$ is a pow of $2$ less that $10^6$ before you use `dft(a)` or `idft(a)`
> 

## fmt.hpp

> fast Mobius transform
