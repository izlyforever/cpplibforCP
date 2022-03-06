#pragma once
#include "fft.hpp"
#include "mod.hpp"
#include "poly.hpp"

// do not use it if T = ModLL
template<typename T>
class PolyBaseFFT : public PolyBase<T> {
 protected:
  PolyBaseFFT mul(const PolyBaseFFT& rhs) const {
    int tot = std::max(1, int((int)this->size() + rhs.size() - 1));
    int sz = 1 << std::__lg(tot * 2 - 1);
    // Must be split to ensure accuracy (or use skill '3 times to 2 times')
    auto A1(*this), A2(*this), B1(rhs), B2(rhs);
    static constexpr int bit = 15, msk = (1LL << bit) - 1;
    for (auto& x : A1) x = int(x) & msk;
    for (auto& x : A2) x = int(x) >> bit;
    for (auto& x : B1) x = int(x) & msk;
    for (auto& x : B2) x = int(x) >> bit;
    std::vector<std::complex<double>> A(sz), B(sz), C(sz);
    for (int i = 0, tSize = (int)this->size(); i < tSize; ++i) {
      A[i] = std::complex<double>((double)A1[i], (double)A2[i]);
    }
    for (int i = 0, rSize = rhs.size(); i < rSize; ++i) {
      B[i] = std::complex<double>((double)B1[i], (double)B2[i]);
    }
    FFT::dft(A); FFT::dft(B);
    C[0] = conj(B[0]);
    for (int i = 1; i < sz; ++i) C[i] = conj(B[sz - i]);
    for (int i = 0; i < sz; ++i) B[i] *= A[i];
    for (int i = 0; i < sz; ++i) C[i] *= A[i];
    FFT::idft(B); FFT::idft(C);
    std::vector<T> ans(tot), A1B2(tot), A1B1(tot);
    // It will lose a lot of precision
    for (int i = 0; i < tot; ++i) {
      A1B2[i] = llround(B[i].imag());
      A1B1[i] = llround(B[i].real() * 0.5 + C[i].real() * 0.5);
      ans[i] = llround(C[i].real() * 0.5 - B[i].real() * 0.5);
    }
    for (auto& x : ans) x <<= bit;
    for (int i = 0; i < tot; ++i) ans[i] += A1B2[i];
    for (auto& x : ans) x <<= bit;
    for (int i = 0; i < tot; ++i) ans[i] += A1B1[i];
    return PolyBaseFFT(std::move(ans));
  }
 public:
  using PolyBase<T>::PolyBase;
};
const constexpr int FFTM = 1e9 + 7;
using PolyFFT = Poly<PolyBaseFFT<MInt<FFTM>>, MInt<FFTM>>;
using PolyFFTDynamic = Poly<PolyBaseFFT<ModInt>, ModInt>;
// The following are not recommented
using PolyFFTLL = Poly<PolyBaseFFT<ModLL>, ModLL>;


// $O(\sqrt{n} \log^2 n)$ base on multi-evaluation
int factorialS(int n, int p) {
  if (n >= p) return 0;
  if (n <= 1) return 1;
  ModInt::setMod(p);
  if (n > p - 1 - n) {
    int ans = ModInt(factorialS(p - 1 - n, p)).inv();
    return (p - n) & 1 ? p - ans : ans;
  }
  ModInt::setMod(p);
  int sn = std::sqrt(n);
  auto A = PolyFFTDynamic::prod(sn);
  std::vector<ModInt> x;
  x.reserve(n  / sn);
  for (int i = sn; i <= n; i += sn) x.emplace_back(i - sn + 1);
  auto y = A.evals(x);
  ModInt r(1);
  for (auto t : y) r *= t;
  for (int i = n / sn * sn + 1; i <= n; ++i) r *= ModInt::raw(i);
  return r;
}
