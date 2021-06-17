#pragma once
#include "mod.hpp"
#include "ntt.hpp"
#include "poly.hpp"

class PolyBaseMFT4 : public PolyBase<ModLL> {
public:
	static inline constexpr int M0 = 595591169, M1 = 645922817, M2 = 897581057, M3 = 998244353;
	static inline constexpr LL M01 = 1LL * M0 * M1, M23 = 1LL * M2 * M3;
	static inline constexpr __int128 M0123 = __int128(M01) * M23;
	static inline constexpr int t01 = 538269027, t23 = 415935157;
	static inline constexpr LL t0123 = 341204425684314487LL;
	static inline NTT<M0> ntt0;
	static inline NTT<M1> ntt1;
	static inline NTT<M2> ntt2;
	static inline NTT<M3> ntt3;
	using PolyBase<ModLL>::PolyBase;
	PolyBaseMFT4 (const PolyBase<ModLL> &x) : PolyBase<ModLL>(x) {}
protected:
	static ModLL crt(int a0, int a1, int a2, int a3) {
		LL ans1 = a0 + LL(a1 - a0) * t01 % M1 * M0;
		LL ans2 = a2 + LL(a3 - a2) * t23 % M3 * M2;
		__int128 ans = ans1 + __int128(ans2 - ans1) * t0123 % M23 * M01;
		if (ans < 0) ans += M0123;
		return ModLL(ans);
	}
	PolyBaseMFT4 mul(const PolyBaseMFT4 &rhs) const {
		int tot = std::max(1, this->size() + rhs.size() - 1);
		int sz = 1 << std::__lg(tot * 2 - 1);
		std::vector<MInt<M0>> a0(sz), b0(sz);
		std::vector<MInt<M1>> a1(sz), b1(sz);
		std::vector<MInt<M2>> a2(sz), b2(sz);
		std::vector<MInt<M3>> a3(sz), b3(sz);
		for (int i = 0, ns = this->size(); i < ns; ++i) {
			LL tmp = this->a[i];
			a0[i] = tmp; a1[i] = tmp; a2[i] = tmp; a3[i] = tmp;
		}
		for (int i = 0, ns = rhs.size(); i < ns; ++i) {
			LL tmp = rhs.a[i];
			b0[i] = tmp; b1[i] = tmp; b2[i] = tmp; b3[i] = tmp;
		}
		ntt0.dft(a0); ntt0.dft(b0);
		ntt1.dft(a1); ntt1.dft(b1);
		ntt2.dft(a2); ntt2.dft(b2);
		ntt3.dft(a3); ntt3.dft(b3);
		for (int i = 0; i < sz; ++i) {
			a0[i] *= b0[i]; a1[i] *= b1[i]; a2[i] *= b2[i]; a3[i] *= b3[i];
		}
		ntt0.idft(a0); ntt1.idft(a1); ntt2.idft(a2); ntt3.idft(a3);
		std::vector<ModLL> ans(tot);
		for (int i = 0; i < tot; ++i) ans[i] = crt(a0[i], a1[i], a2[i], a3[i]);
		return PolyBaseMFT4(ans);
	}
};
// 4-module NTT, Module can be up to 1e14 since N < 1e6 in general
using PolyMFT = Poly<PolyBaseMFT4, ModLL>;

// $O(\sqrt{n} \log n)$ (assume n < 1e12)
LL factorial(LL n, LL p) {
	if (n >= p) return 0;
	if (n <= 1) return 1;
	ModLL::setMod(p);
	if (n > p - 1 - n) {
		LL ans = ModLL(factorial(p - 1 - n, p)).inv();
		return (p - n) & 1 ? p - ans : ans; 
	}
	int s = std::sqrt(n);
	std::vector<ModLL> fac(s + 1), ifac(s + 1), inv(s + 1);
	fac[0] = inv[1] = 1;
	for (int i = 1; i <= s; ++i) fac[i] = fac[i - 1] * ModLL::raw(i);
	ifac[s] = fac[s].inv();
	for (int i = s; i > 0; --i) ifac[i - 1] = ifac[i] * ModLL::raw(i);
	for (int i = 2; i <= s; ++i) inv[i] = inv[p % i] * ModLL::raw(p - p / i);
	// compute $h(m), \cdots, h(m + cnt - 1)$ according to $h(0), h(1), \cdots, h(d)$
	auto solve = [&](std::vector<ModLL> h, ModLL m, int cnt) { // m > h.size()
		int d = h.size() - 1;
		for (int i = 0; i <= d; ++i) {
			h[i] *= ifac[i] * ifac[d - i];
			if ((d - i) & 1) h[i] = -h[i];
		}
		std::vector<ModLL> f(d + cnt);
		auto now = m - ModLL(d);
		for (int i = 0; i < d + cnt; ++i) {
			f[i] = now.inv();
			++now;
		}
		h = (PolyMFT(f) * PolyMFT(h)).a;
		h.resize(d + cnt);
		h = std::vector<ModLL>(h.begin() + d, h.end());
		now = 1;
		for (int i = 0; i <= d; ++i) now *= m - ModLL::raw(i);
		h[0] *= now;
		for (int i = 1; i < cnt; ++i) {
			now *= m + ModLL::raw(i);
			now *= (m + ModLL(i - d - 1)).inv();
			h[i] *= now;
		}
		return h;
	};
	auto invS = ModLL(s).inv();
	std::vector<ModLL> h{1, s + 1};
	for (int bit = std::__lg(s) - 1, d = 1; bit >= 0; --bit) {
		auto nh1 = solve(h, d + 1, d);
		auto nh2 = solve(h, invS * ModLL(d), 2 * d + 1);
		h.insert(h.end(), nh1.begin(), nh1.end());
		d *= 2;
		for (int i = 0; i <= d; ++i) h[i] *= nh2[i];
		if (s >> bit & 1) {
			++d;
			LL tmp = d;
			for (int i = 0; i < d; ++i, tmp += s) h[i] *= ModLL::raw(tmp);
			ModLL last(1), tj = 1LL * s * d;
			for (int i = 0; i < d; ++i) ++tj, last *= tj;
			h.emplace_back(last);
		}
	}
	ModLL ans = 1;
	for (int i = 0; i < s; ++i) ans *= h[i];
	for (LL i = 1LL * s * s + 1; i <= n; ++i) ans *= ModLL::raw(i);
	return ans;
}
// https://vjudge.net/problem/SPOJ-FACTMODP








// // don't use it, use PolyFFT or PolyBaseMFT4
// template<typename T>
// class PolyBaseMFT3 : public PolyBase<T> {
// public:
// 	static inline constexpr int M0 = 469762049, M1 = 998244353, M2 = 1004535809;
// 	static inline constexpr LL M01 = 1LL * M0 * M1;
// 	static inline constexpr int t0 = 554580198, t1 = 395249030;
// 	static inline NTT<M0> ntt0;
// 	static inline NTT<M1> ntt1;
// 	static inline NTT<M2> ntt2;
// 	using PolyBase<T>::PolyBase;
// 	PolyBaseMFT3 (const PolyBase<T> &x) : PolyBase<T>(x) {}
// protected:
// 	static T crt(int a0, int a1, int a2) {
// 		static const T m01(M01);
// 		LL x = (a0 + 1LL * (a1 - a0) * t0 % M1 * M0) % M01;
// 		if (x < 0) x += M01;
// 		LL y = (a2 - x) % M2;
// 		if (y < 0) y += M2;
// 		y = y * t1 % M2;
// 		if (y < 0) y += M01;
// 		return T(x) + T(y) * m01;
// 	}
// 	PolyBaseMFT3 mul(const PolyBaseMFT3 &rhs) const {
// 		int tot = std::max(1, this->size() + rhs.size() - 1);
// 		int sz = 1 << std::__lg(tot * 2 - 1);
// 		std::vector<MInt<M0>> a0(sz), b0(sz);
// 		std::vector<MInt<M1>> a1(sz), b1(sz);
// 		std::vector<MInt<M2>> a2(sz), b2(sz);
// 		for (int i = 0, ns = this->size(); i < ns; ++i) {
// 			int t = this->a[i];
// 			a0[i] = t; a1[i] = t; a2[i] = t;
// 		}
// 		for (int i = 0, ns = rhs.size(); i < ns; ++i) {
// 			int t = rhs.a[i];
// 			b0[i] = t; b1[i] = t; b2[i] = t;
// 		}
// 		ntt0.dft(a0); ntt0.dft(b0);
// 		ntt1.dft(a1); ntt1.dft(b1);
// 		ntt2.dft(a2); ntt2.dft(b2);
// 		for (int i = 0; i < sz; ++i) {
// 			a0[i] *= b0[i]; a1[i] *= b1[i]; a2[i] *= b2[i];
// 		}
// 		ntt0.idft(a0); ntt1.idft(a1); ntt2.idft(a2);
// 		std::vector<T> ans(tot);
// 		for (int i = 0; i < tot; ++i) ans[i] = crt(a0[i], a1[i], a2[i]);
// 		return PolyBaseMFT3(ans);
// 	}
// };
// const constexpr int MTTM = 1e9 + 7;
// using PolyMFT3 = Poly<PolyBaseMFT3<MInt<MTTM>>, MInt<MTTM>>;
// using PolyMFTDynamic = Poly<PolyBaseMFT3<ModInt>, ModInt>;