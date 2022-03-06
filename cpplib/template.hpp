#pragma once
#include <bits/stdc++.h>

template<typename T, typename T2>
using TwiceT = std::enable_if_t<sizeof(T) * 2 == sizeof(T2)>;

template<typename T>
using ArithT = std::enable_if_t<std::is_arithmetic_v<T>>;

template<typename T>
using SignedT = std::enable_if_t<std::is_signed_v<T> || std::is_same_v<__int128_t, T>>;

template<typename T>
using UnsignedT = std::enable_if_t<std::is_unsigned_v<T> || std::is_same_v<__uint128_t, T>>;

template<typename T>
using IntegerT = std::enable_if_t<std::is_integral_v<T> || std::is_same_v<__int128_t, T> || std::is_same_v<__uint128_t, T>>;
