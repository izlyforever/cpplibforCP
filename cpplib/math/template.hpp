#pragma once
#include <bits/stdc++.h>
#include "mod.hpp"

template<typename T, typename T2>
using TwiceT = std::enable_if_t<sizeof(T) * 2 == sizeof(T2)>;

template<typename T>
using IntLongT = std::enable_if_t<
    std::is_same_v<int32_t, T>  ||
    std::is_same_v<uint32_t, T> ||
    std::is_same_v<int64_t, T>  ||
    std::is_same_v<uint64_t, T> ||
    std::is_same_v<__int128, T> ||
    std::is_same_v<__uint128, T>;


template<typename T>
using ArithT = std::enable_if_t<std::is_arithmetic_v<T>>;

template<typename T>
using SignedT = std::enable_if_t<std::is_signed_v<T> || std::is_same_v<__int128, T>>;

template<typename T>
using UnsignedT = std::enable_if_t<std::is_unsigned_v<T> || std::is_same_v<__uint128, T>>;

// valT
template<class T>
struct is_MInt : std::false_type {};

template<int M>
struct is_MInt<MInt<M>> : std::true_type {};

template<class T>
inline constexpr bool is_mint_v = is_MInt<T>::value;

template<typename T>
using ModT = std::enable_if_t<std::is_same_v<ModLL, T> || std::is_same_v<ModInt, T> || is_mint_v<T>>;
