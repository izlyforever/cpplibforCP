#pragma once
#include <bits/stdc++.h>
using LL = long long;

// slow than __builtin_ctz and __builtin_ctzll but funny, you may use -Wno-narrowing when complier
int ctz32(unsigned x) {
  union {
    float f;
    unsigned i;
  } v = {.f = x & ~x + 1};
  return (v.i >> 23) - 127U;
}
int ctz64(unsigned long long x) {
  union {
    double f;
    unsigned long long i;
  } v = {.f = x & ~x + 1};
  return (v.i >> 52) - 1023ULL;
}
// https://xr1s.me/2018/08/23/gcc-builtin-implementation/

// MIT HAKMEM: about two times faster than __builtin_popcount()
int bitCount(unsigned n) {
  unsigned tmp = n - ((n >> 1) & 033333333333U) - ((n >> 2) & 011111111111U);
  return ((tmp + (tmp >> 3)) & 030707070707U) % 63U;
}

// MIT HAKMEM: about two times faster than __builtin_popcountll(), run with 64bit
int bitCountll(unsigned long long n) {
  unsigned long long tmp = n - ((n >> 1) & 0x7777777777777777ULL)
                             - ((n >> 2) & 0x3333333333333333ULL)
                             - ((n >> 3) & 0x1111111111111111ULL);
  return ((tmp + (tmp >> 4)) & 0x0f0f0f0f0f0f0f0fULL) % 255ULL;
}
// https://www.cnblogs.com/lukelouhao/archive/2012/06/12/2546267.html

// faster than bitCount
int bitCountTable(unsigned n) { 
  static int table[256] =  { 
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8, 
  }; 
  return table[n & 0xffUL] + table[(n >> 8) & 0xffUL] +
         table[(n >> 16) & 0xffUL] + table[n >> 24];
}
// slow than bitCountll
int bitCountTableLL(unsigned long long n) {
  return bitCountTable(n >> 32) + bitCountTable(n & 0xffffffffULL);
}
// https://www.cnblogs.com/graphics/archive/2010/06/21/1752421.html

// All below are sightly slow than __builtin_parity and __builtin_parityll
bool parity(unsigned n) { 
  n = n ^ n >> 16;
  n = n ^ n >> 8;
  n = n ^ n >> 4;
  n = n ^ n >> 2;
  return (n ^ n >> 1) & 1U;
}
bool parityll(unsigned long long n) { // slow than parityMIT
  n = n ^ n >> 32;
  n = n ^ n >> 16;
  n = n ^ n >> 8;
  n = n ^ n >> 4;
  n = n ^ n >> 2;
  return (n ^ n >> 1) & 1U;
}
bool parityTable(unsigned n) { // slow than __builtin_parity
  static bool table[256] =  { 
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 
  };
  n = n ^ n >> 16;
  return table[(n ^ n >> 8) & 0xffU];
}
bool parityTablell(unsigned long long n) { // slow than __builtin_parityll
  static bool table[256] =  { 
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 
  };
  n = n ^ n >> 32;
  n = n ^ n >> 16;
  return table[(n ^ n >> 8) & 0xffULL];
}
bool parityMIT(unsigned n) {  // slow than parity
  n = (n ^ n >> 1) & 0x55555555U;
  return (((n ^ n >> 2) & 0x11111111U) % 15U) & 1U;
}
bool parityMITll(unsigned long long n) {
  n = (n ^ n >> 1 ^ n >> 2) & 01111111111111111111111ULL;
  return (((n ^ n >> 3) & 0101010101010101010101ULL) % 63ULL) & 1U;
}

// Handbook of Mathematical Functions by M. Abramowitz and I.A. Stegun, Ed.
// Absolute error <= 6.7e-5
float acosFast(float x) {
  bool flag = (x < 0);
  x = abs(x);
  float now = sqrt(1.0 - x) * (((0.0742610f - 0.0187293f * x) * x - 0.2121144f) * x + 1.5707288f);
  return flag ? 3.14159265358979f - now : now;
}
// Absolute error <= 6.7e-5
float asinFast(float x) {
  bool flag = (x < 0);
  x = abs(x);
  float now = sqrt(1.0 - x) * (((0.0742610f - 0.0187293f * x) * x - 0.2121144f) * x + 1.5707288f);
  return flag ? now - 1.5707963267949f : 1.5707963267949f - now;
}
// https://developer.download.nvidia.cn/cg
