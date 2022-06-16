#include "bits/stdc++.h"
#include "bigInt.hpp"

int main() {
  // freopen("in", "r", stdin);
  std::string a, b;
  std::cin >> a >> b;
  BigInt10 A(a), B(b);
  std::cout << A + B << '\n'; // https://www.luogu.com.cn/problem/P1601
  std::cout << A - B << '\n'; // // https://www.luogu.com.cn/problem/P2142
  std::cout << A * B << '\n'; // https://www.luogu.com.cn/problem/P1303
  std::cout << A / B << '\n'; // https://loj.ac/p/164
  return 0;
}
