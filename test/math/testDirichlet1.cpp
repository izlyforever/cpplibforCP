#include <bits/stdc++.h>
#define clog(x) std::clog << (#x) << " is " << (x) << '\n';
using LL = long long;
#include "../cpplib/math/numberTheory.hpp"

int main() {
	//freopen("in", "r", stdin);
	std::cin.tie(nullptr)->sync_with_stdio(false);
	std::vector<int> a{0, 1, 10, 100, 1000, 10000, 100000};
	std::vector<int> b{0, 1, 1, 1, 1, 1, 1};
	std::vector<int> c{0, 1, -1, -1, 0, -1, 1}; // mu

	DirichletProduct<int>::setLen(int(a.size() - 1));
	DirichletProduct A(a), B(b), C(c);
	auto D = A;
	std::cout << "Fast mobious Test: \n";
	std::cout << A * B << '\n';
	A.mobious();
	std::cout << A << '\n';
	A.mobiousInv();

	std::cout << "\nFast mobious invserse Test: \n";
	std::cout << D * C << '\n';
	D.mobiousInv();
	std::cout << D << '\n';	

	std::cout << "\nFast transpose mobious Test: \n";
	D = A;
	std::cout << A.mulT(B) << '\n';
	A.mobiousT();
	std::cout << A << '\n';
	A.mobiousInvT();

	std::cout << "\nFast transpose mobious invserse Test: \n";
	std::cout << D.mulT(C) << '\n'; // overflow for int occur
	D.mobiousInvT();
	std::cout << D << '\n';	
	
	std::cout << "\nSee origin A: \n";
	std::cout << A;
	
	return 0;
}