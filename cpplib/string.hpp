#pragma once
#include <bits/stdc++.h>

// 小写字母建 Trie，其它字符同理修改 charToInt 函数即可
class Trie {
	using Node = std::array<int, 26>;
	std::vector<Node> nxt;
	// 0 表示没有，1 表示有且没被访问过，2 表示有且被访问过
	std::vector<int> tag;
	void addNode(int fa, int c) {
		nxt[fa][c] = nxt.size();
		nxt.emplace_back(Node());
		tag.emplace_back(0);
	}
	int charToInt(char x) { return x - 'a';}
public:
	Trie() : nxt(1), tag(1) {}
	void insert(std::string s) {
		int p = 0;
		for (auto x : s) {
			int c = x - 'a';
			if (nxt[p][c] == 0) addNode(p, c);
			p = nxt[p][c]; 
		}
		tag[p] = 1;
	}
	int find(std::string s) {
		int p = 0;
		for (auto x : s) {
			int c = charToInt(x);
			if (nxt[p][c] == 0) return 0;
			p = nxt[p][c];
		}
		if (tag[p] != 1) return tag[p];
		tag[p] = 2;
		return 1;
	}
};
// 模板例题：https://www.luogu.com.cn/problem/P2580

// 01 Trie 求异或最大值
class Trie01 {
	using Node = std::array<int, 2>;
	std::vector<Node> ch;
	void addNode(int fa, int c) {
		ch[fa][c] = ch.size();
		ch.emplace_back(Node());
	}
public:
	Trie01() : ch(1) {}
	void insert(int x) {
		for (int i = 30, p = 0; i >= 0; --i) {
			int c = (x >> i) & 1;
			if (ch[p][c] == 0) addNode(p, c);
			p = ch[p][c];
		}
	}
	int getMax(int x) {
		int r = 0;
		for (int i = 30, p = 0; i >= 0; --i) {
			int c = (x >> i) & 1;
			if (ch[p][c ^ 1]) {
				p = ch[p][c ^ 1];
				r |= (1 << i);
			} else {
				p = ch[p][c];
			}
		}
		return r;
	}
	int getAns(const std::vector<int> &a) {
		int r = 0;
		for (auto x : a) {
			insert(x);
			r = std::max(r, getMax(x));
		}
		return r;
	}
};
// 典型例题：https://www.luogu.com.cn/problem/P4551

// 01-Trie (Fusion Tree) 求异或和（支持修改，全局加 1，暂不支持合并）这里求异或和，可以从最低位往最高位建树节省空间。
class FusionTree {
	using Node = std::array<int, 4>;
	std::vector<Node> ch;
	#define lsonFT ch[p][0]
	#define rsonFT ch[p][1]
	// ch[p][2] 表示以 p 为根的子树的大小
	// ch[p][3] 表示以 p 为根的子树的异或值
	void addNode(int p, int c) {
		ch[p][c] = ch.size();
		ch.emplace_back(Node());
	}
	void pushUp(int p) {
		ch[p][3] = 0;
		if (lsonFT) ch[p][3] ^= (ch[lsonFT][3] << 1);
		if (rsonFT) ch[p][3] ^= (ch[rsonFT][3] << 1) | (ch[rsonFT][2] & 1);
	}
	// 注意这里 ch[lson][2] = ch[p][2] - ch[rsonFT] 是延迟更新的。
	void insert(int p, int x) {
		++ch[p][2];
		if (!x) return;
		if (!ch[p][x & 1]) addNode(p, x & 1);
		insert(ch[p][x & 1], x >> 1);
		pushUp(p);
	}
	void erase(int p, int x) {
		--ch[p][2];
		if (!x) return;
		erase(ch[p][x & 1], x >> 1);
		pushUp(p);
	}
	void addAll(int p) {
		if (!ch[p][2]) return;
		if (rsonFT) addAll(rsonFT);
		// 为了进位，先补 0，补 0 的时候记得更新 ch[lsonFT][2]（它延迟更新了）
		if (!lsonFT) addNode(p, 0);
		ch[lsonFT][2] = ch[p][2] - (rsonFT ? ch[rsonFT][2] : 0);
		std::swap(lsonFT, rsonFT);
		pushUp(p);
	}
public:
	FusionTree() : ch(1) {}
	void insert(int x) {
		insert(0, x);
	}
	void erase(int x) {
		erase(0, x);
	}
	void addAll() {
		addAll(0);
	}
	int getVal() {
		return ch[0][3];
	}
};
// 例题：https://www.luogu.com.cn/problem/P6018

// 前缀函数 O(n) 算法
std::vector<int> prefixFunction(std::string s) {
	int n = s.size();
	std::vector<int> p(n);
	for (int i = 1; i < n; ++i) {
		int j = p[i - 1];
		while (j > 0 && s[i] != s[j]) j = p[j - 1];
		if (s[i] == s[j]) ++j;
		p[i] = j;
	}
	return p;
}

// 基于前缀函数的 KMP 算法，返回所有匹配在 t 的首位置
std::vector<int> kmp(std::string s, std::string t) {
	std::vector<int> ans;
	int n = s.size(), m = t.size();
	if (n > m) return ans;
	auto p = prefixFunction(s);
	for (int i = 0, j = 0; i < m; ++i) {
		while (j > 0 && s[j] != t[i]) j = p[j - 1];
		if (s[j] == t[i] && ++j == n) ans.emplace_back(i - n + 1);
	}
	return ans;
}

// 返回 长度为 i 的前缀出现的次数
std::vector<int> countPrefix(std::string s) {
	auto p = prefixFunction(s);
	int n = s.size();
	std::vector<int> ans(n + 1);
	for (auto x : p) ++ans[x];
	for (int i = n - 1; i > 0; --i) ans[p[i - 1]] += ans[i];
	for (int i = 0; i <= n; ++i) ++ans[i];
	return ans;
}
// 返回 s 长度为 i 的前缀在 t 中出现的次数
std::vector<int> countPrefix(std::string s, std::string t) {
	auto p = prefixFunction(s);
	int n = s.size(), m = t.size();
	std::vector<int> ans(n + 1);
	for (int i = 0, j = 0; i < m; ++i) {
		while (j > 0 && t[i] != s[j]) j = p[j - 1];
		if (t[i] == s[j]) ++j;
		++ans[j];
	}
	++ans[0];
	for (int i = n; i > 0; --i) ans[p[i - 1]] += ans[i];
	return ans;
}
// 例题：https://codeforces.com/problemset/problem/432/D

// $z[i]$ 表示 s 和 $s[i, n - 1]$ 的最长公共前缀，约定 $z[0] = 0$
std::vector<int> zFunction(std::string s) {
	int n = s.size();
	std::vector<int> z(n);
	for (int i = 1, l = 0, r = 0; i < n; ++i) {
		if (i <= r && z[i - l] < r - i + 1) {
			z[i] = z[i - l];
		} else {
			int j = std::max(0, r - i + 1);
			while (i + j < n && s[j] == s[i + j]) ++j;
			z[i] = j;
			if (i + j - 1 > r) {
				l = i;
				r = i + j - 1;
			}
		}
	}
	return z;
}

// 基于 Z 函数的 KMP 算法，返回所有匹配在 t 的首位置
std::vector<int> kmpZ(std::string s, std::string t) {
	std::vector<int> ans;
	int n = s.size(), m = t.size();
	if (n > m) return ans;
	auto z = zFunction(s);
	for (int i = 0, l = 0, r = -1; i < m; ++i) {
		if (i > r || z[i - l] >= r - i + 1) {
			int j = std::max(0, r - i + 1);
			while (j < n && i + j < m && s[j] == t[i + j]) ++j;
			if (j == n) ans.emplace_back(i);
			if (i + j - 1 > r) {
				l = i;
				r = i + j - 1;
			}
		}
	}
	return ans;
}

class Automaton {
	static inline constexpr int CHAR = 26;
	using Node = std::array<int, CHAR>;
	std::vector<Node> nxt;
	std::vector<int> cnt, fail, last;
	int charToInt(char x) { return x - 'a';}
	void addNode(int fa, int c) {
		nxt[fa][c] = nxt.size();
		nxt.emplace_back(Node());
		cnt.emplace_back(0);
		fail.emplace_back(0);
		last.emplace_back(0);
	}
public:
	Automaton() : nxt(1), cnt(1), fail(1), last(1) {}
	void insert(std::string s) {
		int p = 0;
		for (auto x : s) {
			int c = charToInt(x);
			if (nxt[p][c] == 0) addNode(p, c);
			p = nxt[p][c];
		}
		++cnt[p];
	}
	void build() {
		std::queue<int> Q;
		for (int c = 0; c < CHAR; ++c) {
			if (nxt[0][c]) Q.push(nxt[0][c]);
		}
		while (!Q.empty()) {
			int p = Q.front(); Q.pop();
			for (int c = 0; c < CHAR; ++c) {
				if (int &q = nxt[p][c]; q != 0) {
					fail[q] = nxt[fail[p]][c];
					Q.push(q);
					// 用作模式匹配时计数的优化
					last[q] = cnt[fail[q]] ? fail[q] : last[fail[q]];
				} else {
					q = nxt[fail[p]][c];
				}
			}
		}
	}
	// 具体写法见题目要求
	int query(std::string s) {
		int p = 0, r = 0;
		auto add = [&](int & x) {
			r += x; x = 0;
		};
		for (auto x : s) {
			int c = charToInt(x);
			p = nxt[p][c];
			if (cnt[p]) add(cnt[p]);
			int q = p;
			while (last[q]) {
				q = last[q];
				if (cnt[q]) add(cnt[q]);
			}
		}
		return r;
	}
};
// 模板例题：https://www.luogu.com.cn/problem/P3808

// O(N) 复杂度计算后缀数组的 SA-IS 算法
// 请确保最后一个元素为 0，且 a 中其它元素都是正整数，且最大值较小。
std::vector<int> SAIS(std::vector<int> a) {
	enum TYPE {L, S};
	int n = a.size() - 1, mx = *std::max_element(a.begin(), a.end()) + 1;
	std::vector<int> SA(n + 1, -1);
	std::vector<int> bucket(mx), lbucket(mx), sbucket(mx);
	for (auto x : a) ++bucket[x];
	for (int i = 1; i < mx; ++i) {
		bucket[i] += bucket[i - 1];
		lbucket[i] = bucket[i - 1];
		sbucket[i] = bucket[i] - 1;
	}
	// 确定 L, S 类型以及 * 型的位置
	std::vector<TYPE> type(n + 1);
	type[n] = S;
	for (int i = n - 1; i >= 0; --i) {
		type[i] = (a[i] < a[i + 1] ? S : (a[i] > a[i + 1] ? L : type[i + 1]));
	}
	// 诱导排序(从 * 型诱导到 L 型、从 L 型诱导到 S 型)
	auto inducedSort = [&]() {
		for (int i = 0; i <= n; ++i) {
			if (SA[i] > 0 && type[SA[i] - 1] == L) {
				SA[lbucket[a[SA[i] - 1]]++] = SA[i] - 1;
			}
		}
		for (int i = 1; i < mx; ++i) {
			sbucket[i] = bucket[i] - 1;
		}
		for (int i = n; i >= 0; --i) {
			if (SA[i] > 0 && type[SA[i] - 1] == S) {
				SA[sbucket[a[SA[i] - 1]]--] = SA[i] - 1;
			}
		}
	};
	// 首先根据诱导排序给出 LMS-prefix 的排序
	std::vector<int> pos;
	for (int i = 1; i <= n; ++i) {
		if (type[i] == S && type[i - 1] == L) {
			pos.emplace_back(i);
		}
	}
	for (auto x : pos) SA[sbucket[a[x]]--] = x;
	inducedSort();
	// 根据 LMS-prefix 的排序给出 LMS-substring 的命名，即得到 S1
	auto isLMSchar = [&](int i) {
		return i > 0 && type[i] == S && type[i - 1] == L;
	};
	auto equalSubstring = [&](int x, int y) {
		do {
			if (a[x] != a[y]) return false;
			++x; ++y;
		} while (!isLMSchar(x) && !isLMSchar(y));
		return a[x] == a[y];
	};
	// 注意到因为 LMS-prefix 排序会导致仅有相邻的 LMS-substring 才可能相等
	std::vector<int> name(n + 1, -1);
	int lx = -1, cnt = 0;
	bool flag = true; // 表示无相同的 LMS-substring
	for (auto x : SA) if (isLMSchar(x)) {
		if (lx >= 0 && !equalSubstring(lx, x)) ++cnt;
		if (lx >= 0 && cnt == name[lx]) flag = false;
		name[x] = cnt;
		lx = x;
	}
	std::vector<int> S1;
	for (auto x : name) if (x != -1) S1.emplace_back(x);
	auto getSA1 = [&]() {
		int n1 = S1.size();
		std::vector<int> SA1(n1);
		for (int i = 0; i < n1; ++i) SA1[S1[i]] = i;
		return SA1;
	};
	auto SA1 = flag ? getSA1() : SAIS(S1);
	// 再次诱导排序，根据 S1 的排序得到 SA
	lbucket[0] = sbucket[0] = 0;
	for (int i = 1; i < mx; ++i) {
		lbucket[i] = bucket[i - 1];
		sbucket[i] = bucket[i] - 1;
	}
	std::fill(SA.begin(), SA.end(), -1);
	// 这里是逆序扫描 SA1，因为 SA 中 S 型桶是倒序的
	for (int i = SA1.size() - 1; i >= 0; --i) {
		SA[sbucket[a[pos[SA1[i]]]]--] = pos[SA1[i]];
	}
	inducedSort();
	return SA;
}

std::vector<int> SAIS(const std::string &s) {
	// s 的字符集为小写字，则可使用下面函数。
	// auto f = [](char x) -> int { return int(x - 'a') + 1;};
	auto f = [](char x) -> int { return int(x) + 1;};
	std::vector<int> a;
	for (auto c : s) a.emplace_back(f(c));
	a.emplace_back(0);
	auto sa = SAIS(a);
	return std::vector<int>(sa.begin() + 1, sa.end());
}

// 最小表示法
template<typename T>
int minPresent(std::vector<T>& a) {
	int k = 0, i = 0, j = 1, n = a.size();
	while (k < n && i < n && j < n) {
		if (a[(i + k) % n] == a[(j + k) % n]) {
			++k;
		} else {
			a[(i + k) % n] > a[(j + k) % n] ? i += k + 1 : j += k + 1;
			if (i == j) ++i;
			k = 0;
		}
	}
	return std::min(i, j);
}

// Lyndon decomposition using Duval algorithm
std::vector<std::string> duval(const std::string &s) {
	std::vector<std::string> r;
	int n = s.size(), i = 0;
	while (i < n) {
		int j = i + 1, k = i;
		while (j < n && s[k] <= s[j]) {
			if (s[k] < s[j]) k = i;
			else ++k;
			++j;
		}
		while (i <= k) {
			r.emplace_back(s.substr(i, j - k));
			i += j - k;
		}
	}
	return r;
}
// 模板例题: https://www.luogu.com.cn/problem/P6114

// 仅仅处理奇数长回文串，这个实现像极了 Z-函数
std::vector<int> Manacher(std::string s) {
	int n = s.size();
	std::vector<int> d(n);
	for (int i = 0, l = 0, r = -1; i < n; ++i) {
		int k = i > r ? 1 : std::min(d[l + r - i], r - i);
		while (k <= i && i + k < n && s[i - k] == s[i + k]) ++k;
		d[i] = k--;
		if (i + k > r) {
			l = i - k;
			r = i + k;
		}
	}
	return d;
}
// 通过中间加相同特殊符号统一成仅需考虑奇数长回文串
