此文档仅为补充，这里就放一个大致目录和一些核心内容的说明（有些英文说太麻烦），可以在我的 [cnblog](https://www.cnblogs.com/izlyforever/p/templateOfCpp.html) 里找到很多核心代码的原理

# 文档

- 使用：编译器需支持 C++17，强烈建议开启 O2 优化
- 分类：数学，数据结构，字符串，图论，几何，杂类
- 约定：以 S 为后缀的算法都是 慢且简单 的算法，？表示暂未实现。
- 下标：默认以 0 开头，目前只有**树状数组和一些树算法**为了方便起见从 1 开头。

## 思想

- 动态规划
- 二分
- 倍增
- 分块
- 分治
- [水涨船高]((https://codeforces.com/blog/entry/58316))
- 决策单调性优化
- Meet in Middle
- Small to Large

## 数学

### 基础模块：primary.hpp

- 模快速幂
- 向上取整和向下取整
- int128 读写（快读，int, long long 也可以使用）
- 二进制快速 gcd，拓展 gcd
- 中国剩余定理 CRT
- 常规二项式系数和模二项式系数（单例）
- Lagrange 插值
- 模自然数方幂和 $O(k)$ 算法
- 模 $N \times N$ 矩阵乘法类（缓存优化）
- MEX（集合中不出现的最小的自然数）
- 快速排序（没带随机选择，别用）
- 预处理最小素因子
- BerlekampMassey（用来找最短递推公式）

### mod.hpp

- MInt: 这个是类模板
- ModInt
- ModLL

### FFT.hpp

### NTT.hpp

### FMT.hpp

快速 Mobius 变换，叫这个名字是因为跟数论函数的 Mobius 变换形式上一致

### 初等数论：numberTheory.hpp

- 快速 $O(n \log \log n)$ 素数筛和慢速线性筛
- 快速计算 $\pi(x)$
- 快速计算第 n 个素数（从 1 开始标号，p[1] = 2，p[0] 无意义）
- Euler 线性筛
- Mobius 线性筛
- 拓展欧拉定理（结论十分简单，证明需要分解素因式）
- min_25 筛 $O(n^{\frac{2}{3}})$ 算法求 Euler 函数前缀和，Mobius（绝对值）前缀和（内含 整除分块）
- 最小素因子线性筛
- 预处理素因子个数（算重/不算重）
- 素因子分解（算重/不算重）
- 求原根
- 大素数 Miller-Rabin 概率判别法
- Pollard-Pho 大整数最大最小素因子分解
- 模素数取 log（BabyStepGaintStep）
- 模素数开根号 $O(log^2 p)$
- 模素数开根号 $O(log p)$ 的 Cipolla 算法
- 模偶素数幂开根号？
- 模奇素数幂开根号？
- 模任意数开根号（先因式分解，看作模素数方开根号，再 CRT 整合）？

### BigInt ？

### 杂类

- 快速暴力 $n$ 个集合中选 $k$ 个，二进制为 1 的表示选择
- Fibonacci 数列
- floorSum：$\displaystyle \sum_{i = 0}^{n - 1} \lfloor \frac{a \cdot i + b}{m} \rfloor$
- sumNum：$\displaystyle \sum_{\sum c_i x_i = m} \frac{(\sum x_i)!}{\prod (x_i !)}$
- decInc: 每次可选择 n 减一 或 m 加一，使得 m 是 n 的倍数的最小次数
- Gauss 消元法浮点数版
- 模 Gauss 消元法
- 线性规划之单纯形算法
- 任意模数多项式乘法 $O(n^{\log_2 3})$ 的 Karatsuba 算法（包括并行版）
- 线性规划
- FirstInRange：求最小的 $x$ 使得 $l \leq a x \mod m \leq r$。类似 exgcd 的处理：求最小非负整数 $x$ 使得 $l \leq ax - m y \leq r$ 等价于 $l \leq (az - m)y - a(yz - x) \leq r$，注意到我们要始终保持 $a < m$，因此当 $2a > m$ 时需要特判一下。转化成 $m - r \leq (m - a) x - m(x - y - 1) \leq m - l$


### 多项式（[多项式全家桶](https://www.luogu.com.cn/training/3015#problems) 已全部 AC）

- 仅包含乘法的四大多项式底层基类分别为：PolyBaseNTT, PolyBaseMFT3(弃用，被后面两个淘汰了), PolyBaseMFT4, PolyBaseFFT
- PolyBaseNTT：基于固定的 NTT-friendly（原根一般为 3）模数快速数论变化（看具体题目，一般为 998244353）
- PolyBaseMFT3：基于三个固定的 NTT-friendly 且原根为 3 的三模数（469762049, 998244353, 1004535809），利用 crt 求解任意模数多项式乘法（已被淘汰，请勿使用）
- PolyBaseMFT4：基于四个固定的 NTT-friendly 且原根为 3 的四模数（595591169, 645922817, 897581057, 998244353），利用 crt 求解任意模数多项式乘法
- PolyBaseFFT：基于 FFT 求解任意模数多项式乘法（需要注意精度）
- 通过模板继承拓展得到全面的多项式类 Poly (加减乘除余，转置乘法，求导，积分，指数，对数，求逆，开方，一点求值，多点求值，快速幂模，内积，一个首一多项式的次方模 $x^n$ 先取对数乘以次数再取指数得到，三角函数，反三角函数)，这个过程学到了很多东西
- 多项式静态函数：$O(n \log^2 n)$ 计算 $\sum_{i = 1}^n \frac{a_i}{1 - b_i}$
- 多项式静态函数：$O(k \log k \log n)$ 求 $k$ 阶常系数递推公式的第 $n$ 项
- 多项式静态函数：模自然数方幂和 $O(k \log k)$ 得到前 $k$ 个答案
- 多项式静态函数：Lagrange 插值：先分治求 $g(x) = \prod(x - x_i)$，再求 $g'(x)$ 在 $x$ 处的多点求值，再分治即可。
- 求阶乘 $n! \mod p$：基于多点求值 $O(\sqrt{n} \log^2 n)$ 求 $\sqrt{n}$ 个点之后暴力
- 求阶乘 $n! \mod p$：min_25 用点求点 $O(\sqrt{n} \log n)$ 求 $\sqrt{n}$ 个点之后暴力

> 无运算的多项式底层基类：PolyBase（standard 在取余时，特别重要不可省略）

#### 使用准则

- 多项式项数 $N < 4 \cdot 10^6$
- $M$ 要是超了 int，那就只能用 ModLL 版本 4 模数 Poly
- 否则，要是 $M$ 不固定就用使用 ModInt 的 FFT 版 Poly
- 否则，当 $M$ 为固定的 NTT-friendly 素数时，使用 NTT 版 Poly
- 否则，使用 MInt 的 FFT 版 Poly

#### 极简版多项式模板（polyS）

由于多项式模板一直扩展，动则 1000+ 行，实在有点搞，所以就搞了一个极简版的。

####  4 次 FFT 原理

我们本来要求 $A \dot B$，然后我们把它拆成 $A = A_1 + 2^{d} A_2, B = B_1 + 2^d B_2$，那么 
$$
AB = A_1 B_1 + (A_1 B_2 + A_2 B_1) 2^d + A_2 B_2 2^{2d}
$$

于是我们只需计算 $(A_1 + i A_2)(B_1 + i B_2)$ 和 $(A_1 + i A_2)(B_1 - i B_2)$ 即可，但是注意到 $dft((B_1 + i B_2))[j] = \overline{dft((B_1 + i B_2))[n - j]}$ 所以本来需要 5 次 FFT，现在只需要 4 次即可。

> 其实本质上，我们可以只做 3.5 次 FFT，因为 2 次 dft 我们可以得到 $A_1, A_2, B_1, B_2$ 的 dft 值，然后我们最后只需 3 次实数版 idft 即可（算作 1.5 次）！所以总的来说是 3.5 次。但是实现的时候也没办法搞 0.5 次，可惜。

#### min_25 用点求点原理（其实可以用下降幂更简洁的处理）

学习资料：[zzqsblog](https://www.cnblogs.com/zzqsblog/p/8408691.html), [bztMinamoto](https://www.cnblogs.com/bztMinamoto/p/10661226.html)

我们令 $s = \sqrt{n}$ 然后 $\displaystyle g_{s}(x) = \sum_{i = 1}^{s}(x + i)$，我们想要得到 $g_s(0), g_s(s), \cdots g_s((s - 1)s)$ 的值。然后 $n! = \prod_{i = 0}^{s - 1} g_s(i s) \cdot \prod_{i = s^2 + 1}^n i$

现在假设我们已经得到了

$$
g_d(0), g_d(s), \cdots g_d(d s)
$$

> 一个 $d$ 次多项式由它在 $d + 1$ 个不同点的取值唯一决定（多于 d + 1 个点也可以）

1. 我们如何求

$$
g_{d + 1}(0), g_{d + 1}(s), \cdots g_{d + 1}((d + 1) s)
$$

注意到 $g_{d + 1}(x) = g_{d}(x) \cdot (x + d + 1)$ 即可 $O(d)$ 计算出 前 $d + 1$ 个，最后一个直接暴力计算即可。

2. 我们如何求

$$
g_{2d}(0), g_{2d}(s), \cdots g_{2d}(2d s)
$$

同样我们注意到 $g_{2d}(x) = g_{d}(x) \cdot g_d(x + d)$ 如果我们设 $h(i) = g_d(i s)$，（那么很关键的一点 $g_d(is + d) = g((d / s + i) s = h(d / s + i)$，卧槽， $d / s$ 在模 p 意义下得到就可以了，而且肯定大于 d，否则矛盾！）那么问题就转化成如何根据一个 $d$ 次多项式的值：$h(0), h(1), \cdots, h(d)$ 求

$$
h(d + 0), \cdots h(d + d)
$$

以及 

$$
h(d / s + 0), \cdots h(d / s + d)
$$

我们不妨对于任意的给定的 $k > d$，先求出

$$
h(k + 0), h(k + 1), \cdots h(k + d) 
$$

注意到根据 Lagrange 插值多项式

$$
\begin{aligned}
h(x) &= \sum_{i = 0}^{d} h(i) \prod_{j = 0, j \neq i}^{d} \frac{x - j}{i - j} \\ 
&= \sum_{i = 0}^d (-1)^{d - i} h(i) \binom{x}{i} \binom{x - i - 1}{d - i} \\
&=  \left(\prod_{i = x - d}^x i \right) \sum_{i = 0}^d \frac{h(i)}{i!(d - i)!(-1)^{d - i}} \cdot \frac{1}{(x - i)}
\end{aligned}
$$

注意到这里的卷积跟我们普通的卷积不一致，左边长度为 $d + 1$ 的多项式乘以右边长度为 $2d + 1$ 的多项式，然后次数为 $d, \cdots 2d$ 这 $d + 1$ 位是有效的。

1. 不能写成除以阶乘的形式，因为 x 有可能很大。
2. 为了保证 $d/s < 2 d$，我们需要使用 Wilson 定理即 $(p - 1)! = -1$

#### 分治 FFT

已知 $f_i = \sum_{j=1}^i f_{i-j} g_j$ 和 $f_0$ 求 $f$

这个显然可以分治来做，其实用生成函数推理可知 $f = \frac{f_0}{1 - g}$

#### 下降幂与点值

设 n 次多项式 $f(x) = \sum_{i = 0}^n b_i x^{\underline{i}}$，则

$$
\frac{f(m)}{m!} = \sum_{i = 0}^m b_i \frac{1}{(m - i)!} \qquad 0 \leq m \leq n
$$

因此 $EGF(f) = b e^x$，反过来也一样  $b = EGF(f) e^{-x}$（注意这里可以简单的多点求值，可以求更多的点）

> 下降幂与连续点值有 $O(n \log n)$ 的转化。而普通多项式跟连续点值却没有，可以认为普通多项式要的连续其实是类似 FFT 那样的连续。但是注意到以连续点求连续点有 $O(n \log n)$ 的做法


#### [Lagrange 反演](https://users.math.msu.edu/users/magyarp/Math880/Lagrange.pdf)

若 $f(x), g(x) \in F[[x]]$ 且 $f(g(x) = x$，则

$$
[x^n] g(x) = \frac{1}{n} [x^{-1}] \frac{1}{f(x)^n}
$$

特别地，若 $f(x) = \frac{x}{\phi(x)}$，则

$$
[x^n] g(x) = \frac{1}{n} [x^{n-1}] \phi(x)^n
$$

### 几何

- 二维凸包
- 旋转卡壳（彻底弄懂原理，之后应用此原理推广应用解决：https://www.luogu.com.cn/problem/P7164）
- 分治法求平面最短距离
- k 维偏序之 bitset 暴力优化 $O(\frac{k n^2}{w})$
- 四边形优化 DP


## 图论

- dfs 序
- Euler 序
- 最近公共祖先 LCA
- 最小生成树 Prim
- 最小树形图 LiuZhu
- 拓扑排序
- Euler 路
- Hamilton 路？（NPC 问题，下次一定）
- 带路径 Floyd
- 最短路 Dijkstra
- 最短路 BellmaFord
- 最短路 SPFA
- 连通分量之 Kosaraju 缩点
- 连通分量之 2-SAT
- 割点割边
- 有向图 S-T 最大流 Dinic $O(n^2 m)$
- 有向图 S-T 最大流的最高标号预流推进算法（HLPP） $O(n^2 \sqrt{m})$ 算法
- 无向图全局最小割 StoerWagner 算法
- 最小费用最大流（势能 Dijkstra）
- 三元环计数


### 三元环计数上界 $O(m \sqrt{m})$

首先定向之后，每条边的出度不会超过 $\sqrt{m}$（妙啊），这是因为
- 原来不超过的必然不超过
- 原来超过的只会连度大于等于它的，这个个数不会超过 $\sqrt{m}$


## 字符串

- Trie 普通字符串版
- Trie01 求两两异或最大值（可在线）
- Trie01(FusionTree) 求异或和（支持修改，全局加 1）
- 前缀函数
- 基于前缀函数的 KMP 算法
- 基于前缀函数求前缀出现次数
- Z-函数
- 基于 Z-函数的 KMP 算法
- AC 自动机
- 后缀数组计算的 O(N) 诱导排序 SA-IS 算法
- 最小表示法
- Lyndon 分解的 Duval 算法
- 处理回文的 Manacher 算法

## 数据结构

- 暴力枚举
- 纠错码
- 离散化
- 并查集 (Disjoint Union Set)
- 树状数组 (FenwickTree)
- 线段树 (Segment Tree)
- 可持续化线段树(Persistable Segment Tree)
- 树状数组套线段树求动态区间第 k 小
- 莫队（线段树一样都是通用的类型，具体问题具体写）
- 最长（严格）递增子序列
- 单调队列
- 单调栈
- 笛卡尔树
- 三维偏序之陈丹琪分治
- 第二分块差值版（在线算法，离线可节省空间）
- 第二分块绝对值版（在线算法，可搞成带修改版本！）