# FI 关于 $\delta$ 的多项式关系搜索：问题定义与工作流程


## 1. 研究目标

本模块的目标是：
- 输入若干 Feynman 积分 $\mathrm{FI}_{\vec{\nu}}$ 的一维级数展开
- 在有限域 $Z_p$ 下，搜索它们之间关于 $\delta$ 的多项式线性关系

目标关系写成：

$$
\sum_{\vec{\nu} \in G}
\left(
c^{\vec{\nu}}_0 + c^{\vec{\nu}}_1 \delta + \cdots + c^{\vec{\nu}}_{m} \delta^{m}
\right)
\mathrm{FI}_{\vec{\nu}}(\delta)
= 0
$$

其中：

- $G$ 是选定的一组积分指标 $\vec{\nu}$
- $m$ 是允许的 $\delta$ 多项式最高次数
- $c_k^{\vec{\nu}} \in Z_p$ 是待求系数
- $\mathrm{FI}_{\vec{\nu}}(\delta)$ 已知到某个一维展开阶数 $d$

通常取

$$
d > m
$$

也就是说，已知的级数阶数高于所假设关系中 $\delta$ 多项式的最高次数。

本模块不负责计算 $\mathrm{FI}_{\vec{\nu}}$ 本身，而是负责从已知的一维级数中搜索上式的非平凡解。

## 2. 输入与输出

### 2.1 输入

搜索问题的输入由四部分组成：

1. 一组目标积分 $G$
   每个元素是一个向量 $\vec{\nu}$。

2. 一组独立边界条件  
   如果主积分个数为 $N_{\mathrm{master}}$，则输入 $N_{\mathrm{master}}$ 组独立边界条件。每组边界条件可取标准基形式

   $$
   (1,0,\dots,0),\quad
   (0,1,\dots,0),\quad
   \dots,\quad
   (0,0,\dots,1)
   $$

   对每一组边界条件，都可以得到同一批 $\mathrm{FI}_{\vec{\nu}}$ 的一维展开。

3. 每个 $\mathrm{FI}_{\vec{\nu}}$ 在每组边界条件下的一维展开  
   记第 $b$ 组边界条件下的展开为：

   $$
   \mathrm{FI}_{\vec{\nu}}^{(b)}(\delta)
   =
   \sum_{r=0}^{d} a^{\vec{\nu},(b)}_r \delta^r
   $$

   其中系数 $a_r^{\vec{\nu},(b)}$ 属于有限域 $Z_p$。

4. 两个截断参数
   包括展开已知到的最高阶 $d$，以及关系式中允许的 $\delta$ 多项式最高次数 $m$。

### 2.2 输出

第一步直接得到的是一个或多个非平凡关系：

$$
\sum_{\vec{\nu} \in G}
\left(
c^{\vec{\nu}}_0 + c^{\vec{\nu}}_1 \delta + \cdots + c^{\vec{\nu}}_{m} \delta^{m}
\right)
\mathrm{FI}_{\vec{\nu}}(\delta)
= 0
$$

更具体地，第一步输出可以采用下列两种等价形式之一：

1. 输出行约化后的线性关系结果  
   例如把某些未知系数写成其它未知系数的线性组合：

   $$
   c_9 = 111 c_1 + 12131 c_2
   $$

   这种形式直接反映零空间中的自由变量与主变量关系。

2. 输出秩亏线性系统的行最简形  
   即给出矩阵消元后的结果，由此可直接读出全部关系。

但这还不是最终目标。

最终目标是得到原始 Feynman 积分之间的约化关系，也就是把更复杂的 $\mathrm{FI}_{\vec{\nu}}$ 写成更简单的 $\mathrm{FI}_{\vec{\mu}}$ 的线性组合。

因此，在第一步得到 $c_k^{\vec{\nu}}$ 的关系之后，还需要再做一次整理：

1. 先取第一步得到的零空间基；
2. 对每一个自由变量，依次取 1，其余自由变量取 0；
3. 这样得到一组具体的 $c_k^{\vec{\nu}}$；
4. 在每一组具体关系中令 $\delta = 1$；
5. 于是得到一批关于原始 $\mathrm{FI}_{\vec{\nu}}$ 的线性关系；
6. 再对这些 $\mathrm{FI}_{\vec{\nu}}$ 关系做一次消元，最终整理成约化形式。

因此，最终输出应当是形如

$$
\mathrm{FI}_{\vec{\nu}_{\mathrm{comp}}}
=
\sum_{\vec{\mu}\ \text{更简单}}
r_{\vec{\nu}_{\mathrm{comp}},\vec{\mu}}
\mathrm{FI}_{\vec{\mu}}
$$

的约化关系。

## 3. 数学结构

### 3.1 未知量

对每个 $\vec{\nu} \in G$ 和每个 $0 \le k \le m$，引入一个未知量：

$$
x_{\vec{\nu},k} = c_k^{\vec{\nu}}
$$

因此未知量总数为：

$$
N_{\mathrm{var}} = |G| \cdot (m + 1)
$$

为使最终输出更接近“复杂量由简单量表示”的形式，需要对这些未知量固定一个复杂度顺序。

对每个 $\vec{\nu}$，定义：

$$
\mathrm{props}(\vec{\nu})=\#\{i\mid \nu_i\neq 0\}
$$

以及

$$
\mathrm{dots}(\vec{\nu})=\sum_i \nu_i-\mathrm{props}(\vec{\nu})
$$

其中：

- $\mathrm{props}(\vec{\nu})$ 表示非零分量数目
- $\mathrm{dots}(\vec{\nu})$ 表示总幂次超过 1 的部分

规定“更简单”的变量满足：

1. $\mathrm{props}(\vec{\nu})$ 更小；
2. 若 $\mathrm{props}(\vec{\nu})$ 相同，则 $\mathrm{dots}(\vec{\nu})$ 更小；
3. 若上述两者都相同，则 $k$ 更小。

因此，为了让消元结果优先把复杂变量写成简单变量的线性组合，变量列顺序应按“从复杂到简单”排列，也就是：

1. $\mathrm{props}(\vec{\nu})$ 大的在前；
2. 若 $\mathrm{props}(\vec{\nu})$ 相同，则 $\mathrm{dots}(\vec{\nu})$ 大的在前；
3. 若上述两者都相同，则 $k$ 大的在前；
4. 若仍然相同，再按任意固定的规则比较 $\vec{\nu}$ 本身即可。

### 3.2 已知量

对每一组边界条件 $b$，每个积分的一维展开已知：

$$
\mathrm{FI}_{\vec{\nu}}^{(b)}(\delta)
=
\sum_{r=0}^{d} a_r^{\vec{\nu},(b)} \delta^r
$$

于是：

$$
\left(
\sum_{k=0}^{m} c_k^{\vec{\nu}} \delta^k
\right)
\mathrm{FI}_{\vec{\nu}}^{(b)}(\delta)
=
\sum_{n=0}^{d}
\left(
\sum_{k=0}^{\min(n,m)} c_k^{\vec{\nu}} a_{n-k}^{\vec{\nu},(b)}
\right)\delta^n
$$

这里只保留到 $\delta^d$，更高阶项忽略。

### 3.3 线性方程组

同一个拟设关系必须对每一组边界条件同时成立。因此，对每个边界条件 $b$，都要求：

$$
\sum_{\vec{\nu} \in G}
\sum_{k=0}^{\min(n,m)}
c_k^{\vec{\nu}} a_{n-k}^{\vec{\nu},(b)} = 0,
\quad n = 0,1,\dots,d
$$

也就是说，每组边界条件都会贡献 $d+1$ 条线性方程。把所有边界条件的方程堆叠起来，就得到一个总的齐次线性系统：

$$
A x = 0
$$

其中：

- 行数：$(d+1)N_{\mathrm{master}}$
- 列数：$|G|(m+1)$
- $A$ 的每一列对应一个未知量 $c_k^{\vec{\nu}}$

若把矩阵的行指标写成 $(b,n)$，列指标写成 $(\vec{\nu},k)$，则矩阵元为：

$$
A_{(b,n),(\vec{\nu},k)} =
\begin{cases}
a_{n-k}^{\vec{\nu},(b)}, & n \ge k \\
0, & n < k
\end{cases}
$$

## 4. 处理流程

给定：

- 一组 $\mathrm{FI}_{\vec{\nu}}$
- 一组独立边界条件
- 截断阶数 $d$
- 多项式次数上限 $m$

整个问题的处理顺序如下：

1. 固定待搜索的积分集合 $G$。

2. 固定关系中允许出现的 $\delta$ 多项式次数上限 $m$。

3. 对每个边界条件 $b$、每个 $\vec{\nu} \in G$，读入

   $$
   \mathrm{FI}_{\vec{\nu}}^{(b)}(\delta)
   =
   \sum_{r=0}^{d} a_r^{\vec{\nu},(b)} \delta^r
   $$

   的截断级数。

4. 以所有 $c_k^{\vec{\nu}}$ 为未知量，按卷积公式构造总矩阵 $A$。

5. 求解齐次线性系统 $Ax=0$ 的零空间。

6. 将零空间结果写成显式关系，或者写成行约化后的矩阵形式。

7. 对每一个自由变量，依次取 1，其余自由变量取 0，从而得到一组具体的 $c_k^{\vec{\nu}}$ 解。

8. 对每一组具体解，将关系中的 $\delta$ 取为 1，于是得到一批关于原始 $\mathrm{FI}_{\vec{\nu}}$ 的线性关系：

   $$
   \sum_{\vec{\nu}\in G}
   \left(
   \sum_{k=0}^{m} c_k^{\vec{\nu}}
   \right)
   \mathrm{FI}_{\vec{\nu}}
   =0
   $$

9. 将这些 $\mathrm{FI}_{\vec{\nu}}$ 关系重新组成线性系统，并按积分复杂度对变量排序。

10. 再做一次消元，把复杂积分写成简单积分的线性组合。

## 5. 与现有 linear 模块的关系

$\texttt{reconstruct/include/linear.hpp}$ 已经提供了线性系统求解功能。

本模块与它的关系是：

- search 模块负责把 $\mathrm{FI}_{\vec{\nu}}$ 的级数数据组织成矩阵 $A$
- linear 模块第一次用于对 $A$ 做高斯消元并给出零空间信息
- search 模块再把零空间基转成 $\mathrm{FI}_{\vec{\nu}}$ 的实际关系
- linear 模块第二次用于对这些 $\mathrm{FI}_{\vec{\nu}}$ 关系做消元，得到最终约化式

因此 search 的核心不是重新实现线性代数，而是：

1. 定义变量顺序
2. 正确构造由所有边界条件共同给出的总矩阵
3. 把解空间结果重新映射回 $c_k^{\vec{\nu}}$
4. 从 $c_k^{\vec{\nu}}$ 的零空间生成 $\mathrm{FI}_{\vec{\nu}}$ 的实际关系
5. 再把这些关系整理成最终约化式
