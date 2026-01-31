# 代码结构与实现文档

## 1. 概述

本文档说明项目的代码结构和实现细节。代码采用C++模板实现，支持在有限域 $\mathbb{Z}_p$ 中进行精确计算。主要功能包括FBI约化、微分方程求解和二维级数运算。

完整的数学理论背景请参考 [`problem_and_workflow.md`](./problem_and_workflow.md) 和论文 [2412.21053v1]。

## 2. 整体架构

### 2.1 文件组织

项目采用头文件和实现文件分离的结构：

- **`include/`**：头文件目录，包含所有类的声明
- **`src/`**：实现文件目录（`.tpp` 文件），包含模板类的实现
- **`test/`**：测试和应用程序
- **`docs/`**：文档目录

### 2.2 核心模块

代码分为以下几个核心模块：

1. **代数基础模块**：
   - `ffType.hpp`：有限域类型定义
   - `rational.hpp`：多项式和有理函数类
   - `series.hpp`：二维幂级数类

2. **FBI结构模块**：
   - `family.hpp`：FBI族（Family）管理
   - `sector.hpp`：Sector分析与分类
   - `fbi_reducer.hpp`：FBI数值约化引擎
   - `fbi_interpolater.hpp`：FBI符号约化（Firefly插值）

3. **微分方程模块**：
   - `de_builder.hpp`：数值点微分方程构建器
   - `de_interpolater.hpp`：符号微分方程构建器（Firefly插值）
   - `diffeq.hpp`：二维微分方程系统求解器

4. **辅助模块**：
   - `parser.hpp`、`lexer.hpp`、`token.hpp`：表达式解析
   - `linear.hpp`：线性代数工具
   - `integratedseries.hpp`：级数积分和输出

## 3. 核心类详解

### 3.1 有限域类型 (`ffType.hpp`)

#### `FlintMod` 类

基于FLINT库实现的有限域 $\mathbb{Z}_p$ 类型。

**核心功能：**
- 质数模运算（加、减、乘、除、幂）
- 模逆运算
- 与整数、字符串的转换

**数学背景：** 在有限域中进行计算可以避免浮点数误差，保证符号计算的精确性。通过中国剩余定理，可以从多个质数下的结果重建有理数解。

**关键方法：**
```cpp
static void set_modulus(int64_t p);  // 设置质数模
FlintMod inverse() const;             // 计算模逆
```

### 3.2 多项式和有理函数 (`rational.hpp`)

#### `Power` 结构体

表示单项式 $X^i Y^j$ 的幂次。

**成员变量：**
- `x_power`：$X$ 的指数
- `y_power`：$Y$ 的指数

#### `Polynomial<T>` 类

表示二变量多项式 $P(X, Y) = \sum_{i,j} c_{ij} X^i Y^j$。

**存储结构：** 使用哈希表 `unordered_map<Power, T>` 存储非零项，实现稀疏存储。

**核心操作：**
- `addMonomial(coeff, power)`：添加单项式
- `getCoeff(x, y)`：获取系数
- `derivativeX()`, `derivativeY()`：偏导数
- `evaluate(X, Y)`：多项式求值

**数学对应：** 对应数学框架中的多项式 $P(X,Y)$，用于表示有理函数的分子和分母。

#### `Rational<T>` 结构体

表示有理函数 $R(X,Y) = P(X,Y) / Q(X,Y)$。

**成员变量：**
- `numerator`：分子多项式 $P$
- `denominator`：分母多项式 $Q$

**核心操作：**
- 有理函数加法、乘法
- 判断是否为常数或空

**数学对应：** 对应数学框架第5.1节，用于表示微分方程系数矩阵的元素和FBI约化系数。

### 3.3 二维级数类 (`series.hpp`)

#### `Series<T>` 类

表示二维幂级数 $f(X,Y) = \sum_{i+j \leq d} c_{ij} X^i Y^j$。

**存储结构：** 
- 一维数组存储系数，按总度数排序
- 存储顺序：$(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), \ldots$
- 索引映射：`getIndex(i, j) = (i+j)(i+j+1)/2 + j`

**核心操作：**
- `getCoeff(i, j)`, `setCoeff(i, j, coeff)`：系数访问
- `operator+(other)`：级数加法
- `operator*(poly)`：级数与多项式乘法
- `operator*(rational)`, `operator/(rational)`：级数与有理函数的乘除

**关键静态方法（优化版本）：**
```cpp
static void mulRat(Series& result, const Series& series, const Rational& rational);
static void divRat(Series& result, const Series& series, const Rational& rational);
```

这些静态方法使用预分配空间，避免临时对象创建，提高性能。

**数学对应：** 
- 对应数学框架第4.3节的级数解表示
- 第5.1节的级数与有理函数运算算法在 `src/series.tpp` 中实现

### 3.4 Sector类 (`sector.hpp`)

#### `Sector<T>` 类

表示一个sector，包含矩阵 $S$ 的分析结果。

**构造参数：**
- `S`：$(N+B) \times (N+B)$ 矩阵（参考数学框架第2.3节）
- `numProps`：传播子数量 $N$
- `numBranch`：分支数量 $B$

**核心成员变量：**
- `S_`：输入矩阵
- `reducedS_`：行简化阶梯形（RREF）
- `rowOperation_`：行变换矩阵（当 `dimNull_==0` 时等于 $S^{-1}$）
- `dimNull_`：零空间维数，对应 $\dim(\text{Null}(S))$
- `candz_`：解向量 $(C_1, \ldots, C_B, z_1, \ldots, z_N)$
- `C_`：$C = \sum_{b=1}^{B} C_b$
- `z0_`：方程右侧的 $z_0$

**核心方法：**
- `getCase()`：返回Case类型（0, 1, 2, 3）
  ```cpp
  if (dimNull_ == 0 && C_ != T(0)) return 0;      // Case 0
  else if (dimNull_ == 0 && C_ == T(0)) return 1; // Case 1
  else if (dimNull_ != 0 && C_ != T(0)) return 2; // Case 2
  else return 3;                                   // Case 3
  ```

**数学对应：**
- 第2.4节：Sector定义与分类
- 第3.2节：求解线性方程组 $S \cdot (C_1, \ldots, C_B, z_1, \ldots, z_N)^T = (z_0, \ldots, z_0, 0, \ldots, 0)^T$

**实现细节（`src/sector.tpp`）：**
1. `rowReduce()`：高斯消元法计算RREF和行变换矩阵
2. `solveCandZ()`：根据 `dimNull_` 求解 $(C_b, z_\alpha, z_0)$
   - 若 `dimNull_ == 0`：直接求解线性方程
   - 若 `dimNull_ > 0`：取 $z_0 = 0$，求零空间的一组基，选择使 $C$ 尽可能非零的解

### 3.5 Family类 (`family.hpp`)

#### `Family<T>` 类

管理FBI族的全局信息，包括所有sector的分析。

**构造参数：**
- `topS`：顶层矩阵 $R(\mathbf{X})$（对应论文中的矩阵 $R$）
- `numProps`：传播子数量 $N$
- `numBranch`：分支数量 $B$

**核心成员变量：**
- `topS_`：顶层 $R$ 矩阵
- `branchIndices_`：每个传播子对应的分支索引
- `sectors_`：所有有效sector的 `Sector<T>` 对象
- `cases_`：每个sector的Case类型（0/1/2/3/-1）
- `masterIdxs_`：所有MFBI对应的sector索引

**核心方法：**
- `getCase(nu)`：获取指数向量 $\vec{\nu}$ 对应的Case类型
- `getSector(nu)`：获取指数向量对应的Sector对象
- `getIndexOfMaster(nu)`：判断 $\vec{\nu}$ 是否为MFBI，若是则返回其索引
- `isMaster(nu)`：判断是否为MFBI

**辅助方法：**
- `nBranch(nu)`：计算 $\vec{\nu}$ 涉及的分支数
- `nProps(nu)`：计算 $\vec{\nu}$ 中非零元素个数
- `nuSum(nu)`：计算 $\sum_i \nu_i$

**数学对应：**
- 第2.2节：分支结构
- 第2.4节：Sector分类
- 第6节第1-2步：Family初始化和Sector分析

**实现细节（`src/family.tpp`）：**

1. `constructBranchIndices()`：根据 `topS_` 确定每个传播子的分支归属
   - 同一分支内的传播子对应 `topS_` 中相同的行/列结构

2. `findSectors()`：枚举并分析所有sector
   - 遍历 $2^N$ 个可能的sector
   - 对每个sector构造子矩阵 $S$（调用 `getSubS()`）
   - 创建 `Sector<T>` 对象并分析
   - 记录Case类型
   - 识别MFBI（Case 0且 $\vec{\nu}$ 为角积分）

### 3.6 FBIReducer类 (`fbi_reducer.hpp`)

#### `FBIReducer<T>` 类

FBI约化引擎，将任意FBI表示为MFBI的线性组合。

**设计特点：**
- 持有 `Family<T>` 对象（非指针），在构造时创建
- 构造时指定工作维度delta，Family使用此delta初始化所有MFBI
- 线程安全：每个FBIReducer实例独立，可安全并行使用

**构造参数：**
- `topS`：数值S矩阵（$(N+B) \times (N+B)$）
- `numProps`：传播子数量 $N$
- `numBranch`：分支数量 $B$
- `delta`：工作维度（MFBI的目标delta值）

**核心方法：**
- `getReductionCoeff(nu, delta)`：获取约化系数向量
  - 参数：要约化的FBI指数 $\vec{\nu}$ 和维度 $\Delta$
  - 返回值：长度为 $M$ 的向量（$M$ 是MFBI数量）
  - 含义：$I_{\vec{\nu}}^{\Delta} = \sum_{k=1}^{M} c_k \cdot I_{\text{MFBI}_k}^{\Delta_0}$
- `getNumMaster()`：返回MFBI数量
- `getMasterNus()`：返回所有MFBI的指数向量列表
- `getMasterDeltas()`：返回所有MFBI的delta值（均为构造时的targetDelta）

**约化实现（`src/fbi_reducer.tpp`）：**

核心函数 `reduceFBI(nu, delta)` 根据Case类型分派：

##### Case 0：`reduceCase0(nu, delta)`

实现数学框架第3.3.1节的约化策略。

**步骤：**

1. **判断是否为MFBI**：
   ```cpp
   int mfbiIndex = family_->getIndexOfMaster(nu);
   if (mfbiIndex != -1) {
       // 返回单位向量
       vector<T> result(numMFBIs, T(0));
       result[mfbiIndex] = T(1);
       return result;
   }
   ```

2. **判断是否为角积分**（所有 $\nu_i = 1$）

3. **维度迁移**（如果不在目标维度）：
   - 使用公式（数学框架第3.3.1节）：
     $$I_{\vec{\nu}}^{\Delta} = \frac{C}{(2\Delta - \nu - B) z_0} I_{\vec{\nu}}^{\Delta-1} - \frac{1}{(2\Delta - \nu - B) z_0} \sum_{\alpha} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}$$
   
4. **IBP约化**（`case0IBP`）：
   - 构造右侧向量包含 $I_{\vec{\nu}-\vec{e}_{\max}}^{\Delta}$ 和 $I_{\vec{\nu}-\vec{e}_{\max}-\vec{e}_i}^{\Delta-1}$
   - 使用 $S^{-1}$ 求解：
     $$\nu_{\max} I_{\vec{\nu}}^{\Delta} = \sum_{j} (S^{-1})_{\text{row}, j} \cdot (\text{rhs})_j$$

**对应公式：** 数学框架第3.2节，方程(13)的递推关系。

##### Case 1：`reduceCase1(nu, delta)`

实现数学框架第3.3.2节。

**公式：**
$$(2\Delta - \nu - B) I_{\vec{\nu}}^{\Delta} = -\sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}$$

**代码实现：**
```cpp
T factor = T(2) * delta - T(nuSum) - T(numBranch);
// 计算右侧和，递归调用 getReductionCoeff(nu-e_alpha, delta-1)
// 结果除以 factor
```

##### Case 2：`reduceCase2(nu, delta)`

实现数学框架第3.3.3节。

**公式：**
$$I_{\vec{\nu}}^{\Delta} = \frac{1}{C} \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta}$$

注意这里维度不变，只降低指数。

##### Case 3：`reduceCase3(nu, delta)`

实现数学框架第3.3.4节。

**公式：**
$$I_{\vec{\nu}}^{\Delta} = -\sum_{\alpha \neq \beta} \frac{z_\alpha}{z_\beta} I_{\vec{\nu}+\vec{e}_\beta-\vec{e}_\alpha}^{\Delta}$$

**代码实现：**
1. 找到第一个非零的 $z_\beta$
2. 对所有 $\alpha \neq \beta$，递归约化 $I_{\vec{\nu}+\vec{e}_\beta-\vec{e}_\alpha}^{\Delta}$

**缓存机制：**
- 使用 `map<CacheKey, vector<T>> cache_` 缓存约化结果
- `CacheKey = tuple<vector<int>, T>`，对应 $(\vec{\nu}, \Delta)$
- 避免重复计算，显著提升性能

### 3.7 FBIInterpolater类 (`fbi_interpolater.hpp`)

#### `FBIInterpolater<T>` 类

使用Firefly插值重构FBI约化系数的符号形式（有理函数）。

**设计特点：**
- 将FBI约化系数从数值形式提升为符号有理函数 $r_k^{(\vec{\nu},\Delta)}(X,Y)$
- 使用Firefly黑盒插值方法在多个数值点求值并重构
- 线程安全：在并行插值中为每个线程创建独立的FBIReducer

**构造参数：**
- `polyTopS`：符号S矩阵（Polynomial<T>的矩阵，关于X,Y）
- `numProps`：传播子数量 $N$
- `numBranch`：分支数量 $B$
- `delta`：工作维度
- `prime`：有限域质数

**核心方法：**
- `getReductionCoeff(fbi_list)`：批量计算FBI约化系数
  - 参数：FBI列表 `vector<pair<vector<int>, T>>`（指数向量和delta）
  - 返回：每个FBI的约化系数向量 `vector<vector<Rational<T>>>`
  - 每个系数是关于$(X,Y)$的有理函数

**实现细节（`src/fbi_interpolater.tpp`）：**

1. **Firefly黑盒类** `FBIBlackBox<T>`：
   - 在数值点$(X,Y)$处：创建FBIReducer，调用数值约化，返回系数
   - 实现`operator()`接口供Firefly调用

2. **插值流程**：
   ```cpp
   // 创建Firefly Reconstructor（2个变量：X和Y）
   Reconstructor<FBIBlackBox<T>> reconst(2, n_threads, blackbox);
   
   // 设置质数并执行重构
   FFInt::set_new_prime(prime_);
   reconst.reconstruct(1);  // 单质数重构
   
   // 获取有理函数结果并转换
   auto ff_results = reconst.get_result_ff();
   ```

3. **线程安全设计**：
   - 黑盒在每次调用时创建新的FBIReducer（传入当前$(X,Y)$）
   - 避免共享状态，实现并行安全

**数学对应：** 数学框架第3节FBI约化，将数值约化提升为符号约化。

### 3.8 DEBuilder类 (`de_builder.hpp`)

#### `DEBuilder<T>` 类

在给定数值点$(X,Y)$构建MFBI的微分方程系数矩阵。

**设计特点：**
- 接受数值S矩阵和导数矩阵
- 在单个数值点计算微分方程系数
- 用于Firefly插值的黑盒求值

**构造参数：**
- `numTopS`：数值S矩阵（$(N+B) \times (N+B)$，在某个$(X,Y)$点求值）
- `numDRdX`：$\partial R/\partial X$的数值（$N \times N$矩阵）
- `numDRdY`：$\partial R/\partial Y$的数值（$N \times N$矩阵）
- `numProps`：传播子数量 $N$
- `numBranch`：分支数量 $B$
- `delta`：工作维度

**核心方法：**
- `buildDEMatrices(AX, AY)`：构建数值微分方程矩阵
  - `AX[i][j]`, `AY[i][j]`：数值系数（类型T）
  - 表示：$\frac{\partial f_i}{\partial X} = \sum_j A_X[i][j] \cdot f_j$

**实现细节（`src/de_builder.tpp`）：**

1. **FBI导数公式**：
   $$\frac{\partial I_{\vec{\nu}}^{\Delta}}{\partial X} = \sum_{i,j} -\frac{1}{2} \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}$$
   
   其中 $\text{factor}_{ij} = \begin{cases} \nu_i \nu_j & i \neq j \\ \nu_i(\nu_i+1) & i = j \end{cases}$

2. **构建流程**：
   - 对每个MFBI计算其对$X$和$Y$的导数（得到FBI的线性组合）
   - 使用FBIReducer将每个FBI约化到MFBI基
   - 组合系数得到矩阵元素$A_X[i][j]$和$A_Y[i][j]$

**数学对应：** 数学框架第4.4节，构建微分方程系数矩阵。

### 3.9 DEInterpolater类 (`de_interpolater.hpp`)

#### `DEInterpolater<T>` 类

使用Firefly插值重构MFBI微分方程系数矩阵的符号形式。

**设计特点：**
- 将微分方程系数从数值形式提升为符号有理函数矩阵
- 使用Firefly黑盒插值：在多个$(X,Y)$点调用DEBuilder求值
- 线程安全：为每个求值点创建独立的DEBuilder

**构造参数：**
- `topS`：符号S矩阵（Polynomial<T>的矩阵）
- `numProps`：传播子数量 $N$
- `numBranch`：分支数量 $B$
- `delta`：工作维度
- `prime`：有限域质数

**核心成员变量：**
- `topS_`：符号S矩阵
- `dRdX_`, `dRdY_`：符号导数矩阵（在构造时计算）
- `numMasterFBI_`：MFBI数量

**核心方法：**
- `buildDEMatrices(AX, AY)`：构建符号微分方程矩阵
  - 返回：$M \times M$有理函数矩阵，$A_X[i][j], A_Y[i][j]$类型为`Rational<T>`

**辅助方法：**
- `evaluateTopS(X,Y)`, `evaluateDRdX(X,Y)`, `evaluateDRdY(X,Y)`：在数值点求值
- `createDEBuilder(X,Y)`：创建DEBuilder对象用于黑盒求值
- `convertToPolynomial()`, `convertToRational()`：Firefly结果转换

**实现细节（`src/de_interpolater.tpp`）：**

1. **Firefly黑盒类** `DEBlackBox<T>`：
   - 在数值点$(X,Y)$：调用`createDEBuilder`，计算数值矩阵，返回扁平化结果
   - 返回向量顺序：先$A_X$的所有元素，后$A_Y$的所有元素

2. **插值流程**：
   ```cpp
   // 创建Reconstructor
   Reconstructor<DEBlackBox<T>> reconst(2, n_threads, blackbox);
   
   // 设置质数并重构
   FFInt::set_new_prime(prime_);
   reconst.reconstruct(1);
   
   // 获取结果并解析为AX和AY矩阵
   auto ff_results = reconst.get_result_ff();
   // 前M*M个结果 -> AX，后M*M个结果 -> AY
   ```

3. **构造函数行为**：
   - 计算符号导数：`dRdX_[i][j] = topS_[i+B][j+B].derivativeX()`
   - 创建临时FBIReducer确定MFBI数量（evaluateTopS(1,1)）
   - **重要**：此步骤可能在Firefly插值后执行，需确保质数正确

**数学对应：** 数学框架第4.4节，构建微分方程系数矩阵 $A_X$ 和 $A_Y$。

### 3.10 微分方程求解类 (`diffeq.hpp`)

#### `DiffSystem<T>` 类

求解二维微分方程系统。

**构造参数：**
- `AX`：$\partial/\partial X$ 方向的系数矩阵 $A_X$
- `AY`：$\partial/\partial Y$ 方向的系数矩阵 $A_Y$

矩阵元素类型为 `Rational<T>`。

**核心方法：**
```cpp
void solve(vector<Series<T>>& result, const vector<T>& f0, int deg);
```

求解微分方程系统：
$$\frac{\partial \mathbf{f}}{\partial X} = A_X \mathbf{f}, \quad \frac{\partial \mathbf{f}}{\partial Y} = A_Y \mathbf{f}$$

给定初值 $\mathbf{f}(0,0) = \mathbf{f}_0$，计算级数解到度数 `deg`。

**数学对应：** 数学框架第4节，微分方程求解。

**实现细节（`src/diffeq.tpp`）：**

核心算法在 `solveStandardDE(result, R1, g1, R2, g2, f0)`：

1. **预处理**：
   - 设 $R_1 = P_1/Q_1$，$R_2 = P_2/Q_2$
   - 验证 $Q_1(0,0) \neq 0$，$Q_2(0,0) \neq 0$

2. **初值设置**：
   ```cpp
   result.setCoeff(0, 0, f0);
   ```

3. **递推计算**（按总度数 $m+n$ 递增）：

   **当 $m > 0$ 时**，使用 $X$ 方向方程（对应数学框架第4.4节公式）：
   $$c_{mn} = \frac{1}{m \cdot Q_1^{(0,0)}} \left[ [Q_1 g_1]_{m-1,n} + \sum P_1^{(k,l)} c_{i,j} - \sum (i+1) Q_1^{(k,l)} c_{i+1,j} \right]$$
   
   其中求和指标满足 $k = m-1-i$，$l = n-j$。

   **当 $m = 0, n > 0$ 时**，使用 $Y$ 方向方程：
   $$c_{0n} = \frac{1}{n \cdot Q_2^{(0,0)}} \left[ [Q_2 g_2]_{0,n-1} + \sum P_2^{(k,l)} c_{i,j} - \sum (j+1) Q_2^{(k,l)} c_{i,j+1} \right]$$

4. **性能优化**：
   - 代码中根据循环范围和多项式单项式数选择迭代策略
   - 若多项式稀疏，直接遍历单项式；否则遍历 $(i,j)$ 范围

**耦合项计算（`getSub`）：**

对于第 $i$ 个方程，计算耦合项：
$$g_1 = Q_1^{(ii)} \sum_{j<i} A_X^{(ij)} f_j, \quad g_2 = Q_2^{(ii)} \sum_{j<i} A_Y^{(ij)} f_j$$

这确保了三角化求解：第 $i$ 个方程只依赖前 $i-1$ 个已求解的函数。

## 4. 计算流程

### 4.1 主程序示例（`test/calcseries.cpp`）

完整的计算流程：

```cpp
// 1. 设置有限域模数
FlintMod::set_modulus(prime);

// 2. 加载数据
auto AX = loadMatrix<FlintMod>("DEX", parser);  // X方向微分方程矩阵
auto AY = loadMatrix<FlintMod>("DEY", parser);  // Y方向微分方程矩阵
auto red_coe = loadMatrix<FlintMod>("coe", parser);  // 约化系数矩阵

// 3. 创建微分方程系统并求解
DiffSystem<FlintMod> diffSys(AX, AY);
vector<Series<FlintMod>> mseries = solveDiffEquations(diffSys, degree);

// 4. 计算所有FBI的级数展开
auto series = computeSeries(red_coe, mseries, degree);

// 5. 积分和输出
auto integrated_series = performIntegration(series, powers, degree);
outputIntegratedSeries(integrated_series, config, ...);
```

**说明：**
- `AX`、`AY`：MFBI的微分方程系数矩阵
- `red_coe`：FBI到MFBI的约化系数矩阵
  - 行索引：各个FBI
  - 列索引：MFBI
  - 元素：有理函数 $r_k^{(\vec{\nu}, \Delta)}(X, Y)$

### 4.2 数据文件格式

**微分方程矩阵（DEX, DEY）：**
```
行数 列数
R(0,0) R(0,1) ...
R(1,0) R(1,1) ...
...
```
每个元素是有理函数的字符串表示，例如：`"(X+Y)/(1-X*Y)"`

**约化系数矩阵（coe）：**
```
FBI数量 MFBI数量
r(0,0) r(0,1) ...
r(1,0) r(1,1) ...
...
```

### 4.3 约化系数的计算

约化系数矩阵通常由单独的程序生成（不在当前代码库中）：

1. 初始化 `Family<T>` 对象
2. 枚举需要计算的所有FBI $(\vec{\nu}_k, \Delta_k)$
3. 对每个FBI调用 `FBIReducer::getReductionCoeff(nu_k, delta_k)`
4. 将约化系数写入文件

## 5. 性能优化策略

### 5.1 缓存机制

**FBIReducer缓存：**
- 约化过程高度递归，相同 $(\vec{\nu}, \Delta)$ 可能多次出现
- 使用 `map` 缓存避免重复计算

**Series预分配：**
- 静态方法 `mulRat`、`divRat` 等使用预分配的结果对象
- 避免临时对象创建和拷贝

### 5.2 稀疏存储

**多项式存储：**
- 使用 `unordered_map` 仅存储非零单项式
- 哈希函数优化（`PowerHash`）

**级数存储：**
- 密集存储（数组），因为通常大部分系数非零
- 按总度数排序，便于递推计算

### 5.3 算法选择

**微分方程求解：**
- 根据多项式稀疏度动态选择遍历策略
- 稀疏多项式：遍历单项式
- 稠密求和范围：遍历 $(i,j)$ 网格

## 6. 类型参数化

所有核心类都是模板类 `template<typename T>`，支持不同的系数类型：

- **`FlintMod`**：有限域 $\mathbb{Z}_p$（主要用途）
- **`int`、`long`**：整数运算（测试用）
- **`mpq_class`**：有理数（高精度）

模板设计使得代码可以灵活适配不同的数值类型。

## 7. 辅助模块

### 7.1 表达式解析（`parser.hpp`, `lexer.hpp`, `token.hpp`）

**功能：** 将字符串表示的有理函数解析为 `Rational<T>` 对象。

**示例：**
```cpp
Parser<FlintMod> parser;
auto rational = parser.parseRational("(X+2*Y)/(1-X*Y)");
```

**词法分析（Lexer）：**
- 识别标记：数字、变量（X, Y）、运算符（+, -, *, /, ^）、括号

**语法分析（Parser）：**
- 递归下降解析
- 构造多项式和有理函数对象

### 7.2 线性代数（`linear.hpp`）

提供矩阵运算支持：
- 矩阵乘法
- 高斯消元
- 行简化阶梯形

被 `Sector` 类用于矩阵 $S$ 的分析。

### 7.3 积分级数（`integratedseries.hpp`）

**功能：** 对级数进行数值积分和格式化输出。

在二圈情形下，需要对 $X, Y$ 参数进行积分：
$$\int_0^1 \int_0^1 f(X, Y) \, dX \, dY$$

使用预计算的幂次积分：
$$\int_0^1 \int_0^1 X^m Y^n \, dX \, dY = \frac{1}{(m+1)(n+1)}$$

## 8. 数学公式与代码对应总结

| 数学概念 | 数学框架章节 | 代码类/函数 | 关键方法 |
|---------|------------|-----------|---------|
| FBI定义 | 2.1 | - | - |
| 矩阵 $S$ | 2.3 | `Sector` | 构造函数 |
| Sector分类 | 2.4 | `Sector::getCase()` | 返回0/1/2/3 |
| 约化关系（IBP） | 3.2, eq.(13) | `FBIReducer::case0IBP()` | 使用 $S^{-1}$ |
| 维度迁移 | 3.2, eq.(14) | `FBIReducer::reduceCase0()` | 维度调整部分 |
| Case 0约化 | 3.3.1 | `FBIReducer::reduceCase0()` | - |
| Case 1约化 | 3.3.2 | `FBIReducer::reduceCase1()` | - |
| Case 2约化 | 3.3.3 | `FBIReducer::reduceCase2()` | - |
| Case 3约化 | 3.3.4 | `FBIReducer::reduceCase3()` | - |
| 微分方程系统 | 4.1 | `DiffSystem` | - |
| 标准方程 | 4.2 | `DiffSystem::solveStandardDE()` | - |
| 系数递推 | 4.4 | `DiffSystem::solveStandardDE()` | 主循环 |
| 级数×有理函数 | 5.1 | `Series::mulRat()`, `divRat()` | - |
| FBI级数展开 | 5.2 | `computeSeries()` (test/) | - |

## 9. 参考文献

[2412.21053v1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
