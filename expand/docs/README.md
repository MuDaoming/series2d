# 项目说明

## 1. 项目目标

计算二圈Feynman积分（FBI）的二维幂级数展开，在质数域 $\mathbb{Z}_p$ 下进行数值计算。

## 2. 新旧项目关系

### 2.1 目录结构

```
series2d_new/
├── backup/              # 旧版代码（参考用）
│   ├── docs/            # 旧版文档
│   ├── include/         # 旧版头文件
│   └── src/             # 旧版实现
├── docs/                # 新版文档
├── include/             # 新版头文件
├── src/                 # 新版实现
└── test/                # 测试代码
```

### 2.2 新旧方法对比

| 方面 | 旧方法 (backup) | 新方法 |
|-----|----------------|--------|
| **约化结果** | 有理函数系数 | 级数系数 |
| **核心步骤** | 先约化成有理函数，后乘级数 | 直接LRR递推级数 |
| **数据流** | FBI → 有理函数向量 → 级数 | FBI → 级数 |
| **效率** | 数值插值重构有理函数开销大 | 直接数值递推更高效 |

### 2.3 类对应关系

| 旧类 (backup) | 新类 | 变化 |
|--------------|------|------|
| `Family<T>` | `Family<T>` | **已实现**，逻辑基本相同 |
| `Sector<T>` | `Sector<T>` | **已实现**，逻辑基本相同 |
| `FBIReducer<T>` | `SeriesSolver<T>` | 待实现，约化直接输出级数 |
| `DEBuilder<T>` | 合并到 `SeriesSolver<T>` | 微分方程和约化统一处理 |
| `Polynomial<T>` | `Polynomial<T>` | 待实现，可参考旧版 |
| `Series<T>` | `Series<T>` | 待实现，可复用旧版 |
| `Rational<T>` | 不需要 | 新方法不使用有理函数 |

## 3. 参考资源

### 3.1 数学原理

- **`backup/docs/problem_and_workflow.md`**：完整的数学推导和计算流程
  - FBI定义（方程11）
  - 矩阵S的构造（方程12）
  - IBP递推关系（方程13）
  - 维度迁移关系（方程14）
  - 四种Case的约化策略
  - 微分方程的构建

### 3.2 约化逻辑

- **`backup/include/fbi_reducer.hpp`** + **`backup/src/fbi_reducer.tpp`**
  - `reduceCase0()` - Case 0约化（IBP + 维度迁移）
  - `reduceCase1()` - Case 1约化
  - `reduceCase2()` - Case 2约化
  - `reduceCase3()` - Case 3约化
  - `isCorner()` - 判断是否为角积分
  - `findMaxIndex()` - 找最大指数

### 3.3 微分方程构建

- **`backup/include/de_builder.hpp`** + **`backup/src/de_builder.tpp`**
  - `computeFBIDerivative()` - 计算FBI对X/Y的导数
  - `buildDEMatrices()` - 构建微分方程矩阵A_X和A_Y

### 3.4 级数运算

- **`backup/include/series.hpp`** + **`backup/src/series.tpp`**
  - 二维幂级数的存储方式
  - `mulPoly()` - 级数与多项式乘法
  - `divPoly()` - 级数除以多项式
  - `getIndex(i,j)` - 二维坐标转一维索引

### 3.5 多项式/有理函数

- **`backup/include/rational.hpp`**：有理函数类（新方法不需要）
- **多项式类**：backup中嵌入在rational中，需要独立实现

## 4. 新版实现状态

### 4.1 已完成

- [x] `Family<T>` - 管理所有sector
- [x] `Sector<T>` - 处理单个sector的S矩阵、计算C和z
- [x] `isZero<T>()` / `normalize<T>()` - 类型萃取辅助函数

### 4.2 待实现

- [ ] `Polynomial<T>` - 二变量多项式类
- [ ] `Series<T>` - 二维幂级数类
- [ ] `SeriesSolver<T>` - 级数求解器（核心）

## 5. 关键公式速查

### 5.1 LRR递推

$$g_{pq} = \frac{1}{D_{00}} \left( \sum_i [N_i \cdot f_i]_{pq} - \sum_{(a,b)\neq(0,0)} D_{ab} \cdot g_{p-a,q-b} \right)$$

### 5.2 微分方程

$$\frac{\partial I_{\vec{\nu}}^{\Delta}}{\partial X} = \sum_{i,j} \left(-\frac{1}{2}\right) \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}$$

### 5.3 维度迁移（Case 0）

$$C \cdot I_{\vec{\nu}}^{\Delta-1} = (2\Delta - \nu - B) z_0 \cdot I_{\vec{\nu}}^{\Delta} + \sum_{\alpha} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}$$

## 6. 文档索引

| 文档 | 内容 |
|-----|------|
| `docs/README.md` | 项目说明（本文档） |
| `docs/code_structure.md` | 代码架构和类设计 |
| `docs/problem_and_workflow.md` | 问题定义和计算流程 |
| `docs/series_solver.md` | SeriesSolver设计文档 |
| `backup/docs/problem_and_workflow.md` | 完整的数学原理（参考） |
