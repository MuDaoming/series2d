# 代码结构与实现文档

## 1. 概述

本文档说明新版项目的代码结构和实现细节。新版本采用**符号计算+数值计算**混合方法：
- **符号阶段**：使用GiNaC计算多项式/有理函数形式的C和z
- **数值阶段**：逐阶递推计算级数系数

完整的数学理论背景请参考 [`problem_and_workflow.md`](./problem_and_workflow.md) 和论文 [2412.21053v1]。

## 2. 已实现的核心类

### 2.1 Sector类 `Sector<T>`

**文件**：`include/sector.hpp`, `src/sector.tpp`

处理单个sector的S矩阵，进行RREF行化简，计算C和z系数。

**模板参数**：
- `T` - 数值类型，支持：
  - `GiNaC::ex`（符号计算）
  - 数值类型（`int`, `double`, 有限域类型等）

**核心功能**：
- `rowReduce()` - 行化简为RREF形式
- `solveCandZ()` - 求解C和z系数
- `getCase()` - 返回Case类型（0/1/2/3）

```cpp
template<typename T>
class Sector { 
public:
    // 构造函数：自动调用rowReduce()和solveCandZ()
    Sector(const std::vector<std::vector<T>>& S, int numProps, int numBranch);
    
    int getCase() const;                // Case类型 (0/1/2/3)
    const std::vector<std::vector<T>>& getInvS() const;  // S的逆（dimNull=0时）

private:
    // 核心算法
    void rowReduce();                   // RREF行化简
    void solveCandZ();                  // 求解C和z系数
    std::vector<std::vector<T>> findNullSpace() const;   // 零空间基
    std::vector<T> solveLinear(const std::vector<T>& b) const;  // 解线性方程

    // 成员变量
    std::vector<std::vector<T>> S_;           // 输入的(N+B)×(N+B)矩阵
    int numProps_;                            // 传播子数N
    int numBranch_;                           // 分支数B
    std::vector<std::vector<T>> reducedS_;    // RREF后的矩阵
    std::vector<std::vector<T>> rowOperation_;// 行变换矩阵（dimNull=0时为逆）
    int dimNull_;                             // 零空间维度
    T z0_;                                    // z_0
    std::vector<T> candz_;                    // (C_1,...,C_B,z_1,...,z_N)
    T C_;                                     // C = C_1 + ... + C_B
};
```

**Case分类**：
| Case | dimNull | C | 说明 |
|------|---------|---|------|
| 0 | =0 | ≠0 | 包含主积分 |
| 1 | =0 | =0 | 无主积分 |
| 2 | >0 | ≠0 | 无主积分 |
| 3 | >0 | =0 | 无主积分 |

**类型萃取与辅助函数**（同文件）：

使用`if constexpr`和类型萃取来区分GiNaC::ex和数值类型：

```cpp
// 类型萃取：判断是否为GiNaC::ex
template<typename T>
struct is_ginac_ex : std::is_same<T, GiNaC::ex> {};

// 零判断函数
template<typename T>
inline bool isZero(const T& val) {
    if constexpr (is_ginac_ex<T>::value) {
        return GiNaC::normal(val).is_zero();  // 符号类型
    } else {
        return val == T(0);  // 数值类型
    }
}

// 化简函数
template<typename T>
inline T normalize(const T& val) {
    if constexpr (is_ginac_ex<T>::value) {
        return GiNaC::normal(val);  // 约分有理函数
    } else {
        return val;  // 数值类型直接返回
    }
}
```

### 2.2 Family类 `Family<T>`

**文件**：`include/family.hpp`, `src/family.tpp`

管理FBI族的全局信息，枚举所有有效sector。

**核心功能**：
- `constructBranchIndices()` - 从topS构造分支索引
- `findSectors()` - 枚举所有有效sector
- `getSubS()` - 构造子矩阵

```cpp
template<typename T>
class Family {
public:
    // 构造函数：自动调用constructBranchIndices()和findSectors()
    Family(const std::vector<std::vector<T>>& topS, int numProps, int numBranch);
    
    // 获取Sector
    const Sector<T>* getSector(std::vector<int> nu) const;
    const Sector<T>* getSectorByIdx(int idx) const;
    
    // 索引转换
    std::vector<int> secvecFromIdx(int n) const;    // 索引 -> secvec
    int idxFromSecvec(const std::vector<int>& secvec) const;  // secvec -> 索引
    
    // 辅助函数
    bool isMaster(const std::vector<int>& nu) const;
    void setMasterDelta(T delta);        // 设置所有主积分的delta

private:
    // 核心算法
    void constructBranchIndices();       // 从topS构造分支索引
    void findSectors();                  // 枚举所有有效sector
    bool getSubS(const std::vector<int>& nu, std::vector<std::vector<T>>& subS) const;

    // 成员变量
    std::vector<std::vector<T>> topS_;   // 顶层S矩阵
    int numProps_;                       // 传播子数N
    int numBranch_;                      // 分支数B
    std::vector<int> branchIndices_;     // 每个传播子对应的分支索引
    std::vector<Sector<T>> sectors_;     // 所有有效sector
    std::vector<int> sectorIdxs_;        // sector索引列表
    std::vector<int> cases_;             // 每个sector的Case（按idx索引）
    std::vector<int> masterIdxs_;        // 主积分的sector索引列表
    std::vector<T> masterDeltas_;        // 主积分的delta参数（默认T(0)）
};
```

## 3. 文件组织

```
series2d_new/
├── docs/
│   ├── problem_and_workflow.md    # 问题定义与流程
│   └── code_structure.md          # 代码架构（本文档）
├── include/
│   ├── sector.hpp                 # Sector类 + isZero/normalize
│   └── family.hpp                 # Family类
└── src/
    ├── sector.tpp                 # Sector实现
    └── family.tpp                 # Family实现
```

## 4. 使用示例

### 4.1 符号计算（使用GiNaC::ex）

```cpp
#include <ginac/ginac.h>
#include "family.hpp"

using namespace GiNaC;

int main() {
    symbol X("X"), Y("Y");
    
    // 构造S矩阵（包含X,Y的多项式）
    std::vector<std::vector<ex>> topS = {
        {0, 0, 0, 1, 1, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 0},
        // ... 9x9矩阵
    };
    
    // 创建Family
    Family<ex> family(topS, 6, 3);  // numProps=6, numBranch=3
    
    // 遍历所有sector
    for (const auto& sector : family.getSectors()) {
        std::cout << "C = " << sector.getCSum() << "\n";
        for (int i = 0; i < sector.getNumProps(); ++i) {
            std::cout << "z[" << i << "] = " << sector.getZ(i) << "\n";
        }
    }
    
    return 0;
}
```

### 4.2 编译命令

```bash
g++ -std=c++17 -I./include -o program main.cpp -lginac -lcln
```

## 5. 待实现组件

### 5.1 第二阶段：级数与多项式

- [ ] `Polynomial<T>` - 二变量多项式类
- [ ] `Series<T>` - 二维幂级数类
- [ ] 输入文件解析器

### 5.2 第三阶段：逐阶递推

- [ ] `FBISystem<T>` - 逐阶递推求解器
- [ ] 微分方程递推
- [ ] IBP约化递推（四种Case）

### 5.3 第四阶段：集成

- [ ] Mathematica预计算脚本导出
- [ ] 性能优化

## 6. 参考文献

[2412.21053v1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
