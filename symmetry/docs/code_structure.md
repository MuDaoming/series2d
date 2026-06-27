# FBI Sector Symmetry：代码结构与实现文档

## 1. 概述

本模块实现 [`problem_solution.md`](./problem_solution.md) 定义的 FBI sector 参数置换 symmetry 搜索。

实现原则：

1. 直接读取 `expand` 当前的 $S(X,Y)$。
2. 显式枚举 branch 置换 $\tau$。
3. 对固定 $\tau$ 使用 Pak 型前缀剪枝规范化传播子排列 $\sigma$。
4. 按 canonical key 构造 symmetry orbit。
5. 所有输出关系执行精确 GiNaC 恒等式验证。

模块独立放在 `symmetry/`，不修改 `expand` 的约化和级数代码。

## 2. 文件组织

```text
symmetry/
├── docs/
│   ├── problem_solution.md
│   └── code_structure.md
├── include/
│   ├── symmetry_types.hpp
│   ├── branch_transform.hpp
│   ├── sector_canonicalizer.hpp
│   ├── symmetry_finder.hpp
│   └── symmetry_formatter.hpp
├── src/
│   ├── branch_transform.tpp
│   ├── sector_canonicalizer.tpp
│   ├── symmetry_finder.tpp
│   └── symmetry_formatter.tpp
├── tools/
│   └── fbi_sector_symmetry/
│       ├── Makefile
│       └── main.cpp
└── tests/
    └── scarecrow/
        ├── Makefile
        └── test_scarecrow.cpp
```

## 3. 数学对象到代码实体

| Phase A 对象 | 代码实体 |
|---|---|
| sector $s$ | `SectorId` |
| branch 置换 $\tau$ | `Permutation branchPermutation` |
| 传播子映射 $\sigma$ | `std::vector<int> sourceToTarget` |
| $g_\tau(X,Y)$ | `BranchCoordinateTransform` |
| sector canonical key | `std::string canonicalKey` |
| $s\to\operatorname{Can}(s)$ | `CanonicalWitness` |
| symmetry orbit | `SymmetryOrbit` |
| $s\to t$ 精确关系 | `SectorMapping` |

## 4. 核心数据类型

### 4.1 `SectorId`

**文件**：`include/symmetry_types.hpp`

```cpp
struct SectorId {
    std::vector<int> bits;
};
```

公共辅助函数：

```cpp
int sectorIndex(const SectorId& sector);
std::string sectorToString(const SectorId& sector);
std::vector<int> activePropagators(const SectorId& sector);
```

### 4.2 `CanonicalWitness`

对应 `problem_solution.md` 第 5 节。

```cpp
struct CanonicalWitness {
    std::vector<int> branchPermutation; // new branch b <- old tau[b]
    std::vector<int> orderedProps;       // canonical positions -> original prop
    std::string key;
};
```

`orderedProps` 明确保存规范矩阵每个位置来自哪个原传播子，避免在恢复 sector 间映射时混淆置换方向。

### 4.3 `SectorCanonicalForm`

```cpp
struct SectorCanonicalForm {
    SectorId sector;
    std::string key;
    std::vector<CanonicalWitness> witnesses;
};
```

所有达到同一最小 key 的 witness 都保留，用于恢复完整 symmetry 和 sector automorphism。

### 4.4 `SectorMapping`

对应 `problem_solution.md` 第 4、7 节。

```cpp
struct SectorMapping {
    SectorId source;
    SectorId target;
    std::vector<int> sourceToTarget;     // full N-sized map, -1 for inactive
    std::vector<int> branchPermutation;
    GiNaC::ex transformedX;
    GiNaC::ex transformedY;
    bool verified = false;
};
```

### 4.5 `SymmetryOrbit`

```cpp
struct SymmetryOrbit {
    SectorId representative;
    std::vector<SectorId> members;
    std::vector<SectorMapping> mappingsToRepresentative;
};
```

## 5. Branch 坐标变换

### 5.1 `BranchTransformBuilder`

**文件**：`branch_transform.hpp`

实现 `problem_solution.md` 第 2.2–2.3 节。

```cpp
class BranchTransformBuilder {
public:
    BranchTransformBuilder(
        const GiNaC::symbol& X,
        const GiNaC::symbol& Y,
        const GiNaC::ex& shiftA,
        const GiNaC::ex& shiftB);

    std::pair<GiNaC::ex, GiNaC::ex> inducedTransform(
        const std::vector<int>& tau) const;

    std::vector<std::vector<int>> allBranchPermutations() const;
};
```

前置条件：

- 第一版要求 `tau.size() == 3`。

后置条件：

- 返回精确的 $g_\tau(X,Y)$；
- 表达式经过 `normal`。

## 6. Sector Canonicalizer

### 6.1 `SectorCanonicalizer`

**文件**：`sector_canonicalizer.hpp`

实现 `problem_solution.md` 第 5 节。

```cpp
class SectorCanonicalizer {
public:
    SectorCanonicalizer(
        const std::vector<std::vector<GiNaC::ex>>& R,
        const std::vector<int>& branchOfProp,
        const GiNaC::symbol& X,
        const GiNaC::symbol& Y,
        const GiNaC::ex& shiftA,
        const GiNaC::ex& shiftB);

    SectorCanonicalForm canonicalize(const SectorId& sector) const;
};
```

### 6.2 固定 $\tau$ 的 Pak 搜索

关键私有状态：

```cpp
struct PartialCandidate {
    std::vector<int> orderedProps;
    std::vector<bool> used;
    std::vector<std::string> prefixTokens;
};
```

关键伪代码：

```text
candidates = {empty candidate}
for canonical position pos:
    destinationBranch = canonicalBranchSequence[pos]
    sourceBranch = tau[destinationBranch]

    expanded = {}
    for candidate in candidates:
        for unused prop in sourceBranch:
            next = candidate + prop
            append R[prop, previousProps] and R[prop, prop] tokens
            expanded += next

    bestPrefix = minimum full prefix in expanded
    candidates = all entries whose prefix == bestPrefix
```

`canonicalBranchSequence` 按 branch 0、1、2 排列，每个 branch 内放置该 sector 在对应源 branch 中的活跃传播子。

该过程依据 `problem_solution.md` 第 5.2 节，与完整枚举相容 $\sigma$ 后取最小序列等价。

### 6.3 表达式序列化

```cpp
std::string expressionKey(const GiNaC::ex& value);
```

处理流程：

1. `normal(value)`；
2. 对分子分母分别 `expand`；
3. 使用 GiNaC 确定输出生成字符串。

canonical key 由 token 长度前缀编码后连接，避免简单分隔符冲突。

## 7. Symmetry Finder

### 7.1 `SymmetryFinder`

**文件**：`symmetry_finder.hpp`

实现 `problem_solution.md` 第 6–8 节。

```cpp
class SymmetryFinder {
public:
    SymmetryFinder(
        const std::vector<std::vector<GiNaC::ex>>& topS,
        int numProps,
        int numBranches,
        const GiNaC::symbol& X,
        const GiNaC::symbol& Y,
        const GiNaC::ex& shiftA,
        const GiNaC::ex& shiftB);

    std::vector<SectorId> enumerateValidSectors() const;
    std::vector<SymmetryOrbit> findOrbits() const;
    bool verify(const SectorMapping& mapping) const;
};
```

构造函数：

1. 验证 `topS` 维度为 `B+N`；
2. 从 incidence block 提取 `branchOfProp`；
3. 提取右下角 $R$；
4. 创建 `SectorCanonicalizer`。

### 7.2 恢复 sector 间映射

若 source witness 的 `orderedProps[k]` 与 target witness 的 `orderedProps[k]` 占据相同 canonical position，则定义

```text
sourceToTarget[sourceProp] = targetProp
```

两个 witness 中的 branch 重标记共同确定 source 到 target 的 branch 置换。实现通过逐 branch 检查传播子映射的归属直接恢复该置换，而不依赖容易出错的置换复合公式。

随后调用 `verify()`。

### 7.3 代表选择

同一 key 下按 `sectorIndex` 降序排列，第一个 sector 为代表。该规则对应 `problem_solution.md` 第 6 节。

## 8. IO 与工具

### 8.1 `fbi_sector_symmetry`

```bash
fbi_sector_symmetry <S_path> <config_path> <output_path>
```

工具复用 `expand/include/io.hpp`：

- `parseMatrixFile` 读取 $S(X,Y)$；
- `parseConfigFile` 读取 `B,N,p,a,b`。

`a,b` 在 config 中是有限域元素。工具使用有理重构恢复小有理数，并检查重构值模 $p$ 后等于原值。若重构失败则报错，不能把有限域代表元直接当作有理平移量。

输出包含：

```text
[orbits]
[sector_mappings]
[automorphisms]
```

每条 mapping 输出 source、target、`sigma`、`tau`、`Xmap` 和 `Ymap`。

## 9. 验证策略

### 9.1 Branch transform

检查六个 $S_3$ 置换：

- 恒等置换给出 $(X,Y)$；
- 交换 branch 2、3 给出 $(X,1-Y-2b)$；
- 每个 $g_\tau$ 映回 branch 参数后与直接置换结果一致。

### 9.2 Pak canonicalization

构造小型带 branch 的对称矩阵：

- 与完整枚举全部相容 $\sigma$ 的最小 key 比较；
- 检查并列 witness 可恢复 automorphism。

### 9.3 精确 mapping

对每条输出关系重新检查完整矩阵恒等式，不依赖 canonical key。

### 9.4 Scarecrow 回归

输入：

```text
expand/runs/scarecrow/S
expand/runs/scarecrow/config
```

检查：

1. 枚举 sector 数与 `B=3,N=4` 的 branch 覆盖条件一致；
2. 所有输出 mapping 的 `verified=true`；
3. orbit 成员无重复且覆盖全部有效 sector；
4. 打印实际找到的不同 sector 关系供领域审查。

## 10. 实现顺序

1. `symmetry_types` 和 branch transform。
2. expression key 和固定 $\tau$ 的 Pak canonicalizer。
3. sector 枚举、orbit 和 mapping 恢复。
4. 精确 verifier。
5. CLI formatter。
6. 小型单元测试。
7. scarecrow 回归。

## 11. 参考

1. [`problem_solution.md`](./problem_solution.md)
2. `expand/include/io.hpp`
3. `expand/include/family.hpp`
4. A. Pak, arXiv:1111.0868.
