# 代码规范

## 命名规范

| 类型 | 规范 | 示例 |
|------|------|------|
| 文件名 | snake_case，头文件 `.hpp`，实现 `.tpp` | `series_solver.hpp` |
| 类名 | PascalCase | `SeriesSolver` |
| 成员变量 | camelCase + 下划线后缀 | `numProps_`, `targetDeg_` |
| 成员函数 | camelCase | `getCase()`, `solveLRRAtDeg()` |
| 模板参数 | 大写缩写 | `RT`, `PT`, `ST` |
| 局部变量 | camelCase | `nuSum`, `maxIndex` |

## 模板参数说明

| 参数 | 符号计算 | 数值计算 | 用途 |
|------|----------|----------|------|
| `RT` | `GiNaC::ex` | `Rational<FlintMod>` | 有理函数 |
| `PT` | `GiNaC::ex` | `Polynomial<FlintMod>` | 多项式 |
| `ST` | `GiNaC::ex` | `FlintMod` | 标量 |

## 新增文件时

1. 遵循 snake_case 命名
2. 头文件放 `include/`，实现放 `src/`
3. 更新 `docs/code_structure.md` 的文件组织部分

## 详细规范

完整的代码规范请参考 [`docs/code_structure.md`](docs/code_structure.md) 第2节。
