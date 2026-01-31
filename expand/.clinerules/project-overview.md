# 项目概述

## 项目目标

本项目用于计算二圈Feynman积分（FBI）的二维幂级数展开，在质数域 $\mathbb{Z}_p$ 下使用**逐阶递推**方法。

## 核心文档

开始工作前，请先了解以下文档：

| 文档 | 内容 | 何时参考 |
|------|------|----------|
| [`docs/problem_and_workflow.md`](docs/problem_and_workflow.md) | 问题定义、约化公式、计算流程 | 了解算法逻辑 |
| [`docs/code_structure.md`](docs/code_structure.md) | 代码规范、类结构、核心实现 | 了解代码架构 |
| [`docs/todo.md`](docs/todo.md) | 待办事项和进度 | 查看当前任务 |

## 核心类

- **FlintMod**：有限域元素
- **Polynomial**：二变量多项式
- **Series**：二维幂级数
- **Sector**：处理单个sector的RREF和C/z计算
- **Family**：管理所有sector
- **SeriesSolver**：级数求解器（核心算法）

## 目录结构

```
expand/
├── docs/           # 文档
├── include/        # 头文件（.hpp）
├── src/            # 实现文件（.tpp）
└── test/           # 测试代码
```
