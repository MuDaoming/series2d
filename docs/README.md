# Docs 目录管理规范

## 1. 概述

本文档定义 docs 目录的组织方式，用于系统性地记录项目中的领域知识，以及问题定义与解决方案。

```
docs/
├── background/           # 领域知识、术语定义
└── problem_solution/     # 问题定义与数学/逻辑层面的解决方案
```

## 2. 文件命名规范

### 2.1 格式

```
<前缀><编号>_<简短描述>.md
```

编号按创建顺序递增，不体现层级关系。

### 2.2 前缀定义

| 前缀 | 含义 | 示例 |
|------|------|------|
| B | Background - 背景知识 | `B001_vector_space.md` |
| PS | Problem & Solution - 问题与解决方案（合并在同一文件） | `PS001_reduce_2d_series.md` |

## 3. 文件状态

### 3.1 状态类型

| 状态 | 含义 |
|------|------|
| `active` | 当前使用中 |
| `deprecated` | 已废弃 |
| `paused` | 暂停 |
| `completed` | 已完成 |

### 3.2 Front Matter 字段

```yaml
---
id: PS001
title: Reduce 2D Series
status: active
created: 2024-01-15
updated: 2024-02-20
depends:
  - B001
  - PS002
parts:
  - PS001.1  # 子问题，见 xxx 文件
  - PS001.2
notes: |
  2024-02-20: 修订了符号约定
---
```

- **depends**：本文档依赖的前置文档
- **parts**：本文档引用的子部分文档（用于拆分大问题/大方案）

## 4. 引用规范

使用相对路径引用其他文档：

```markdown
基于 [PS001: Reduce 2D Series](./problem_solution/PS001.md) 的定义。
```

## 5. 各目录内容定义

### 5.1 background/

存放阅读项目文档前必须了解的领域知识。术语定义、领域基础概念、本项目使用的数学符号约定等。可以包含 Markdown 文件和 PDF 文件。基础知识长期有效，通常不标记状态。

### 5.2 problem_solution/

Problem & Solution。问题定义与数学/逻辑层面的解决方案放在同一文件中。

**文件内容：** 问题的背景与动机、形式化定义、约束条件、解决方案的原理推导、数学证明（如果需要）、正确性论证、适用范围与局限性。

一个问题可以拆分为多个子问题，通过 `parts` 字段引用。多个 PS 可以组合使用，通过 `depends` 字段表达依赖。
