# 旧工作归档

2026-11-07 将原根目录的六个 J 相关目录整体移动到这里，保留原名、内部相对结构、未提交文件、数值结果、日志、编译产物以及外部程序副本：

| 原目录 | 归档位置 | 内容 |
| --- | --- | --- |
| `j_integral/` | `j_integral/` | 初期 J 展开、FI 层验证、真空与 DB 测试 |
| `j_integral_bootstrap/` | `j_integral_bootstrap/` | AB 原包及接入试验、DP 全局结果、审计与天河记录 |
| `j_integral_global/` | `j_integral_global/` | 全局数值 J-IBP、物理端点与完整 DE |
| `j_integral_methods/` | `j_integral_methods/` | 固定点 z-IBP、Hermite 等方法对比 |
| `j_integral_sixpoint/` | `j_integral_sixpoint/` | 两圈六点例子及完整 J DE |
| `j_integral_x1/` | `j_integral_x1/` | 另一种分支参数组织方式的独立试验 |

这些目录是历史记录，不能视为新的生产入口。历史日志、JSON 中的绝对路径、提交时哈希和天河文件路径按原样保存，不伪造为新位置。旧脚本中若写死原根目录路径，重放时需显式调整；归档不承诺所有历史命令直接可执行。

移动仅发生于本机。此前天河作业使用其独立部署目录，不读取此处本机源码；本次不移动或改写那些部署与作业。

已验证的参考数据在 `j_integral_bootstrap/reference/` 和 `j_integral_sixpoint/reference/`。DP 的旧 AB 运行含外层逐块调用产生的额外检查，不能作为正确入口的性能结果；修正与数值点审计也保存在归档内。
