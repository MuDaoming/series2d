# DP 全 sector dot<=2 复核

从公共 DP 物理家族输入重新生成全部保留三个 branch、总 dot<=2 的 2952 个目标，不提供旧固定基底，由程序重新选取。

结果仍为 496 维：481 个普通 J、15 个边界函数，所有 DE 闭合。原先记录的 62 个 J>FI 原始 sector 中，每个 sector 的普通 J 数都与这次复核相同。

这是在固定运动学、$d=26/7,p=1000003,t=17$ 的有限域数值 IBP 检查，不是任意高 dot 的完备性证明。没有为机制研究重新做 DP 全局 Blade 约化；FI 计数仍使用已有 maximal-cut 清单，并在 register 中明确标记跨 sector 代表问题。

`request.json` 是完整输入，`j_result.json.gz` 是公共接口输出的无损压缩存档（本地运行仍输出 j_result.json）。该全局检查与轻量第二点验证分开，避免运行 `verify_points.py` 时无意重复较大的计算。
