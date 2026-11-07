# DP 原包阅读与可复现诊断

数学对象、上游代码流程及实测解释见 [package_review.md](../docs/package_review.md)。所有脚本只写 `results/dp_diagnosis/`，不修改 canonical 结果、主脚本或上游源码。最终盲点 `dp_verification/` 不在这些脚本的数据源中。

从项目根目录运行：

```bash
python j_integral_bootstrap/diagnosis/validate_new_points.py
```

这会验证运行时冻结的 canonical 已知 DE/R 与当前 `dp_more_full/points/`，包含隐含零及基底恒等约化。检查点快照、SHA256、实际点列表和不一致明细写入 `new_points_validation.json`。每次重跑覆盖本诊断目录中的上次快照与报告；如需保留历史请先复制该目录。

## 脚本用途

| 脚本 | 输入与用途 | 主要输出 |
| --- | --- | --- |
| `check_inputs.py` | 冻结 snapshot 的已知非零 DE、cut 对角和 exceptional 对角，对已有学习点逐项检查 | `input_validation.json` |
| `check_known_zeros.py` | 同一 snapshot 已知 DE 的隐含零检查 | `known_zeros.json` |
| `validate_new_points.py` | 当前 canonical DE/R 对新增完整点逐项验证，包括零；不修改原结果 | `new_points_validation.json`、两个冻结快照 |
| `profile_block.py` | 复制当前主检查点，在首个未解块隔离默认 R 和额外 R=4、12；成本包含主驱动输入与极点学习 | `profile_block.json`、`snapshot/` |
| `capture_block.py` | 同上，保存首个实际 AB 输入 A/B/包络/fixed/训练点作为胶囊，并单独检验 fixed 恒等式相容性 | `capsule.json`、`capture_profile.json` |
| `profile_capsule.py` | 直接重放已冻结胶囊，无需读取后续变化的主检查点；R 对照的推荐重放入口 | `profile_capsule.json` |
| `normalize_capsule.py` | 历史隔离：固定包络只改变共同标量规范；包含一次很宽松非端点试验，耗时较大 | `normalization.json` |
| `check_new_coefficients.py` | pole_route 已重建真实 C 对胶囊包络逐条检验分母整除和 cap | `new_coefficients_envelope.json` |
| `new_fixed.py` | 合并新增真实 C，以精确分母/cap 固定，检验是否与原包解析 ansatz 相容 | `new_fixed.json` |
| `residues.py` | A22 非端点线性因子，一阶极点残差的低负整数本征值 | `negative_residue_eigenvalues.json` |
| `integer_residues.py` | 同上，整数本征值范围扩大至 −15…60 | `integer_residue_eigenvalues.json` |
| `top_integer_residues.py` | A11 非端点线性因子的对应测试 | `top_integer_residue_eigenvalues.json` |
| `nonlinear_residues.py` | A22 非线性因子的残差矩阵按剩余域展开后检查整数本征值 | `nonlinear_residues.json` |

`profile_block.py` 和 `capture_block.py` 会跟随当前主驱动、当前检查点和当前 pole-learning 数据；当主检查点继续推进时，它们不一定再捕获旧的六维块。旧实验应使用 `capsule.json` 加 `profile_capsule.py` 重放。所有计时报告保留原包每个 FactorBlock 的 metadata，分别对应 $N_{A_{22},Q}$、$V_{A_{22}}$、$P_{A_{11}}$、Hom、支撑、坐标和 global cap。

原包回归测试未改动：`upstream_regression.log` 是四项 public PL 测试，`extension_regression.log` 是两项 structural extension 测试，全部通过。前者包括原包自带的迟发共振反例。

## 当前证据与未证明事项

已知 A/B 与原学习点一致；隐含零也一致。冻结六维块的 R、R+4、R+12 产生相同 7175 维 ansatz 且实际点值不相容。增加端点规范幂未解决；对所有其余候选因子统一放松后，只得到严重欠定的相容问题，不能当作根因证明。

已有新增真 C 子集满足旧包络，并与解析条件相容；尚未补齐的条目仍可能存在包络问题。整数残差检查覆盖指定范围，不能证明无限导数完备性，不能替代 J 的端点区域分析。未得到足以唯一归因“缺极点”“有限 R”“解析规范”的直接真 C 反例。

对新增完整点的当前冻结验证为：339 条已知 DE、2404 个已知 R 对象，在 t=58…62 全部一致；DE 840720 分量，R 5961920 分量，包括隐含零。冻结 SHA256 见 JSON 和 `docs/package_review.md`。这些新点可以参与后续学习，故它们不是最终隔离盲点；完整 DP 是否完成仍由 canonical 状态及最终盲验证决定。

## 已解决的历史六维 × 47 维阻塞

新增两个纯代数诊断入口：`restricted_envelope.py` 重建历史候选并检验真分母整除；`restricted_repair.py` 对真 C 逐因子做精确主部成员检查，并用原 65 点验证单因子升幂的修复。输出 `results/dp_diagnosis/restricted/`。结论和复现细节见 [RESTRICTED_DP_RESOLUTION.md](RESTRICTED_DP_RESOLUTION.md)：$t=1/2$ 极点幂少一阶，原 65 点实际上足够；没有读取 fresh 验证点。
