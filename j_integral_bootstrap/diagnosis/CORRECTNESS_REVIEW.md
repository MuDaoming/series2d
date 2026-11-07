# 运行边界正确性审查与修复

范围：`bootstrap_global.py`、`bootstrap_reductions.py`、`native_fit.py`、两个 native verifier。审查不读取 fresh R 数据，不修改科学结果或 DP checkpoint。

## 已复现并修复

1. **约化完成状态只比较集合大小。** 原 `len(done)==len(targets)` 会将未列为取样目标的基底主对象与缺失目标相抵。例如基底为 A、B，目标只有 A、T；done={A,B} 时虽然没有 T 的约化，也会标成 complete。现在直接计算 `missing_targets`；有缺失就不完成。反向情形（全部目标已完成但 done 还含额外主对象）也不会再被误判为未完成。
2. **空验证可产生成功证书。** 两个 verifier 原本接受空 points 并输出 all_pass=true、0 检查；缺少 status 也默认 complete。现在 status 必须显式为 complete，points 必须是非空数组；约化 verifier 另外要求 tested_objects 非空。原有逐分量检查仍包括隐含零，没有缩减验证范围。
3. **有理恒等式的节点只按整数去重。** 原 `[0,p]` 被当成两个点，实际同为 $\mathbb F_p$ 的 0。cap=1 的固定零函数可以因此被错误拟合成 x。现在先模 p 规范化，再去重，cap+1 计数对应真正不同的域元素。

修改前的直接反例输出保留于 `results/dp_diagnosis/review/counterexamples.json`。`review_reproductions.py` 是修改前的反例脚本；修复后的可重复回归是下列 unittest。

## 回归

- `test_native_fit.py`：新增同余节点不算两个点，足够真正不同的节点仍正常固定。
- `test_completion_guards.py`：覆盖“集合等大但缺目标”和“目标全部完成但含额外主对象”。
- `test_verify_connection.py`：有效系数、隐含零错误、缺失 status、空/非数组 manifest。
- `test_verify_reductions.py`：既有主对象恒等式、零目标及错误系数回归；新增缺失 status、空 manifest、空目标集合。
- 既有 `test_reuse_unique_fit.py`、`test_fixed_fit.py` 也一起运行，确认原生唯一解复用和固定有理条目行为保持正确。

native verifier 重新编译时使用构建锁与临时可执行文件，编译成功后原子替换二进制，并记录源码、依赖及二进制 SHA256。不修改任何既有科学验证报告。

## 仍建议补充的输入一致性检查

审查时 `bootstrap_global.py` 的 checkpoint 恢复只检查 prime 与 master 名称，不检查 t0；diagonal 输入也只检查 prime。`bootstrap_reductions.py` 的 checkpoint 同样没有核对 t0。因而若有人复用同名工作目录但改变展开中心，旧 numerator_tau/denominator_tau 可能按新中心解释。最小检查是改变 t0 后重用 checkpoint，应在恢复前拒绝；物理输入变化还应有 family/kinematics 指纹。本次授权的三个已实证修复没有顺带修改这些入口；当前 DP 的输入中心固定，不能据此推断现有科学结果有误。

本次回归结果：拟合相关 5 项、完成状态 2 项、native verifier 3 项，共 10 项全部通过；相关 Python 文件 py_compile 通过。两个可执行文件的当前源码/依赖 SHA256、二进制 SHA256 均与 build receipt 一致。日志在 `results/dp_diagnosis/review/{fit,completion,verifier}_regression.log`。
