# 原 t,z 参数化下，直接处理 J 约化与 DE 的方法比较

此目录比较不同算法；继续使用原参数化

$$
X_1=(1-z)t,\qquad X_2=(1-z)(1-t),\qquad X_3=z.
$$

不改动 `j_integral_x1` 的试验，也不覆盖原 `j_integral/results`。核心实现均为通用 C++，传播子、运动学和维度仍从原 JSON 输入读取。

已实现并比较四条路线：原 Taylor 展开搜索、固定 t 点的结构化 z-IBP 后重建、固定 t 点的总导数 ansatz 后重建、保留符号 t 的直接总导数 ansatz。对象定义和推导见 [problem_solution](docs/problem_solution.md)，结果与局限见 [testing](docs/testing.md)，机器可读记录见 [summary](results/summary.json)。

```bash
make -C j_integral_methods/native -j2 all baseline test
j_integral_methods/build/direct_j j_integral/examples/dp_topcut.json \
  j_integral_methods/results/dp/direct.json --method hermite --max-dot 2
j_integral_methods/build/verify_j_points j_integral/examples/dp_topcut.json \
  j_integral_methods/results/dp/direct.json j_integral_methods/results/dp/direct_points.json
```

直接程序同时输出 J 的约化和 DE，**不读取旧约化、旧 DE 或 Taylor 系数**。两种取样路线都可更改 `--max-points`；`--degree-limit` 控制降极点核关系或 ansatz 范围，具体含义见文档。

完整串行比较及第二有限域回归：

```bash
python3 j_integral_methods/tests/run_experiments.py --repeats 3
```

这个 Python 文件只负责启动 C++、检查结果和收集计时，不实现积分约化。
