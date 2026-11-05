# 固定 X1、先积分 X2 的独立试验

本目录单独测试参数组织方式，不覆盖 `j_integral/results` 的原 t,z 测试。使用同一套通用 C++ FBI 引擎，坐标映射和 Jacobian 从本目录输入读取；编译产物、文档、输入和结果均存放在本目录。

$$
X_1=x,\qquad X_2=(1-x)y,\qquad X_3=(1-x)(1-y),\qquad x,y\in[0,1].
$$

程序内部变量名 `t,z` 分别表示这里的 x,y。定义、边界和比较方法见 [problem_solution](docs/problem_solution.md)，实测结果见 [testing](docs/testing.md)。

```bash
make -C j_integral_x1/native all test
j_integral_x1/build/j_series j_integral_x1/examples/dp_topcut.json \
  j_integral_x1/results/dp_topcut/series.json --max-dot 2 --boundaries --order 200 --verify-all
j_integral_x1/build/j_search j_integral_x1/results/dp_topcut/series.json \
  j_integral_x1/results/dp_topcut/search.json --degree 8 --train 120 --de-degree 60
```
