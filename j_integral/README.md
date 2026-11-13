# J 积分

本项目把两圈 FI 的 FBI 表示先对一个分支参数积分，得到 $J(t)$，再处理最后的 $t$ 积分。当前新的生产代码只做第一重 $J$ 的数值 IBP。

- [src/IBP/README.md](src/IBP/README.md)：当前唯一的数值 IBP 入口、数学定义、输入输出及运行方法。
- [tests/README.md](tests/README.md)：DP、e6使用原包AB命令及新IBP接口的测试。
- [issues/README.md](issues/README.md)：单独定义、研究待解决的问题。
- [archive/README.md](archive/README.md)：所有旧实现、算法试验、计算结果和运行记录。归档代码不参与新程序编译。

目录结构：

```text
j_integral/
  README.md
  src/
    IBP/
      README.md
      main.cpp
      Makefile
      kernel/       # 通用 FBI 与 z-IBP 算法
      examples/     # 不同物理家族的输入；程序中不写拓扑或运动学
      tests/
  issues/
    e6_sector_0110101/
  archive/
    j_integral/
    j_integral_bootstrap/
    j_integral_global/
    j_integral_methods/
    j_integral_sixpoint/
    j_integral_x1/
```
