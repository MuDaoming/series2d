26-04-08

- 1. 对于项目中的db例子(实际只有六个传播子且topsector没有主积分)测试了maximal-cut, 用一个sector的series确实能找到此sector的关系
- 2. 对于dp测试maximal-cut(参数: deg=200, dot2且3branch积分), 发现一些sector找不到关系, 如topsector, 连0阶关系也即被积函数关系都没有找到
- 3. 原则上可以在expand的LRR层面就实现maximal-cut
- 4. 目前有一个地方做复杂了, 建立某个积分的微分方程时, R矩阵应该只取当前sector的R矩阵, 现在取了整个R, 会产生更高sector的积分, 这是不必要的