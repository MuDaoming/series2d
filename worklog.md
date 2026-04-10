26-04-08

- 1. 对于项目中的db例子(实际只有六个传播子且topsector没有主积分)测试了maximal-cut, 用一个sector的series确实能找到此sector的关系. <span style="color:green;">这告诉我们: maximal-cut应该是能行的</span>
- 2. 对于dp测试maximal-cut(参数: deg=200, dot2且3branch积分), 发现一些sector找不到关系, 如topsector, 连0阶关系也即被积函数关系都没有找到. 后面找到原因了, 是因为G太小了, 只有topsector的dot2的45个积分. 我用所有topsector的dot3的165积分就能够搜索到0阶关系. 尽管用这些关系最后解出来FI的在此sector的主积分是22个dot2的积分, 但是需要注意, 用于求解FI的FI relations都是dot3的relation. <span style="color:green;">这告诉我们: 用来搜索主积分的集合应该大于主积分, 譬如主积分是dot2积分的真子集, 可能dot2积分的集合也搜不出来, 要用dot3积分的集合搜索才行.</span>
- 3. <span style="color:blue;">原则上可以在expand的LRR层面就实现maximal-cut, 应该会快不少(吗?), 只在当前sector约化和求解方程</span>
- 4. 目前有一个地方做复杂了, 建立某个积分的微分方程时, R矩阵应该只取当前sector的R矩阵, 现在取了整个R, 会产生更高sector的积分, 这是不必要的. <span style="color:blue;">未来应该修改代码,使得建立方程和约化一样只在当前sector做, 尤其是如果要做maximal-cut</span>
- 5. dp(deg=20,dot3点165个topsector积分)的例子看出来. 只有20个方程, m=0时有165个变量, 能搜出多项式关系(没有检验正确性, 但是能求解出来FI), m=1时有165*2个变量, 也能搜出多项式关系, 但是多项式关系$\delta$取1得到的FI关系没有能够求解FI. 需要想一下为什么.<span style="color:blue;">方程的个数$deg*N_{bc}$与变量的个数$(m+1)*|G|$满足什么条件时能保证求解出来的多项式关系一定没有错?</span>
  
26-04-08

- 1. 给`search_poly_relation`添加了检验机制, 最后一阶会用来检验求解出的关系是否正确. 用新的正确的poly_relation_searcher发现dot3其实也搜不出来.

- 2. dp的bc1求解到了dot3(165个)1000阶, 也没有搜到任何关系.