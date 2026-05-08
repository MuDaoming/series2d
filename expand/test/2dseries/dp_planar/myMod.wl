(*myMod能够处理Expand之后的X,Y,s的多项式, 将其有理数系数转化到有限域*)
(*如果要比较两个有理式在有限域下是否相同, 记得先把分母的常数归一化为1*)
p = 2305843009213693951;
Clear[myMod];
Attributes[myMod] = {Listable};
myMod[a_., p_] /; IntegerQ[a] := Mod[a, p];
myMod[Rational[a_., b_.], p_] /; IntegerQ[a] && IntegerQ[b] := 
  Mod[Mod[a, p] PowerMod[b, -1, p], p];
myMod[x_Plus, p_] := Plus @@ (myMod[#, p] & /@ List @@ x);
myMod[Times[x___, X^n_., y___], p_] := myMod[Times[x, y], p] X^n;
myMod[Times[x___, Y^n_., y___], p_] := myMod[Times[x, y], p] Y^n;
myMod[Times[x___, s^n_., y___], p_] := myMod[Times[x, y], p] s^n;
myMod[X^n_., p] := X^n;
myMod[Y^n_., p] := Y^n;
myMod[s^n_., p] := s^n;