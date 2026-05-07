SetDirectory[DirectoryName[$InputFileName]];

Get["/home/mudaoming/Documents/Packages/calcloop-main/FBI/FBIReduce.wl"];

replace = {Global`FBIReduce`$D -> d, X[1] -> X, X[2] -> Y, X[3] -> Z};
var = {X, Y, Z};
varrep = {X -> X, Y -> Y (1 - X), Z -> 1 - X - Y (1 - X)};
a = 12/37; b = 3/29;
transrep = {X -> X + a, Y -> Y + b};
p = 2305843009213693951;

(* 1) define family *)
fam = "vac2loop";
loopmom = {l1, l2};
extmom = {};
numrep = {msq -> 1, d -> 7/13};
sprep = {p1^2 -> 0, p2^2 -> msq, p1 p2 -> s};
pdlist = {l1^2 - msq, l2^2 - msq, (l1 + l2)^2 - msq};
Global`FBIReduce`FBIFamilyDefine[fam, pdlist, loopmom, extmom, sprep];

myMod[x_Integer, pp_] := Mod[x, pp];
myMod[x_Rational, pp_] := Mod[Mod[Numerator[x], pp] * PowerMod[Mod[Denominator[x], pp], -1, pp], pp];
myMod[x_?NumericQ, pp_] := Mod[x, pp];
myMod[x_Symbol, pp_] := x;
myMod[x_List, pp_] := (myMod[#, pp] & /@ x);
myMod[x_Power, pp_] := myMod[x[[1]], pp]^x[[2]];
myMod[x_Times, pp_] := Times @@ (myMod[#, pp] & /@ List @@ x);
myMod[x_Plus, pp_] := Plus @@ (myMod[#, pp] & /@ List @@ x);

(* 2) output p and Z_p a,b,S *)
S = (((Global`FBIReduce`FamilyInfo[fam])["SR"][[1]] /. replace) /. numrep);
Sshifted = Expand[(S /. varrep) /. transrep];
Smod = myMod[Sshifted, p] // Expand;
Put[Smod, "S"];

cfg = StringRiffle[{"B = 3", "N = 3", "p = " <> ToString[p], "d = " <> ToString[myMod[(d /. numrep), p]], "a = " <> ToString[myMod[a, p]], "b = " <> ToString[myMod[b, p]]}, "\n"];
Export["config", cfg <> "\n", "Text"];

(* 3) output Z_p MMA DE matrices *)
de = (Global`FBIReduce`genBranchDE[fam, {X[1], X[2], X[3]}] /. replace) // Factor;
master = (Global`FBIReduce`KMasterInfo[fam])[[All, 3]];
U = X Y + Y Z + Z X;
fac = Table[U^(Total[master[[i]]] - ((2 + 1)/2) d), {i, 1, Length@master}];
Table[
  de[[i]] = DiagonalMatrix@Table[D[fac[[j]], var[[i]]] / fac[[j]], {j, 1, Length@master}] +
            Table[de[[i, j, k]] fac[[j]]/fac[[k]], {j, 1, Length@master}, {k, 1, Length@master}],
  {i, 1, Length@var}
];
bDE = <|Thread[var -> de]|>;
DE = <||>;
DE[X] = (Sum[(bDE[v] /. varrep) D[(v /. varrep), X], {v, var}] /. numrep /. transrep) // Factor;
DE[Y] = (Sum[(bDE[v] /. varrep) D[(v /. varrep), Y], {v, var}] /. numrep /. transrep) // Factor;

modRat[expr_] := Module[{t, num, den, dn},
  t = Together[expr];
  num = PolynomialMod[Expand[Numerator[t]], p];
  den = PolynomialMod[Expand[Denominator[t]], p];
  dn = Coefficient[den, X, 0] // Coefficient[#, Y, 0] &;
  If[dn === 0, Return[Indeterminate]];
  num = PolynomialMod[num * PowerMod[dn, -1, p], p];
  den = PolynomialMod[den * PowerMod[dn, -1, p], p];
  If[den === 1, num, num/den]
];
Put[Map[modRat, DE[X], {2}], "AX_mma_modp"];
Put[Map[modRat, DE[Y], {2}], "AY_mma_modp"];
Print["done vac"];
