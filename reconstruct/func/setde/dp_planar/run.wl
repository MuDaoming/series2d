SetDirectory[DirectoryName[$InputFileName]];

Get["/home/mudaoming/Documents/Packages/calcloop-main/FBI/FBIReduce.wl"];

replace = {Global`FBIReduce`$D -> d, X[1] -> X, X[2] -> Y, X[3] -> Z};
varrep = {X -> X, Y -> Y (1 - X), Z -> 1 - X - Y (1 - X)};
a = 12/37; b = 3/29;
transrep = {X -> X + a, Y -> Y + b};
p = 2305843009213693951;

(* 1) define family *)
fam = "dpplanar";
loopmom = {k1, k2};
extmom = {p1, p2, p3, p4};
sprep = {
  p1^2 -> 1,
  p1 p2 -> -5/2,
  p1 p3 -> -2,
  p1 p4 -> 31/2,
  p2^2 -> 2,
  p2 p3 -> -9/2,
  p2 p4 -> -5/2,
  p3^2 -> 1,
  p3 p4 -> -13/2,
  p4^2 -> 2
};
pdlist = ({k1, k1 + p1, k1 + p1 + p2, k1 - k2, k2, k2 + p1 + p2, k2 + p1 + p2 + p3, k2 + p1 + p2 + p3 + p4}^2);
numrep = {d -> 4 - 2*(1/7)};
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

cfg = StringRiffle[{"B = 3", "N = 8", "p = " <> ToString[p], "d = " <> ToString[myMod[(d /. numrep), p]], "a = " <> ToString[myMod[a, p]], "b = " <> ToString[myMod[b, p]]}, "\n"];
Export["config", cfg <> "\n", "Text"];

Print["done dp_planar S/config"];
