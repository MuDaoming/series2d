(* DP: only maximal cut, dot-only MI ordering, small target batches. *)
base=DirectoryName[DirectoryName[$InputFileName]];out=FileNameJoin[{base,"results","dp_fi_ibp"}];
If[!DirectoryQ[out],CreateDirectory[out,CreateIntermediateDirectories->True]];
Get["/home/mudaoming/Documents/Packages/blade/blade-main/Blade.wl"];
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
Blade`Private`CacheDirectory[key_]:=FileNameJoin[{out,"blade_cache",ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}];
$BLWorkingDirectory=FileNameJoin[{out,"blade_cache"}];BLNthreads=2;BLSSolveNthreads=2;
BLSetReducerOptions["IntegralOrdering"->2];
prop={l1^2-1,(l1+p1)^2-1,(l1-p2)^2-1,l2^2-1,(l2-p3)^2,(l2+p4)^2,
 (l1+l2+p1-p3)^2-1,(l1+l2-p2+p4)^2-1,
 (l1+p3)^2,(l1+p4)^2,(l2+p1)^2};
sp={p1^2->0,p2^2->0,p3^2->1,p4^2->1,p1 p2->1/2,p1 p3->-1/4,p1 p4->-1/6,p2 p3->-1/3,p2 p4->3/8,p3 p4->-7/8};
top=Join[ConstantArray[1,8],{0,0,0}];
BLFamilyDefine[jdpcutcheck,26/7,prop,{l1,l2},{p1,p2,p3,p4},{},sp,top,{},top];
masters=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];
Put[masters,FileNameJoin[{out,"blade_masters.wl"}]];
If[!ListQ[masters]||masters==={},Print["INVALID_MASTERS=",masters];Exit[1]];
rankFree=And@@(MatchQ[#,BL[_,{_Integer..}]]&&And@@Thread[Take[#[[2]],8]>=1]&&Drop[#[[2]],8]=={0,0,0}&/@masters);
If[!rankFree,Print["RANK_MASTER_DETECTED=",masters];Exit[1]];
Export[FileNameJoin[{out,"blade_master_info.json"}],<|"dimension"->"26/7","prime"->1000003,"master_indices"->(Take[#[[2]],8]&/@masters),"master_count"->Length[masters],"max_master_dot"->Max[(Total[Take[#[[2]],8]]-8)&/@masters],"rank_free"->rankFree,"cut"->top|>,"RawJSON"];
Print["DOT_ONLY_MASTER_COUNT=",Length[masters]];
(* First batch is only corner and the eight single dots; never the full family. *)
indices=Join[{ConstantArray[1,8]},Table[ConstantArray[1,8]+UnitVector[8,i],{i,8}]];
target=BL[jdpcutcheck,Join[#,{0,0,0}]]&/@indices;
res=BLReduce[Join[target,masters],masters,"single_dot_cut","BladeMode"->None];
If[!ListQ[res]||Length[res]!=Length[target]+Length[masters],Print["REDUCTION_FAILED"];Exit[1]];
res=res/.Thread[Take[res,-Length[masters]]->masters];
Put[Thread[Join[target,masters]->res],FileNameJoin[{out,"blade_rules_single_dot.wl"}]];
p=1000003;modr[x_]:=Mod[Numerator[Together[x]]PowerMod[Mod[Denominator[Together[x]],p],-1,p],p];
labels=("FI{"<>StringRiffle[ToString/@Take[#[[2]],8],","]<>"}")&/@masters;
records=MapThread[Function[{nu,v},<|"indices"->nu,"coefficients"->AssociationThread[labels,modr/@(Coefficient[Expand[v],#]&/@masters)],"unresolved"->(Expand[v/.Thread[masters->ConstantArray[0,Length[masters]]]]=!=0)|>],{indices,Take[res,Length[indices]]}];
Export[FileNameJoin[{out,"blade_modp_single_dot.json"}],<|"prime"->p,"dimension"->"26/7","master_indices"->(Take[#[[2]],8]&/@masters),"rank_free"->True,"cut"->top,"records"->records|>,"RawJSON"];
Print["BLADE_CUT_BATCH_DONE"];
