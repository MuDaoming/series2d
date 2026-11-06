(* This entry performs maximal-cut master discovery only. No BLReduce or DE. *)
$HistoryLength=0;
out=DirectoryName[$InputFileName];SetDirectory[out];
fd=Import[FileNameJoin[{out,"input.json"}],"RawJSON"];
Get[FileNameJoin[{fd["blade_root"],"Blade.wl"}]];
BLNthreads=1;BLSSolveNthreads=1;
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
With[{dir=FileNameJoin[{out,"cache"}]},Blade`Private`CacheDirectory[key_]:=FileNameJoin[{dir,ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}]];
BLSetReducerOptions["IntegralOrdering"->2];
fam=ToExpression["miCount"<>fd["name"]];
prop=ToExpression/@fd["propagators"];loops=ToExpression/@fd["loops"];legs=ToExpression/@fd["external_momenta"];
sp=(ToExpression[#[[1]]]->ToExpression[#[[2]]])&/@fd["scalar_products"];
top=fd["top"];n=fd["physical_count"];
(* Omitting cut means every sector is considered, each on its own maximal cut. *)
BLFamilyDefine[fam,ToExpression[fd["dimension"]],prop,loops,legs,{},sp,top,{}];
actual=Blade`Private`CacheDirectory["sectormappings"];
If[!StringStartsQ[actual,FileNameJoin[{out,"cache"}]],Print["CACHE_ISOLATION_FAILED"];Exit[2]];
Print["MI_ONLY_ENTRY|",fd["name"],"|",actual];
start=AbsoluteTime[];raw=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];
uniquePath=FileNameJoin[{actual,"results","rawmasters"}];
If[!FileExistsQ[uniquePath]||!ListQ[raw],Print["MAXIMALCUT_DISCOVERY_FAILED"];Exit[3]];
unique=Get[uniquePath];
If[!ListQ[unique]||!And@@(MatchQ[#,BL[_,{_Integer..}]]&/@Join[raw,unique]),Print["INVALID_MASTER_RECORD"];Exit[4]];
indices=#[[2]]&/@unique;rawindices=#[[2]]&/@raw;
topQ[v_]:=And@@Thread[Take[v,n]>0];
sectors=Counts[(StringJoin[ToString/@(Boole[#>0]&/@Take[#,n])])&/@indices];
result=<|"method"->"Blade BLMaximalCutMasters only; all sectors, symmetry-mapped basis", "family"->fd["name"],"dimension"->fd["dimension"],"physical_propagators"->n,"auxiliary_propagators"->(Length[top]-n),"master_count"->Length[unique],"raw_nomap_count"->Length[raw],"topsector_master_count"->Count[topQ/@indices,True],"topsector_raw_nomap_count"->Count[topQ/@rawindices,True],"master_indices"->indices,"topsector_master_indices"->Select[indices,topQ],"sector_counts"->sectors,"rank_free"->And@@(And@@Thread[#>=0]&/@indices),"seconds"->(AbsoluteTime[]-start),"actual_cache"->actual|>;
Export[FileNameJoin[{out,"masters.json"}],result,"RawJSON"];Put[unique,FileNameJoin[{out,"masters.wl"}]];
Print["MI_ONLY_DONE|",fd["name"],"|total=",Length[unique],"|top=",result["topsector_master_count"]];Exit[0];
