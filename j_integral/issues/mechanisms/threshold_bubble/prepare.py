import json,itertools
from pathlib import Path
out=Path(__file__).resolve().parent;targets=[]
for tail in [[1,0],[0,1],[1,1]]:
 nu=[1,1]+tail;active=[i for i,v in enumerate(nu) if v]
 for degree in range(3):
  for dots in itertools.combinations_with_replacement(active,degree):
   q=nu.copy()
   for i in dots:q[i]+=1
   targets.append(q)
f={'schema_version':1,'family':{'name':'massive_massless_threshold_bubble','loops':['l1','l2'],'external_momenta':['Q'],'gram':[['2']],'propagators':[{'loop':[1,0],'shift':[0],'mass2':1},{'loop':[0,1],'shift':[0],'mass2':0},{'loop':[1,1],'shift':[0],'mass2':1},{'loop':[1,1],'shift':[1],'mass2':1}]},'j':{'coordinates':'tz','branches':[0,1,2,2],'boundary':{'z0':'fbi_analytic','z1':'regulated_zero'}},'point':{'prime':1000003,'dimension':'26/7','t':'17','parameters':{}},'request':{'targets':targets,'derivatives':'basis'},'options':{'completion':16,'germ_order':16}}
(out/'request.json').write_text(json.dumps(f,indent=2)+'\n')
text='''$HistoryLength=0;out=DirectoryName[$InputFileName];SetDirectory[out];
Get["/home/mudaoming/Documents/Packages/blade/blade-main/Blade.wl"];
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
Blade`Private`CacheDirectory[key_]:=FileNameJoin[{out,"cache",ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}];
BLNthreads=1;BLSSolveNthreads=1;BLSetReducerOptions["IntegralOrdering"->2];
prop={l1^2-1,l2^2,(l1+l2)^2-1,(l1+l2+Q)^2-1,(l1+Q)^2};top={1,1,1,1,0};
BLFamilyDefine[thresholdbubble,26/7,prop,{l1,l2},{Q},{},{Q^2->2},top,{}];
masters=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];Put[masters,"masters.wl"];
target={BL[thresholdbubble,{1,1,1,1,0}],BL[thresholdbubble,{1,1,0,2,0}],BL[thresholdbubble,{1,0,1,1,0}]};
res=BLReduce[Join[target,masters],masters,"threshold_verify","BladeMode"->None];
If[!ListQ[res]||Length[res]!=Length[target]+Length[masters],Print["REDUCE_FAILED"];Exit[2]];
res=res/.Thread[Take[res,-Length[masters]]->masters];
d=26/7;identity=Together[res[[1]]+res[[2]]/(d-3)-(d-2)/(2(d-3))*res[[3]]];
If[identity=!=0,Print["IDENTITY_FAILED=",identity];Exit[3]];
Export["blade_checks.json",<|"dimension"->"26/7","Q2"->2,"masters"->ToString[masters,InputForm],"topsector_master_count"->Count[(And@@Thread[Take[#[[2]],4]>0]&/@masters),True],"identity_exactly_zero"->True,"reductions"->ToString[Take[res,3],InputForm]|>,"RawJSON"];Print["THRESHOLD_IDENTITY_PASS"];Exit[];
'''
(out/'verify.wl').write_text(text)
