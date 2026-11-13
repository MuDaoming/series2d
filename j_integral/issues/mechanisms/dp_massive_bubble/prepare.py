import json,itertools
from pathlib import Path
out=Path(__file__).resolve().parent
for label,mass in [('equal','1'),('split','5/4')]:
 d=out/label;d.mkdir(exist_ok=True);targets=[]
 for first in [[1,0],[0,1],[1,1]]:
  nu=first+[1,1];active=[i for i,v in enumerate(nu) if v]
  for degree in range(3):
   for dots in itertools.combinations_with_replacement(active,degree):
    q=nu.copy()
    for i in dots:q[i]+=1
    targets.append(q)
 f={'schema_version':1,'family':{'name':'dp_massive_bubble_'+label,'loops':['l1','l2'],'external_momenta':['P','Q'],'gram':[['1','5/24'],['5/24','1']],'propagators':[{'loop':[1,0],'shift':[1,0],'mass2':1},{'loop':[1,0],'shift':[0,0],'mass2':1},{'loop':[0,1],'shift':[0,0],'mass2':mass},{'loop':[1,1],'shift':[0,1],'mass2':1}]},'j':{'coordinates':'tz','branches':[0,0,1,2],'boundary':{'z0':'fbi_analytic','z1':'regulated_zero'}},'point':{'prime':1000003,'dimension':'26/7','t':'17','parameters':{}},'request':{'targets':targets,'derivatives':'basis'},'options':{'completion':16,'germ_order':16}}
 (d/'request.json').write_text(json.dumps(f,indent=2)+'\n')
 text='''$HistoryLength=0;out=DirectoryName[$InputFileName];SetDirectory[out];
Get["/home/mudaoming/Documents/Packages/blade/blade-main/Blade.wl"];
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
Blade`Private`CacheDirectory[key_]:=FileNameJoin[{out,"cache",ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}];
BLNthreads=1;BLSSolveNthreads=1;BLSetReducerOptions["IntegralOrdering"->2];
'''
 text+='prop={(l1+P)^2-1,l1^2-1,l2^2-('+mass+'),(l1+l2+Q)^2-1,(l1+Q)^2,(l2+P)^2,(l2+Q)^2};sp={P^2->1,Q^2->1,P Q->5/24};top={1,1,1,1,0,0,0};\n'
 text+='BLFamilyDefine[dpbubble'+label+',26/7,prop,{l1,l2},{P,Q},{},sp,top,{},top];\n'
 text+='start=AbsoluteTime[];mi=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];Put[mi,"masters.wl"];Export["result.json",<|"masters"->ToString[mi,InputForm],"count"->Length[mi],"seconds"->(AbsoluteTime[]-start)|>,"RawJSON"];Print["RESULT=",mi];Exit[];\n'
 (d/'count.wl').write_text(text)
