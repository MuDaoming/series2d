import json,itertools
from pathlib import Path
import sympy as s
root=Path(__file__).resolve().parents[4];out=Path(__file__).resolve().parent
f=json.load(open(root/'j_integral/tests/e6/input/family_request.json'));G=s.Matrix([[s.Rational(x) for x in row] for row in f['family']['gram']]);V=s.Matrix([[0,0,1,0,0],[0,0,0,0,1],[1,1,-1,-1,0],[1,0,0,-1,0]]);gram=V*G*V.T
for label,p2,k2 in [('both_null',0,0),('p_off',s.Rational(2,7),0),('both_off',s.Rational(2,7),s.Rational(3,11))]:
 d=out/label;d.mkdir(exist_ok=True);g=gram.copy();g[0,0]=p2;g[1,1]=k2
 props=[{'loop':l,'shift':q,'mass2':0} for l,q in [([1,0],[0,0,0,0]),([1,0],[1,0,0,0]),([0,1],[0,0,0,0]),([0,1],[0,1,0,0]),([0,1],[0,0,1,0]),([1,1],[0,0,0,1])]]
 targets=[]
 for bits in itertools.product([0,1],repeat=5):
  if not any(bits[:2]) or not any(bits[2:]):continue
  nu=list(bits)+[1];targets.append(nu)
  active=[i for i,v in enumerate(nu) if v]
  for order in range(1,2 if label=='both_null' else 3):
   for dots in itertools.combinations_with_replacement(active,order):
    other=nu.copy()
    for i in dots:other[i]+=1
    targets.append(other)
 request={'schema_version':1,'family':{'name':'six_line_'+label,'loops':['l1','l2'],'external_momenta':['p','k','q','r'],'gram':[[str(g[i,j]) for j in range(4)] for i in range(4)],'propagators':props},'j':{'coordinates':'tz','branches':[0,0,1,1,1,2],'boundary':{'z0':'fbi_analytic','z1':'regulated_zero'}},'point':{'prime':1000003,'dimension':'26/7','t':'17','parameters':{}},'request':{'targets':targets,'derivatives':'basis'},'options':{'completion':16,'germ_order':16}}
 (d/'request.json').write_text(json.dumps(request,indent=2)+'\n')
 names=['p','k','q','r'];sp=[f'{names[i]} {names[j]}->{g[i,j]}' for i in range(4) for j in range(i,4)]
 text='''$HistoryLength=0;out=DirectoryName[$InputFileName];SetDirectory[out];
Get["/home/mudaoming/Documents/Packages/blade/blade-main/Blade.wl"];
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
Blade`Private`CacheDirectory[key_]:=FileNameJoin[{out,"cache",ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}];
BLNthreads=1;BLSSolveNthreads=1;BLSetReducerOptions["IntegralOrdering"->2];
prop={l1^2,(l1+p)^2,l2^2,(l2+k)^2,(l2+q)^2,(l1+l2+r)^2,(l1+k)^2,(l1+q)^2,(l1+r)^2,(l2+p)^2,(l2+r)^2};
top={1,1,1,1,1,1,0,0,0,0,0};
'''
 text+='sp={'+','.join(sp)+'};\n'
 text+='BLFamilyDefine[six'+label.replace('_','')+',26/7,prop,{l1,l2},{p,k,q,r},{},sp,top,{},top];\n'
 text+='start=AbsoluteTime[];mi=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];Put[mi,"masters.wl"];Export["result.json",<|"masters"->ToString[mi,InputForm],"count"->Length[mi],"seconds"->(AbsoluteTime[]-start)|>,"RawJSON"];Print["RESULT=",mi];Exit[];\n'
 (d/'count.wl').write_text(text)
