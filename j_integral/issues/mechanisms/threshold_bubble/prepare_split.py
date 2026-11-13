import json
from pathlib import Path
base=Path(__file__).resolve().parent;out=base/'split_parent';out.mkdir(exist_ok=True)
f=json.load(open(base/'request.json'));f['family']['propagators'][2]['mass2']='5/4';f['family']['name']='bubble_threshold_mismatch';(out/'request.json').write_text(json.dumps(f,indent=2)+'\n')
text=(base/'verify.wl').read_text().split('BLFamilyDefine[')[0].replace('(l1+l2)^2-1','(l1+l2)^2-5/4')
text+='BLFamilyDefine[thresholdsplit,26/7,prop,{l1,l2},{Q},{},{Q^2->2},top,{},top];\n'
text+='start=AbsoluteTime[];mi=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];Export["result.json",<|"count"->Length[mi],"masters"->ToString[mi,InputForm],"seconds"->(AbsoluteTime[]-start)|>,"RawJSON"];Print["RESULT=",mi];Exit[];\n'
(out/'count.wl').write_text(text)
