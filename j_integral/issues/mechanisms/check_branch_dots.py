"""Exact parameter-weight identity checked against every coefficient of numerical reductions."""
import json
from pathlib import Path
root=Path(__file__).resolve().parent;rows=[]
for d in sorted((root/'branch_choice').iterdir()):
 if not d.is_dir() or not (d/'request.json').exists():continue
 f=json.load(open(d/'request.json'));r=json.load(open(d/'j_result.json'));p=r['point']['prime'];t=int(r['point']['t']);br=f['j']['branches'];n=len(br);acc={}
 for i,b in enumerate(br):
  if b==2:continue
  nu=[1]*n;nu[i]+=1;name='J{'+','.join(map(str,nu))+'}';factor=(1-t) if b==0 else -t
  for key,c in r['reductions'][name].items():acc[key]=(acc.get(key,0)+factor*c)%p
 assert not any(acc.values()),(d,acc)
 rows.append({'case':d.name,'branch_groups':[ [i+1 for i,b in enumerate(br) if b==j] for j in range(3)],'coefficient_identity_checked':True})
(root/'branch_dot_checks.json').write_text(json.dumps(rows,indent=2)+'\n');print('All',len(rows),'branch dot identities PASS')
