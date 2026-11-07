import json
from pathlib import Path
r=Path(__file__).resolve().parents[1]/'results';d=json.loads((r/'dp_diagnosis/snapshot/connection.json').read_text());known=set(d['known_rows']);n=0;bad=[]
for directory in ['dp_full','dp_extra_de','dp_exceptional']:
 for path in (r/directory/'points').glob('point*.json'):
  point=json.loads(path.read_text())
  for m,row in point['DE_values'].items():
   if m not in known:continue
   rhs=d['de'][m]['rhs'];n+=len(d['masters'])-len(rhs)
   for j,v in row.items():
    if v and j not in rhs:bad.append([directory,path.name,m,j,v])
out={'implicit_zero_checks':n,'mismatches':bad};(r/'dp_diagnosis/known_zeros.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
