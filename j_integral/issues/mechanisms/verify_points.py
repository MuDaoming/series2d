"""Recheck top J counts and DE closure at a second prime and t point."""
import json,subprocess,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[3];base=Path(__file__).resolve().parent
paths=sorted(set(base.glob('*/*/request.json'))|set(base.glob('*/request.json')));rows=[]
def top_count(r):return sum(n.startswith('J') and all(int(q)>0 for q in n.split('{')[1][:-1].split(',')) for n in r['basis'])
for path in paths:
 old=json.load(open(path.with_name('j_result.json')));f=json.load(open(path));f['point']['prime']=1000033;f['point']['t']='19'
 with tempfile.TemporaryDirectory(prefix='j_mechanism_') as tmp:
  tmp=Path(tmp);(tmp/'request.json').write_text(json.dumps(f))
  p=subprocess.run([str(root/'j_integral/src/IBP/build/j-ibp'),'--input',str(tmp/'request.json'),'--output',str(tmp/'result.json')],capture_output=True,text=True)
  assert p.returncode==0,(str(path),p.stderr)
  new=json.load(open(tmp/'result.json'));assert not new['unclosed_derivatives'];assert top_count(new)==top_count(old),(path,top_count(new),top_count(old))
  row={'case':str(path.parent.relative_to(base)),'prime':1000033,'t':19,'top_J_count':top_count(new),'total_J_system':len(new['basis']),'top_count_stable':True,'DE_closed':True};rows.append(row);print(row,flush=True)
(base/'second_point_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
