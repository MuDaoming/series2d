import itertools,json,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[4];out=Path(__file__).resolve().parent
f=json.load(open(root/'j_integral/issues/mechanisms/six_line_null_legs/both_null/request.json'))
# Stable dot<=2 input for all branch orders.
targets=[]
for bits in itertools.product([0,1],repeat=5):
 if not any(bits[:2]) or not any(bits[2:]):continue
 nu=list(bits)+[1];active=[i for i,v in enumerate(nu) if v]
 for order in range(3):
  for dots in itertools.combinations_with_replacement(active,order):
   q=nu.copy()
   for i in dots:q[i]+=1
   targets.append(q)
f['request']['targets']=targets;original=f['j']['branches'];summary=[]
for perm in itertools.permutations(range(3)):
 d=out/(''.join(map(str,perm)));d.mkdir(exist_ok=True);q=json.loads(json.dumps(f));q['j']['branches']=[perm[b] for b in original];(d/'request.json').write_text(json.dumps(q,indent=2)+'\n')
 with open(d/'run.log','w') as log:p=subprocess.run([str(root/'j_integral/src/IBP/build/j-ibp'),'--input',str(d/'request.json'),'--output',str(d/'j_result.json')],stdout=log,stderr=log)
 if not (d/'j_result.json').exists():raise RuntimeError((perm,p.returncode))
 r=json.load(open(d/'j_result.json'));top=[n for n in r['basis'] if n.startswith('J') and all(int(x)>0 for x in n.split('{')[1][:-1].split(','))]
 row={'old_to_new_branches':perm,'J_top':len(top),'total_J_system':len(r['basis']),'unclosed':r['unclosed_derivatives'],'returncode':p.returncode};summary.append(row);print(row,flush=True)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
