"""Validate all 53 DP rows and learned-name-only e6 task closures."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from restricted_numeric import SectorRestrictedNumericBackend,closure_tasks,indices
from upstream_j import ROOT
OUT=ROOT/'results/dp_official/restricted_backend';reports=[]
def check(backend,req,points):
 session=backend.prepare(req)
 for t in [17,18]:
  full=json.loads((Path(points)/f'point{t}.json').read_text());x=(t-backend.center)%backend.prime;r=session.evaluate(x)
  expected=[[full['DE_values'][target['master']].get(m,0) for m in req['masters']] for target in req['targets']]
  assert r['values']==expected
  reports.append(dict(task=req['task'],t=t,rows=len(req['targets']),columns=len(req['masters']),compared=len(req['targets'])*len(req['masters']),event=backend.events[-1]))
 session.close()
cfg=ROOT.parent/'j_integral_global/examples';names=json.loads((cfg/'dp496_basis.json').read_text());support={0,1,2,3,6,7};retained=[m for m in names if all(not v or i in support for i,v in enumerate(indices(m)[1]))];assert len(retained)==53
b=SectorRestrictedNumericBackend(cfg/'dp_physical.json',cfg/'dp496_basis.json',cfg/'dp_targets_dot2.json',ROOT.parent/'j_integral_global/build/global_point',OUT/'dp_all53',point_directories=[ROOT/'results/dp_full/points'])
req=dict(mode='de',prime=b.prime,masters=names,task='all53',targets=[dict(kind='derivative',variable='tau',master=m) for m in retained]);check(b,req,ROOT/'results/dp_full/points');assert all(r['event']['route']=='restricted' for r in reports)
case=ROOT.parent/'j_integral_sixpoint';bundle=json.loads((case/'results/upstream/bundle.json').read_text());tasks=closure_tasks(bundle);names=json.loads((case/'inputs/basis.json').read_text());family=json.loads((case/'inputs/family.json').read_text());sectors={s['id']:s for s in bundle['definition']['sectors']}
choices=[]
for task,closure in tasks.items():
 active={i for m in closure for i,v in enumerate(indices(m)[1]) if v}
 if 4<=len(active)<len(family['propagators']) and len(closure)>3:choices.append((len(closure),task))
for _,task in [min(choices),max(choices)]:
 b=SectorRestrictedNumericBackend(case/'inputs/family.json',case/'inputs/basis.json',ROOT/'official_route/cases/sixpoint_basis_targets.json',ROOT.parent/'j_integral_global/build/global_point',OUT/('e6closure_'+task),point_directories=[case/'results/samples/points'],closure_by_task=tasks)
 req=dict(mode='de',prime=b.prime,masters=list(reversed(names)),task=task,targets=[dict(kind='derivative',variable='tau',master=m) for m in sectors[task]['masters']]);check(b,req,case/'results/samples/points')
summary=dict(status='passed',comparisons=sum(r['compared'] for r in reports),reports=reports);(OUT/'closure_regression.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
