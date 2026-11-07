"""Actual sector-restricted worker vs saved full numerical outputs, all columns."""
import sys,json,itertools
from pathlib import Path
from time import perf_counter
sys.path.insert(0,str(Path(__file__).resolve().parent))
from restricted_numeric import SectorRestrictedNumericBackend,indices
from upstream_j import ROOT
OUT=ROOT/'results/dp_official/restricted_backend';OUT.mkdir(parents=True,exist_ok=True)

def run_case(name,family,basis,targets,points):
    names=json.loads(Path(basis).read_text());cfg=json.loads(Path(family).read_text());n=len(cfg['propagators'])
    backend=SectorRestrictedNumericBackend(family,basis,targets,ROOT.parent/'j_integral_global/build/global_point',OUT/name,point_directories=[points],gauge=3)
    groups={}
    for m in names:
        _,nu=indices(m);support=tuple(i for i,v in enumerate(nu) if v)
        if {cfg['propagators'][i]['branch'] for i in support}=={0,1,2}:groups.setdefault(support,[]).append(m)
    selected=[]
    for size in range(3,n+1):
        candidates=[s for s in groups if len(s)==size]
        if candidates:selected.append(sorted(candidates,key=lambda s:(-len(groups[s]),s))[0])
    if name=='dp':
        known=(0,1,2,3,6,7)
        if known in groups and known not in selected:selected.append(known)
    # Test additional high supports as well as the widest/largest blocks.
    if name=='e6':
        selected.extend(s for s in sorted(groups) if len(s)==5 and s not in selected)
    reports=[]
    for support in selected:
        req=dict(mode='de',prime=cfg['prime'],masters=list(reversed(names)),targets=[dict(kind='derivative',variable='tau',master=m) for m in groups[support]],task='physical-support-'+str(support))
        session=backend.prepare(req)
        for t in [17,18]:
            full=json.loads((Path(points)/f'point{t}.json').read_text());x=(t-backend.center)%backend.prime;start=perf_counter();result=session.evaluate(x)
            shift=3*(pow(t,-1,backend.prime)+pow(t-1,-1,backend.prime))%backend.prime
            expected=[[(full['DE_values'][target['master']].get(m,0)+(shift if m==target['master'] else 0))%backend.prime for m in req['masters']] for target in req['targets']]
            assert result['values']==expected,(name,support,t)
            reports.append(dict(case=name,support=support,t=t,target_rows=len(req['targets']),columns=len(names),compared=len(req['targets'])*len(names),seconds=perf_counter()-start,event=backend.events[-1]))
        session.close()
    # Unsupported branch restriction must explicitly route to full, not shrink.
    m=next((m for m in names if len({cfg['propagators'][i]['branch'] for i,v in enumerate(indices(m)[1]) if v})<3),None)
    if m:
        backend.prepare(dict(mode='de',prime=cfg['prime'],masters=names,targets=[dict(kind='derivative',variable='tau',master=m)],task='fallback-branch-check'))
        assert backend.events[-1]['route']=='full'
    return reports

reports=[]
reports+=run_case('dp',ROOT.parent/'j_integral_global/examples/dp_physical.json',ROOT.parent/'j_integral_global/examples/dp496_basis.json',ROOT.parent/'j_integral_global/examples/dp_targets_dot2.json',ROOT/'results/dp_full/points')
reports+=run_case('e6',ROOT.parent/'j_integral_sixpoint/inputs/family.json',ROOT.parent/'j_integral_sixpoint/inputs/basis.json',ROOT/'official_route/cases/sixpoint_basis_targets.json',ROOT.parent/'j_integral_sixpoint/results/samples/points')
summary=dict(status='passed',comparisons=sum(r['compared'] for r in reports),reports=reports)
(OUT/'regression.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
