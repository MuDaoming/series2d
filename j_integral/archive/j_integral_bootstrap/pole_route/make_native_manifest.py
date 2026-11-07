import json,sys,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';O=R/'dp_pole_route'
kind=sys.argv[1];extra='dp_extra_de' if kind=='de' else 'dp_extra_reductions';points={}
for dr in ['dp_full',extra,'dp_more_full']:
 for path in sorted((R/dr/'points').glob('point*.json')):points[str(path.resolve())]=True
sources={'de':['dp_complete_run/connection.json','dp_cuts/reconstruction.json','dp_exceptional/diagonal33.json','dp_driven/partial_33.json','dp_extra_de/partial41.json'],
'reductions':['dp_reduction_run/reduction.json','dp_full/reductions_partial33.json','dp_extra_reductions/partial41.json']}[kind]
for stem in (['shared_de','shared_de_pending_v2','preconditioned_de','preconditioned_de_all','preconditioned_de_three'] if kind=='de' else ['shared_reductions','shared_reductions_v2','shared_reductions_v3','preconditioned_reductions','preconditioned_reductions_all','preconditioned_reductions_three']):
 p=O/f'{stem}.json'
 if p.exists():sources.append(str(p.relative_to(R)))
for p in O.glob(f'native_{kind}*.json'):
 if re.fullmatch(r'native_'+kind+r'\d+',p.stem):sources.append(str(p.relative_to(R)))
manifest={'prime':1000003,'t0':'2/5','masters':json.load(open(R/'dp_complete_run/connection.json'))['masters'],'holdout':3,'section':'DE_values' if kind=='de' else 'reduction_values','points':list(points),'known':[str((R/s).resolve()) for s in sources]}
p=O/f'native_{kind}{len(points)}_manifest.json';p.write_text(json.dumps(manifest));print(p,len(points))
