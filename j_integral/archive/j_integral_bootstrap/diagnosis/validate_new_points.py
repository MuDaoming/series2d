"""Validate frozen canonical known DE/R, including zeros, at new learning points."""
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from reference_space import RatFunc
from fractions import Fraction
out=ROOT/'results/dp_diagnosis';results=ROOT/'results';start=time.perf_counter()
datasets={}
for label,file,knownkey,section in [('DE','dp_complete_run/connection.json','known_rows','de'),('R','dp_reduction_run/reduction.json','known_targets','reductions')]:
 raw=(results/file).read_bytes();(out/('validation_snapshot_'+label+'.json')).write_bytes(raw);d=json.loads(raw);p=d['prime'];f=Fraction(d['t0']);center=f.numerator*pow(f.denominator,-1,p)%p
 rows={}
 for m in d[knownkey]:
  source=d[section].get(m)
  if source is None:
   if label!='R' or m not in d['masters']:raise ValueError('nonmaster missing reduction '+m)
   source={m:{'numerator_tau':[1],'denominator_tau':[1]}}
  source=source['rhs'] if label=='DE' else source
  rows[m]={j:RatFunc.make(v['numerator_tau'],v['denominator_tau'],p) for j,v in source.items()}
 blocks={m:{'index':i,'top':b['top'],'route':b['route']} for i,b in enumerate(d.get('blocks',[])) for m in b.get('top',[])}
 datasets[label]={'rows':rows,'blocks':blocks,'snapshot_sha256':hashlib.sha256(raw).hexdigest(),'snapshot_path':file,'checks':0,'implicit_zero_checks':0,'mismatches':[],'first_failed_entries':{},'checked_points':[],'known_objects':len(rows)}
paths=sorted((results/'dp_more_full/points').glob('point*.json'),key=lambda x:int(json.loads(x.read_text())['t_point']))
for path in paths:
 point=json.loads(path.read_text());t=int(point['t_point']);x=(t-center)%p
 if point['unclosed_derivatives']:raise ValueError('unclosed point')
 for label,dataset in datasets.items():
  expected=point['DE_values'] if label=='DE' else point['reduction_values'];dataset['checked_points'].append(t)
  for m,row in dataset['rows'].items():
   if m not in expected:raise ValueError('missing actual object '+m)
   truth=expected[m];dataset['checks']+=len(point['masters_at_point']);dataset['implicit_zero_checks']+=len(point['masters_at_point'])-len(row)
   for j in set(row)|set(truth):
    try:actual=row[j].evaluate(x) if j in row else 0
    except ZeroDivisionError:actual='pole'
    target=truth.get(j,0)
    if actual!=target:
     failure={'point':t,'row':m,'column':j,'checkpoint_value':actual,'native_value':target,'block':dataset['blocks'].get(m),'numerator_degree':len(row[j].num)-1 if j in row else None,'denominator_degree':len(row[j].den)-1 if j in row else None}
     dataset['mismatches'].append(failure);dataset['first_failed_entries'].setdefault(m+' -> '+j,failure)
 report={label:{k:v for k,v in d.items() if k not in ('rows','blocks')} for label,d in datasets.items()};report['seconds']=time.perf_counter()-start
 (out/'new_points_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps({'point':t,'DE_mismatches':len(datasets['DE']['mismatches']),'R_mismatches':len(datasets['R']['mismatches'])}),flush=True)
