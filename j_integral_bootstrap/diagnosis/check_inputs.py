import sys,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from reference_space import RatFunc
from fractions import Fraction
r=ROOT/'results';cfg=json.loads((r/'dp_full/input/family.json').read_text());p=cfg['prime'];f=Fraction(cfg['t0']);center=f.numerator*pow(f.denominator,-1,p)%p
sources={}
for name,path in [('known',r/'dp_diagnosis/snapshot/connection.json'),('cut',r/'dp_cuts/reconstruction.json'),('exceptional',r/'dp_exceptional/diagonal33.json')]:
 d=json.loads(path.read_text());sources[name]={(m,j):RatFunc.make(v['numerator_tau'],v['denominator_tau'],p) for m,row in d['de'].items() for j,v in row['rhs'].items()}
report={n:{'checks':0,'mismatches':[],'pole_points':0,'entries':len(v)} for n,v in sources.items()};begin=time.perf_counter()
for directory in ['dp_full','dp_extra_de','dp_exceptional']:
 for path in sorted((r/directory/'points').glob('point*.json')):
  d=json.loads(path.read_text());x=(int(d['t_point'])-center)%p
  for name,values in sources.items():
   for (m,j),rf in values.items():
    if m not in d['DE_values']:continue
    try:actual=rf.evaluate(x)
    except ZeroDivisionError:report[name]['pole_points']+=1;continue
    expected=d['DE_values'][m].get(j,0);report[name]['checks']+=1
    if actual!=expected:report[name]['mismatches'].append([directory,path.name,m,j,actual,expected])
report['seconds']=time.perf_counter()-begin
(r/'dp_diagnosis/input_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
