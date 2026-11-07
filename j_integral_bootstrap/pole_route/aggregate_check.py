"""Merge supplemental coefficients and independently check every available nonblind value."""
import json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from rational_learning import make_learner
VERSION=sys.argv[1] if len(sys.argv)>1 else 'v2'
R=ROOT/'results';O=R/'dp_pole_route';L=make_learner(1000003,200001);start=time.monotonic()
sources={'de':['shared_de','shared_de_pending_v2','preconditioned_de','preconditioned_de_all','preconditioned_de_three','native_de43'],
 'reductions':['shared_reductions','shared_reductions_v2','shared_reductions_v3','preconditioned_reductions','preconditioned_reductions_all','preconditioned_reductions_three','native_reductions43']}
connection=json.load(open(R/'dp_complete_run/connection.json'));out={'prime':1000003,'t0':'2/5','masters':connection['masters'],'de':{},'reductions':{}}
report={'sources':{},'de':{},'reductions':{}};points={}
for dr in ['dp_full','dp_extra_de','dp_extra_reductions','dp_more_full']:
 for fn in (R/dr/'points').glob('point*.json'):
  d=json.load(open(fn));t=int(d['t_point']);v=points.setdefault(t,{'DE_values':{},'reduction_values':{}})
  for sec in v:v[sec].update(d[sec])
for section,files in sources.items():
 for fn in files:
  path=O/f'{fn}.json'
  if not path.exists():continue
  d=json.load(open(path));report['sources'][fn]={'coefficients':d.get('new_count'),'seconds':d.get('seconds'),'internal_check_points':d.get('internal_holdout',3 if fn.endswith('three') or fn.startswith('shared') else 8)}
  for row,vs in d[section].items():
   if section=='de':vs=vs['rhs']
   dest=out[section].setdefault(row,{'rhs':{}} if section=='de' else {})
   if section=='de':dest=dest['rhs']
   for col,co in vs.items():
    if col in dest and co!=dest[col]:raise RuntimeError(('duplicate mismatch',row,col,fn))
    dest[col]=co
 sec='DE_values' if section=='de' else 'reduction_values';total=checks=0;failed=[];hist={}
 for row,vs in out[section].items():
  if section=='de':vs=vs['rhs']
  ts=[t for t in points if row in points[t][sec]]
  for col,c in vs.items():
   total+=1;hist[len(ts)]=hist.get(len(ts),0)+1;n=L.poly(c['numerator_tau']);den=L.poly(c['denominator_tau'])
   for t in ts:
    x=(t-200001)%1000003;v=points[t][sec][row].get(col,0);checks+=1
    if not int(den(x)) or int(n(x))!=v*int(den(x))%1000003:failed.append([row,col,t])
 report[section]={'coefficients':total,'pointwise_checks':checks,'points_per_coefficient':hist,'failures':failed}
 if failed:raise RuntimeError(report[section])
for section,sec,basefiles in [('de','DE_values',['dp_complete_run/connection.json','dp_driven/partial_33.json','dp_extra_de/partial41.json','dp_cuts/reconstruction.json','dp_exceptional/diagonal33.json']),('reductions','reduction_values',['dp_reduction_run/reduction.json','dp_full/reductions_partial33.json','dp_extra_reductions/partial41.json'])]:
 known={};allcols={};maxnum=maxden=0
 for fn in basefiles:
  for row,vs in json.load(open(R/fn))[section].items():
   if section=='de':vs=vs['rhs']
   known.setdefault(row,{}).update(vs)
 for row,vs in out[section].items():
  if section=='de':vs=vs['rhs']
  known.setdefault(row,{}).update(vs)
 for t,pv in points.items():
  for row,vs in pv[sec].items():allcols.setdefault(row,set()).update(c for c,v in vs.items() if v)
 missing={row:sorted(cols-set(known.get(row,{}))) for row,cols in allcols.items()}
 for vs in known.values():
  for co in vs.values():maxnum=max(maxnum,len(co['numerator_tau'])-1);maxden=max(maxden,len(co['denominator_tau'])-1)
 report[section]['coverage']={'total_nonzero_columns':sum(map(len,allcols.values())),'missing_coefficients':sum(map(len,missing.values())),'complete_rows':sum(not v for v in missing.values()),'total_rows':len(missing),'highest_numerator_degree':maxnum,'highest_denominator_degree':maxden}
 (O/f'{section}_remaining_{VERSION}.json').write_text(json.dumps({r:v for r,v in missing.items() if v}))
report['seconds']=time.monotonic()-start;report['blind_verification_read']=False

for name,data in [(f'supplemental_{VERSION}.json',out),(f'supplemental_{VERSION}_checks.json',report)]:
 path=O/name
 if path.exists():raise RuntimeError('immutable output exists: '+str(path))
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data));tmp.replace(path)
print(json.dumps(report))
