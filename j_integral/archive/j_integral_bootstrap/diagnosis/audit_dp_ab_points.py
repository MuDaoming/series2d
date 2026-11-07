"""Audit fitted, per-sector checked, final checked and unnecessary DP samples."""
import argparse,json,re,collections
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--connection',required=True);p.add_argument('--numeric-files',required=True);p.add_argument('--events',required=True);p.add_argument('--output',required=True);a=p.parse_args()
d=json.loads(Path(a.connection).read_text());prime=d['prime']
from fractions import Fraction
q=Fraction(d['t0']);center=q.numerator*pow(q.denominator,-1,prime)%prime
training=set().union(*(set(c['training_points']) for c in d['checks']))
sector_checks=set().union(*(set(c['holdout_points']) for c in d['checks']))
final=set(d['global_check']['holdout_points']);normal=training|sector_checks|final
normal_t={(x+center)%prime for x in normal};counts=collections.Counter();points=set()
for line in Path(a.numeric_files).read_text().splitlines():
 m=re.search(r'point(\d+)\.json$',line)
 if not m:continue
 t=int(m[1]);points.add(t);kind='restricted' if '/sectors/' in line else 'full'
 counts['normal' if t in normal_t else 'prefix',kind]+=1
records=[json.loads(s) for s in Path(a.events).read_text().splitlines()];events=collections.Counter()
for r in records:
 kind='prefix' if r.get('task','').startswith('__prefix_de_holdouts__') else 'final' if r.get('task','').startswith('__complete_de_holdouts__') else 'sector'
 events[kind,r['route'],r.get('stage','evaluate')]+=1
result=dict(training_points_tau=sorted(training),sector_holdout_points_tau=sorted(sector_checks),final_holdout_points_tau=sorted(final),training_distinct=len(training),training_and_sector_checks_distinct=len(training|sector_checks),normal_total_distinct=len(normal),normal_points_t=sorted(normal_t),normal_full_successes=counts['normal','full'],normal_restricted_successes=counts['normal','restricted'],prefix_full_successes=counts['prefix','full'],prefix_restricted_successes=counts['prefix','restricted'],all_successes=sum(counts.values()),all_distinct_t=len(points),prefix_distinct_t=len(points-normal_t),events={'/'.join(k):v for k,v in events.items()},cause='run_upstream per-sector partial solve invokes fresh prefix holdouts after each of148 incomplete prefixes;20 normal points vs296 extra prefix points')
assert len(normal)==20 and len(training)==16 and len(training|sector_checks)==18
assert sum(counts.values())==422 and len(points-normal_t)==2*(len(d['checks'])-1)
assert counts['normal','full']==20 and counts['normal','restricted']==106
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if not isinstance(v,(dict,list))},indent=2))
