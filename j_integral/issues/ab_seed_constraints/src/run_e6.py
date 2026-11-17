"""Small reproducible AB-before-z-IBP tests; all heavy algebra stays in existing kernels."""
import argparse,gzip,json,os,subprocess,time
from pathlib import Path
here=Path(__file__).resolve().parents[1];root=here.parents[2]
parser=argparse.ArgumentParser();parser.add_argument('--case',choices=['original63','original_small','separated63'],required=True);args=parser.parse_args()
base=root/'j_integral/tests/e6/fi_vs_j_dot2'
request=json.loads((base/'input'/('j_separated_common_prime.json' if args.case=='separated63' else 'j_original_common_prime.json')).read_text())
if args.case=='original_small':request['point'].update(prime=1000003,t='23')
work=here/'build'/args.case;work.mkdir(parents=True,exist_ok=True);(work/'request.json').write_text(json.dumps(request)+'\n')
exe=here/'build/experiment';times={}
def run(name,command):
 started=time.perf_counter()
 with (work/(name+'.log')).open('w') as log:result=subprocess.run(list(map(str,command)),stdout=log,stderr=log)
 times[name]=time.perf_counter()-started
 if result.returncode:raise RuntimeError(f'{name} failed: {work/(name+".log")}')
 print(args.case,name,round(times[name],3),'s',flush=True)
run('snapshot',[exe,'snapshot',work/'request.json',work/'snapshot.json'])
run('module',['python3',here/'src/build_module.py',work/'snapshot.json',work/'module','--physical-z0'])
run('reference',[exe,'reduce',work/'snapshot.json',work/'reference.json','reference'])
run('image_ab',[exe,'image',work/'snapshot.json',work/'image_ab.json',work/'module/ab_seeds.json'])
run('image_trivial',[exe,'image',work/'snapshot.json',work/'image_trivial.json',work/'module/trivial_seeds.json'])
module=json.loads((work/'module/module_summary.json').read_text());ref=json.loads((work/'reference.json').read_text())
summary={'case':args.case,'point':request['point'],'targets':len(request['request']['targets']),'J_basis_size':len(request['request']['basis']),'primitive_parameters_before':module['ambient_parameters'],'primitive_parameters_AB':module['AB_parameters'],'physical_rows_excluded':len(module['unrepresented_physical_rows']),'timings_wall':times,'reference':{k:v for k,v in ref.items() if k not in ['basis','reductions','unresolved']},'candidates':{}}
for name in ['image_ab','image_trivial']:
 data=json.loads((work/(name+'.json')).read_text());equal=data['reductions']==ref['reductions'];summary['candidates'][name]={k:v for k,v in data.items() if k not in ['reductions','unresolved']};summary['candidates'][name].update(all_target_coefficients_equal=equal,unresolved_count=len(data['unresolved']))
 if not equal or data['unresolved'] or data['unrepresented_relations'] or not data['basis_independent']:raise RuntimeError('candidate failed exact reduction verification')
if module['unrepresented_physical_rows']:raise RuntimeError('physical densities excluded by module')
if args.case=='original63':
 expected=json.load(gzip.open(base/'reference/j_original_common_prime.json.gz','rt'))
 summary['baseline_matches_saved_public_result']=all(expected['reductions'][k]==v for k,v in ref['reductions'].items())
 if not summary['baseline_matches_saved_public_result']:raise RuntimeError('baseline differs from saved public result')
for candidate in ['image_ab','image_trivial']:
 with gzip.open(here/'results'/(args.case+'_'+candidate+'.json.gz'),'wb') as f:f.write((work/(candidate+'.json')).read_bytes())
summary['passed']=True
(here/'results'/(args.case+'.json')).write_text(json.dumps(summary,indent=2)+'\n')
(here/'results'/(args.case+'_module.json')).write_text(json.dumps(module,indent=2)+'\n')
with gzip.open(here/'results'/(args.case+'_reductions.json.gz'),'wt') as f:json.dump(ref['reductions'],f,separators=(',',':'))
print(json.dumps(summary,indent=2),flush=True)
