import json,random,subprocess,time,resource,concurrent.futures,sys
from pathlib import Path
base=Path(__file__).resolve().parent;root=base.parents[3];separated='--separated' in sys.argv;out=base/('j_direct_separated' if separated else 'j_direct');out.mkdir(exist_ok=True);(out/'points').mkdir(exist_ok=True);(out/'logs').mkdir(exist_ok=True);p=1000003;center=2*pow(5,-1,p)%p
source=json.load(open(root/('j_integral/issues/mechanisms/branch_factorization/e6_full/request_blade_basis.json' if separated else 'j_integral/tests/e6/input/family_request.json')));source['request']['targets']=json.load(open(base/'input/common_targets.json'));source['request']['derivatives']=[]
basis=source['request']['basis'];names=['J{'+','.join(map(str,nu))+'}' for nu in source['request']['targets']]
rng=random.Random(20261115);points=[]
while len(points)<520:
 t=rng.randrange(2,p)
 if t not in points:points.append(t)
def sample(t):
 path=out/'points'/f'point{t}.json'
 if path.exists():return path
 req=json.loads(json.dumps(source));req['point']['t']=str(t);inp=out/'points'/f'input{t}.json';inp.write_text(json.dumps(req,separators=(',',':')))
 started=time.monotonic()
 with open(out/'logs'/f'point{t}.log','w') as log:r=subprocess.run([str(root/'j_integral/src/IBP/build/j-ibp'),'--input',str(inp),'--output',str(path)],stdout=log,stderr=log)
 if r.returncode:raise RuntimeError(f'IBP failed t={t}, {r.returncode}')
 result=json.load(open(path));assert result['basis']==basis and result['status']=='complete'
 record={'operation':'numeric_IBP','t':t,'wall_seconds':time.monotonic()-started,'diagnostics':result['diagnostics']}
 with open(out/'events.jsonl','a') as f:f.write(json.dumps(record)+'\n')
 return path
for train in ([2,4,6,8,10,12,16,24,32] if separated else [32,48,64,80,96,112,128,144,160,192,224,256,320,384,448,512]):
 count=train+4
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:files=list(pool.map(sample,points[:count]))
 config={'prime':p,'xs':[(t-center)%p for t in points[:count]],'files':[str(f) for f in files],'targets':names,'basis':basis,'train':train};(out/'fit_input.json').write_text(json.dumps(config,separators=(',',':')))
 subprocess.run([str(base/'reconstruct'),str(out/'fit_input.json'),str(out/'fit_output.json')],check=True)
 fit=json.load(open(out/'fit_output.json'));print({k:v for k,v in fit.items() if k!='matrix'},flush=True)
 if fit['passed']:
  stats={k:v for k,v in fit.items() if k!='matrix'};stats.update(actual_IBP_calls=count,prime=p,targets=len(names),basis_size=len(basis),variable='tau=t-2/5',anchor_cost_included=True)
  (out/'statistics.json').write_text(json.dumps(stats,indent=2)+'\n');(out/'matrix.json').write_text(json.dumps(fit['matrix'],separators=(',',':'))+'\n');break
else:raise RuntimeError('reconstruction did not close')
