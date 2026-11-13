"""Verify compressed evidence without rerunning IBP, using only Python stdlib."""
import gzip,json,subprocess,sys
from pathlib import Path
b=Path(__file__).resolve().parent
load=lambda name:json.loads(gzip.decompress((b/'reference'/name).read_bytes()))
def evaluate(poly,x,p):
 v=0
 for c in reversed(poly):v=(v*x+c)%p
 return v
for tag,source in [('original','j_original_common_prime.json'),('separated','j_separated_common_prime.json')]:
 direct=load(f'j_{tag}_direct_matrix.json.gz');ab=load(f'j_{tag}_AB_result.json.gz');source_basis=json.load(open(b/source))['basis'] if (b/source).exists() else load(source+'.gz')['basis'];ab_basis=ab['context']['masters'];columns={name:i for i,name in enumerate(ab_basis)}
 mapped={(i,columns[source_basis[j]]):pair for i,j,pair in direct['entries']};other={(i,j):pair for i,j,pair in ab['matrix']['entries']};assert mapped==other
 holdouts=load(f'j_{tag}_holdouts.json.gz');p=1000003
 for point in holdouts:
  predicted={}
  for i,j,(num,den) in direct['entries']:
   value=evaluate(num,point['x'],p)*pow(evaluate(den,point['x'],p),-1,p)%p
   if value:predicted.setdefault(i,{})[source_basis[j]]=value
  targets=json.load(open(b/'input/common_targets.json'))
  for i,nu in enumerate(targets):assert predicted.get(i,{})==point['reductions']['J{'+','.join(map(str,nu))+'}']
 print(tag,'AB/direct exact rational entries and',len(holdouts),'independent holdouts PASS')
for tag,expected in [('original',148),('separated',14)]:
 events=[json.loads(x) for x in gzip.decompress((b/'reference'/f'j_{tag}_direct_events.jsonl.gz').read_bytes()).decode().splitlines()];assert len(events)==len({e['t'] for e in events})==expected
for tag,expected in [('original',17),('separated',10)]:
 events=[json.loads(x) for x in gzip.decompress((b/'reference'/f'j_{tag}_AB_events.jsonl.gz').read_bytes()).decode().splitlines()];calls=[e for e in events if e['operation']=='IBP-complete'];assert len(calls)==len({e['t'] for e in calls})==expected
subprocess.run([sys.executable,str(b/'verify_fi_j.py')],check=True)
print('Actual IBP counts and full FI/J numeric comparison PASS')
