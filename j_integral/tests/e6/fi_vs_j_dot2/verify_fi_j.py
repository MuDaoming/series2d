import json
from fractions import Fraction
from pathlib import Path
b=Path(__file__).resolve().parent
import gzip
def saved(path,reference):
 return json.load(open(path)) if path.exists() else json.loads(gzip.decompress((b/'reference'/reference).read_bytes()))
fi=saved(b/'fi/run_rank1/reductions.json','fi_reductions.json.gz');j=saved(b/'j_separated_common_prime.json','j_separated_common_prime.json.gz');p=fi['prime'];assert j['point']['prime']==p;t=int(j['point']['t']);d=Fraction(j['point']['dimension']);h=d.numerator*pow(2*d.denominator,-1,p)%p
alias=(1,0,1,0,0,1,0);canon=(0,1,0,0,1,0,1)
def indices(label):return tuple(map(int,label.split('{')[1][:-1].split(',')))
def ratio(nu):
 A=sum(nu[:2]);C=sum(nu[5:]);B=1
 for k in range(1,A+C-1):B=B*(2*h-2-k)%p
 for n in [A,C]:
  for k in range(1,n):B=B*pow((h-1-k)%p,-1,p)%p
 return pow(t,1-C,p)*pow((1-t)%p,1-A,p)*pow(B,-1,p)%p
masters=[tuple(v[:7]) for v in fi['masters']];common=set(map(tuple,json.load(open(b/'input/common_targets.json'))));checks=0;zeros=0
for nu,values in zip(fi['targets'],fi['values']):
 nu=tuple(nu)
 if nu not in common:assert not any(values);zeros+=1;continue
 rhs={}
 for master,c in zip(masters,values):
  if master==alias:master=canon
  if c:rhs[master]=(rhs.get(master,0)+c)%p
 rhs={k:v for k,v in rhs.items() if v}
 label='J{'+','.join(map(str,nu))+'}';jr={}
 for name,c in j['reductions'][label].items():
  master=indices(name);v=c*ratio(master)*pow(ratio(nu),-1,p)%p
  if v:jr[master]=v
 assert jr==rhs,(nu,jr,rhs)
 checks+=1
report={'prime':p,'common_three_branch_targets_checked':checks,'additional_zero_targets_checked':zeros,'raw_Blade_master_count':len(masters),'after_explicit_duplicate_merge':len(set(canon if m==alias else m for m in masters)),'all_FI_vs_J_reductions_equal':True}
(b/'verification_fi_j.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
