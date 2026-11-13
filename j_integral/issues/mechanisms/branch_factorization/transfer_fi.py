"""Turn separated J reductions into FI reductions and check t independence."""
import gzip,json
from pathlib import Path
from fractions import Fraction
base=Path(__file__).resolve().parent;root=base.parents[3]
def nu(label):return list(map(int,label.split('{')[1][:-1].split(',')))
def convert(native):
 p=native['point']['prime'];t=int(native['point']['t']);d=Fraction(native['point']['dimension']);h=d.numerator*pow(2*d.denominator,-1,p)%p
 def ratio(label):
  v=nu(label);A=sum(v[:2]);C=sum(v[5:]);b=1
  for k in range(1,A+C-1):b=b*(2*h-2-k)%p
  for length in [A,C]:
   for k in range(1,length):b=b*pow((h-1-k)%p,-1,p)%p
  return pow(t,1-C,p)*pow((1-t)%p,1-A,p)*pow(b,-1,p)%p
 out={}
 for label,row in native['reductions'].items():
  if not label.startswith('J'):continue
  factor=pow(ratio(label),-1,p);result={}
  for master,c in row.items():
   assert master.startswith('J'),master
   value=c*ratio(master)*factor%p
   if value:result['I'+master[1:]]=value
  out['I'+label[1:]]=result
 return out
first=json.load(open(base/'e6_full/j_blade_basis.json'));other=json.load(open(base/'e6_full/j_blade_basis_t23.json'));a=convert(first);b=convert(other);assert a==b
report={'prime':first['point']['prime'],'dimension':first['point']['dimension'],'t_checks':[17,23],'targets_checked':len(a),'all_coefficients_t_independent':True}
path=root/'j_integral/issues/mechanisms/six_line_null_legs/both_null/reduction_verification.json'
if path.exists():
 blade=json.load(open(path));active=[0,1,2,3,4,6];mask=[1,1,1,1,1,0,1];masters=blade['master_indices'];checks=[]
 for record in blade['records']:
  assert not record['unresolved'];v=[0]*7
  for i,k in zip(active,record['indices']):v[i]=k
  label='I{'+','.join(map(str,v))+'}';top={}
  for m,c in a[label].items():
   q=nu(m)
   if [int(x>0) for x in q]==mask:top[tuple(q[i] for i in active)]=c
  expected={tuple(v):c for v,c in zip(masters,record['coefficients']) if c}
  assert top==expected,(label,top,expected)
  checks.append(record['indices'])
 report['independent_Blade_sixprop_maximal_cut_targets']=checks
payload={'metadata':report,'basis':['I'+n[1:] for n in first['basis']],'reductions':a}
(base/'e6_full/fi_dot2_modp.json.gz').write_bytes(gzip.compress((json.dumps(payload,separators=(',',':'))+'\n').encode()))
(base/'transfer_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
