import sys,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from reference_space import RatFunc,RatMat
from analytic_bootstrap.rational import factor_powers,build_local_expansion
from analytic_bootstrap.fp import Kd
from flint import nmod_mat
out=ROOT/'results/dp_diagnosis';c=json.loads((out/'capsule.json').read_text());p=c['prime'];B=RatMat({(i,j):RatFunc.make(v['numerator_tau'],v['denominator_tau'],p) for i,j,v in c['B']},tuple(c['B_shape']),p);qs={}
for ij,r in B.entries.items():
 for q,k in factor_powers(r.den,p):
  if len(q)>2:qs.setdefault(q,{})[ij]=k
reports=[]
for q,powers in qs.items():
 t=time.perf_counter()
 if max(powers.values())!=1:reports.append({'q':q,'higher_pole':max(powers.values())});continue
 exp=build_local_expansion(B,q,Kd(p,q),0,powers);m=exp.coeffs[-1];ch=nmod_mat(m.table(),p).charpoly();eig=[i for i in range(-15,61) if ch(i%p)==0]
 rec={'q':q,'dimension':m.nrows,'small_integer_eigenvalues':eig,'seconds':time.perf_counter()-t};reports.append(rec);(out/'nonlinear_residues.json').write_text(json.dumps(reports,indent=2));print(json.dumps(rec),flush=True)
