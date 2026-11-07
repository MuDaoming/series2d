import sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from reference_space import RatFunc
from analytic_bootstrap.rational import factor_powers
from analytic_bootstrap.poly import poly_divmod,poly_eval
c=json.loads((ROOT/'results/dp_diagnosis/capsule.json').read_text());p=c['prime'];report=[]
for i,j,v in c['A']:
 if i!=j:continue
 rf=RatFunc.make(v['numerator_tau'],v['denominator_tau'],p)
 for q,k in factor_powers(rf.den,p):
  if len(q)!=2:continue
  root=-q[0]*pow(q[1],-1,p)%p
  den,_=poly_divmod(rf.den,q,p)
  if k!=1:continue
  residue=poly_eval(rf.num,root,p)*pow(poly_eval(den,root,p),-1,p)%p
  if residue>=p-30 and q not in [(200000,1),(200001,1)]:report.append({'i':i,'master':c['masters'][i],'q':q,'residue':residue-p,'boundary':not c['masters'][i].startswith('J{')})
(ROOT/'results/dp_diagnosis/top_negative_diagonal_residues.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
from flint import nmod_mat
B={(i,j):RatFunc.make(v['numerator_tau'],v['denominator_tau'],p) for i,j,v in c['A']};qs={}
for ij,r in B.items():
 for q,k in factor_powers(r.den,p):
  if len(q)==2 and q not in [(200000,1),(200001,1)]:qs.setdefault(q,[]).append((ij,k))
negative=[];skipped=[]
for q,entries in qs.items():
 if max(k for ij,k in entries)>1:skipped.append(q);continue
 mat=[[0]*6 for _ in range(6)];root=-q[0]*pow(q[1],-1,p)%p
 for (i,j),k in entries:
  r=B[i,j];den,_=poly_divmod(r.den,q,p);mat[i][j]=poly_eval(r.num,root,p)*pow(poly_eval(den,root,p),-1,p)%p
 m=nmod_mat(mat,p);ch=m.charpoly()
 for a in range(-60,16):
  if ch((-a)%p)==0:negative.append({'q':q,'integer_eigenvalue':-a})
out={'integer_residues':negative,'higher_pole_skipped':skipped,'factors_tested':len(qs)-len(skipped)};(ROOT/'results/dp_diagnosis/top_integer_residue_eigenvalues.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
