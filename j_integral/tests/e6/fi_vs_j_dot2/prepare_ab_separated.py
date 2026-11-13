import json,gzip,sys
from pathlib import Path
from flint import nmod_poly
from analytic_bootstrap import RatFunc
base=Path(__file__).resolve().parent;root=base.parents[3];out=base/'j_ab_separated';out.mkdir(exist_ok=True);p=1000003;center=2*pow(5,-1,p)%p;h=26*pow(14,-1,p)%p
point=json.load(open(base/'j_separated_anchor_point.json'));basis=point['basis'];ts=17;targets=json.load(open(base/'input/common_targets.json'));names=['J{'+','.join(map(str,nu))+'}' for nu in targets]
def inds(name):return list(map(int,name.split('{')[1][:-1].split(',')))
def counts(nu):return sum(nu[:2]),sum(nu[5:])
t=nmod_poly([center,1],p);omt=nmod_poly([(1-center)%p,p-1],p);de=[];entries=[];idx={m:i for i,m in enumerate(basis)}
for i,m in enumerate(basis):
 A,C=counts(inds(m));r=RatFunc.make([(h-C-1+3)%p],[center,1],p)+RatFunc.make([(-h+A+1-3)%p],[(1-center)%p,p-1],p);de.append([i,i,[r.num,r.den]])
for i,name in enumerate(names):
 A,C=counts(inds(name))
 for master,value in point['reductions'][name].items():
  Am,Cm=counts(inds(master));a=Cm-C;b=Am-A;num=nmod_poly([value*pow(ts,-a,p)*pow((1-ts)%p,-b,p)%p],p);den=nmod_poly([1],p)
  if a>=0:num*=t**a
  else:den*=t**(-a)
  if b>=0:num*=omt**b
  else:den*=omt**(-b)
  r=RatFunc.make(list(map(int,num)),list(map(int,den)),p);entries.append([i,idx[master],[r.num,r.den]])
anchor={'prime':p,'assignments':{},'masters':basis,'de':{'shape':[61,61],'entries':de},'reductions':{'dot2':{'shape':[1111,61],'entries':entries}}};(out/'anchor.json.gz').write_bytes(gzip.compress(json.dumps(anchor,separators=(',',':')).encode()))
definition={'family':{'name':'e6_J_separated','normalization':'M=[t(t-1)]^3 J; tau=t-2/5; integrated middle branch'},'variable':'tau','sectors':[{'id':'block'+str(i),'masters':[m],'depends_on':[]} for i,m in enumerate(basis)],'reductions':{'dot2':names}};(out/'definition.json').write_text(json.dumps(definition,indent=2)+'\n');diag={'prime':p,'assignments':{},'blocks':{'block'+str(i):{'masters':[m],'matrix':{'shape':[1,1],'entries':[[0,0,de[i][2]]]}} for i,m in enumerate(basis)}};(out/'diagonals.json').write_text(json.dumps(diag,separators=(',',':'))+'\n')
source=json.load(open(base/'input/j_separated_anchor_point.json'));source['request']['derivatives']='basis';(out/'request.json').write_text(json.dumps(source,indent=2)+'\n');norm={'variable':'tau','t0':'2/5','gauge_power':3,'physical_basis':basis};(out/'normalization.json').write_text(json.dumps(norm,indent=2)+'\n')
worker=str(base/'backend_worker.py');config={'anchor_command':[sys.executable,worker,'anchor','--input',str(out/'anchor.json.gz')],'cut_command':[sys.executable,worker,'cut','--input',str(out/'diagonals.json')],'numeric_command':[sys.executable,worker,'numeric','--input',str(out/'request.json'),'--normalization',str(out/'normalization.json'),'--binary',str(root/'j_integral/src/IBP/build/j-ibp'),'--cache',str(out/'numeric')],'learn_profiles':True,'sampling':{'max_samples':64,'holdouts':2,'shared_de_samples':True,'global_holdouts':2,'fit_backend':'stream'}};(out/'backend.json').write_text(json.dumps(config,indent=2)+'\n');print('Separated anchor: constructed from one numerical point and proved t-factorization')
