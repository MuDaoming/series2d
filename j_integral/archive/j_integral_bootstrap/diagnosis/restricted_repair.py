"""Exact principal-part membership of true C in old and minimally repaired AB spaces."""
import runpy,json,time
from pathlib import Path
s=runpy.run_path(str(Path(__file__).with_name('restricted_envelope.py')));globals().update({k:v for k,v in s.items() if not k.startswith('__')})
from reference_space import BootstrapSpec,RatMat,build_pl_ansatz
from analytic_bootstrap.ansatz import _native_conditions_eliminate
from analytic_bootstrap.sparse import iter_sparse_rows
from analytic_bootstrap.rational import factor_valuation
from flint import nmod_poly
pos={m:i for i,m in enumerate(sub)};toppos={m:i for i,m in enumerate(top)}
A={(toppos[m],toppos[j]):v for (m,j),v in D.items() if m in toppos and j in toppos};B={(pos[m],pos[j]):v for (m,j),v in D.items() if m in pos and j in pos}
h=ep**3;shift=RatFunc.make([int(v) for v in h.derivative()],[int(v) for v in h],p);zero=RatFunc.make([0],[1],p)
for i in range(len(top)):A[i,i]=A.get((i,i),zero)+shift
for i in range(len(sub)):B[i,i]=B.get((i,i),zero)+shift
C={(toppos[m],pos[j]):v for (m,j),v in D.items() if m in toppos and j in pos}
ds=[];caps=[]
for m in top:
 dr=[];cr=[]
 for j in sub:
  if (m,j) not in D:dr.append([1]);cr.append(-1)
  elif (m,j)in learned:dr.append(learned[m,j].den);cr.append(len(learned[m,j].num)-1)
  else:
   den=L.monic((qs*cols[j]*ep)**2);dr.append([int(v) for v in den]);cr.append(int(den.degree())+3)
 ds.append(dr);caps.append(cr)

def check(a):
 results=[]
 for b in a.blocks:
  q=tuple(b.factor);qp=nmod_poly(list(q),p);principal={};missing=[];dest={tuple(coord):i for i,coord in enumerate(b.coordinates)}
  for (i,j),r in C.items():
   v=factor_valuation(r.den,q,p)
   if not v:continue
   mod=qp**v;den=nmod_poly(r.den,p)//mod;g,inv,_=den.xgcd(mod)
   if len(g)!=1 or int(g[0])!=1:raise ArithmeticError('non-coprime denominator')
   numerator=(nmod_poly(r.num,p)*inv)%mod
   for k in range(v):
    numerator,rem=divmod(numerator,qp)
    for cp,val in enumerate(rem):
     val=int(val)
     if not val:continue
     coord=(i,j,v-k,cp)
     if coord not in dest:missing.append({'row':top[i],'col':sub[j],'pole':v-k,'power':cp,'value':val})
     else:principal[dest[coord]]=val
  rows={i:dict(row) for i,row in enumerate(iter_sparse_rows(b.basis))};pivots=_native_conditions_eliminate(p,b.basis.ncols,rows)
  if pivots is None:raise RuntimeError('native kernel unavailable')
  residual=dict(principal)
  for c,row in sorted(pivots.items()):
   z=residual.get(c,0)
   if z:
    for j,v in row.items():
     w=(residual.get(j,0)-z*v)%p
     if w:residual[j]=w
     else:residual.pop(j,None)
  rec={'factor':q,'basis_dimension':b.basis.nrows,'coordinates':b.basis.ncols,'missing_coordinates':missing,'residual_nonzero':len(residual)}
  results.append(rec)
 return results
reports=[]
for repair in [False,True,'all-unresolved']:
 den=json.loads(json.dumps(ds));cap=json.loads(json.dumps(caps));changes=[]
 if repair:
  q=[700002,1]
  entries=report['settings'][-1]['failures'] if repair is True else [{'row':m,'col':j} for m in top for j in sub if (m,j) in D and (m,j) not in learned]
  for entry in entries:
   i,j=toppos[entry['row']],pos[entry['col']];den[i][j]=[int(v) for v in L.poly(den[i][j])*L.poly(q)];cap[i][j]+=1;changes.append([entry['row'],entry['col']])
 spec=BootstrapSpec.from_inputs(family={'source':'restricted-historical-diagnosis'},variable='tau',assignments={},masters=sub+top,target=None,de={'A':RatMat(A,(6,6),p),'B':RatMat(B,(47,47),p)},denominators=den,numerator_caps=cap,prime=p,mode='de')
 t=time.perf_counter();a=build_pl_ansatz(spec,backend='PL-Hasse',laurent_mode='structural-refined-all');checks=check(a)
 rec={'repair':repair,'changed_entries':changes,'ambient_naive':sum(v+1 for row in cap for v in row),'ansatz_dimension':a.dimension,'seconds':time.perf_counter()-t,'factor_checks':checks,'all_membership_pass':all(not c['missing_coordinates'] and not c['residual_nonzero'] for c in checks)};
 if repair=='all-unresolved':
  from native_fit import joint_fit
  names=dict(zip(mapping['original_masters'],mapping['compressed_masters']))
  points=[json.loads((r/f'dp_restricted/samples/points/point{t}.json').read_text()) for t in range(17,82)]
  samples=[{'prime':p,'x':(int(point['t_point'])-center)%p,'values':[[point['DE_values'][names[m]].get(names[j],0) for j in sub] for m in top]} for point in points]
  fixed={(toppos[m],pos[j]):v for (m,j),v in learned.items() if m in toppos and j in pos and (m,j) in D}
  fit=joint_fit(a,samples[:-3],fixed,[v['x'] for v in samples]);rec['fit']={k:v for k,v in fit.items() if k!='parameters'}
  if fit['parameters'] is not None:
   answer=a.assemble(fit['parameters']);rec['holdouts_pass']=all(a.evaluate(v['x'],fit['parameters'])==v['values'] for v in samples[-3:]);rec['exact_functions_match']=all((answer[i][j]-C.get((i,j),zero)).is_zero() for i in range(6) for j in range(47))
 reports.append(rec);(out/'repair_membership.json').write_text(json.dumps(reports,indent=2));print(json.dumps({k:v for k,v in rec.items() if k!='factor_checks'}),flush=True)
