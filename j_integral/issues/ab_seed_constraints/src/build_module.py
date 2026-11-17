"""Build the finite-point numerator module with the unchanged AB public API."""
import argparse,json,time
from pathlib import Path
from analytic_bootstrap import BootstrapSpec,RatMat,RatFunc,build_reduction_ansatz
from analytic_bootstrap.poly import poly_mul,poly_divmod,poly_gcd
from analytic_bootstrap.ansatz import SparseConditions

parser=argparse.ArgumentParser();parser.add_argument('snapshot',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--physical-z0',action='store_true');parser.add_argument('--normalize-z0',action='store_true');parser.add_argument('--depth-margin',type=int,default=4);args=parser.parse_args()
started=time.perf_counter();data=json.loads(args.snapshot.read_text());p=data['prime'];rank=len(data['masters'])
def rf(x):return RatFunc.make(x['n'],x['d'],p)
def div(a,b):
 q,r=poly_divmod(a,b,p)
 if any(r):raise ValueError('inexact polynomial quotient')
 return q
weights=[RatFunc.make([0]*int(order)+[1],[1],p) if args.normalize_z0 else RatFunc.make([1],[1],p) for order in data['orders']]
rows={name:[rf(x)*weights[j] for j,x in enumerate(row)] for name,row in data['rows'].items()}
entries={(i,j):rf(x)*weights[j]/weights[i] for i,row in enumerate(data['connection']) for j,x in enumerate(row) if any(x['n'])}
for i in range(rank):
 if args.normalize_z0:
  entries[i,i]=entries.get((i,i),RatFunc.make([0],[1],p))-weights[i].derivative()/weights[i]
connection=RatMat(entries,(rank,rank),p)
dens=[];caps=[]
envelope_rows=list(rows.values())+[[rf(x)*weights[j] for j,x in enumerate(row)] for row in data.get('reference_seeds',[])]
for j in range(rank):
 den=[1]
 for row in envelope_rows:den=poly_mul(den,div(row[j].den,poly_gcd(den,row[j].den,p)),p)
 cap=max(len(poly_mul(row[j].num,div(den,row[j].den),p))-1 for row in envelope_rows)
 # Constant primitives occur in the current resonant Hermite completion.
 cap=max(cap,len(den)-1)
 dens.append(den);caps.append(cap)
offsets=[];ambient=0
for cap in caps:offsets.append(ambient);ambient+=cap+1
spec=BootstrapSpec.from_inputs(family={'name':'FBI-integrand-'+str(data['source'])},variable='z',assignments={},target=['primitive'],masters=data['masters'],de=connection,denominators=[dens],numerator_caps=[caps],prime=p,mode='reduction')
args.output.mkdir(parents=True,exist_ok=True)
print('ambient numerator parameters',ambient,'factors',len(spec.envelope.factors),flush=True)
def progress(v):print('AB factor',v,flush=True)
ansatz=build_reduction_ansatz(spec,depth_margin=args.depth_margin,progress=progress)
if args.physical_z0:
 from analytic_bootstrap import Fp,FpMatrix,Kd,annihilator_module
 from analytic_bootstrap.jets import JetResult
 from analytic_bootstrap.ansatz import FactorBlock,RationalAnsatz,impose_caps
 from analytic_bootstrap.linalg import nullspace
 blocks=[];offset=0
 for block in ansatz.blocks:
  if block.factor==(0,1):
   k=max(coord[2] for coord in block.coordinates)
   constants=sorted({key for power in range(k) for j in range(rank) for key in data['z0_physical_jets'][power+(data['orders'][j] if args.normalize_z0 else 0)][j]})
   matrix=FpMatrix.zeros(Fp(p),len(constants),rank*k)
   for a,key in enumerate(constants):
    for power in range(k):
     for j in range(rank):
      source=data['z0_physical_jets'][power+(data['orders'][j] if args.normalize_z0 else 0)][j]
      matrix.m[a,power*rank+j]=source.get(key,0)
   jet=JetResult(kd=Kd(p,(0,1)),n=rank,k=k,depth=16,jets=matrix,method='physical-FBI-analytic-endpoint')
   full=annihilator_module(jet,backend='dense')
   columns=[(k-pole)*rank+j for a,j,pole,cp in block.coordinates]
   forbidden=sorted(set(range(rank*k))-set(columns))
   if forbidden:full=nullspace(full.submatrix(cols=forbidden).transpose())@full
   newbasis=full.submatrix(cols=columns)
   block=FactorBlock(block.factor,block.coordinates,newbasis,offset,{'source':'unchanged AB annihilator_module on existing physical FBI endpoint jets'})
  block.offset=offset;offset+=block.basis.nrows;blocks.append(block)
 ansatz=RationalAnsatz(spec,blocks,spec.envelope.regular_records,impose_caps(spec,blocks),'finite-window-unconfirmed',{'physical_z0':True,'other_finite_factors':'unchanged AB build_reduction_ansatz'})
module_seconds=time.perf_counter()-started
print('AB parameters',ansatz.dimension,'seconds',module_seconds,flush=True)
def encode(row):return [{'n':x.num,'d':x.den} for x in row]
def vector(row):
 v={}
 for j,x in enumerate(row):
  n=poly_mul(x.num,div(dens[j],x.den),p)
  if len(n)>caps[j]+1:raise ValueError('numerator cap exceeded')
  for k,c in enumerate(n):
   if c:v[offsets[j]+k]=c
 return v
span=SparseConditions(p,ambient);allowed=[]
for i in range(ansatz.dimension):
 parameters=[0]*ansatz.dimension;parameters[i]=1
 row=ansatz.assemble(parameters)[0];span.add(vector(row));allowed.append(encode([x/weights[j] for j,x in enumerate(row)]))
if span.rank!=ansatz.dimension:raise ValueError('assembled module basis is dependent')
def remainder(v):
 v=dict(v)
 for pivot,row in sorted(span.pivots.items()):
  a=v.get(pivot,0)
  if a:
   for c,x in row.items():
    value=(v.get(c,0)-a*x)%p
    if value:v[c]=value
    else:v.pop(c,None)
 return v
unrepresented=[name for name,row in rows.items() if remainder(vector(row))]
trivial=[]
for j in range(rank):
 for k in range(caps[j]+1):
  row=[{'n':[],'d':[1]} for _ in range(rank)];row[j]={'n':[0]*k+[1],'d':dens[j]};trivial.append(encode([rf(x)/weights[j] for j,x in enumerate(row)]))
(args.output/'ab_seeds.json').write_text(json.dumps({'rows':allowed},separators=(',',':'))+'\n')
(args.output/'trivial_seeds.json').write_text(json.dumps({'rows':trivial},separators=(',',':'))+'\n')
summary={'physical_z0':args.physical_z0,'normalized_z0':args.normalize_z0,'prime':p,'point':data['point'],'rank':rank,'ambient_parameters':ambient,'AB_parameters':ansatz.dimension,'denominator_degrees':[len(d)-1 for d in dens],'numerator_caps':caps,'finite_factors':[list(f.factor) for f in spec.envelope.factors],'unrepresented_physical_rows':unrepresented,'tested_physical_rows':len(rows),'depth_margin':args.depth_margin,'module_seconds':module_seconds,'total_seconds':time.perf_counter()-started,'AB_status':ansatz.status,'AB_diagnostics':ansatz.diagnostics,'condition':'finite analytic branch has no negative integer powers; no infinity-region constraints'}
(args.output/'module_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('unrepresented physical rows',len(unrepresented),'total seconds',summary['total_seconds'],flush=True)
