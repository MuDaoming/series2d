"""Prepare explicit AB inputs from a saved physical J connection; never solve AB."""
import json,gzip,hashlib
from pathlib import Path
from fractions import Fraction
from analytic_bootstrap import RatFunc,RatMat
from analytic_bootstrap.upstream.model import family_definition,matrix_record
ROOT=Path(__file__).resolve().parent;ARCHIVE=ROOT.parent/'archive'
def load(path):
    if path.suffix=='.gz':
        with gzip.open(path,'rt') as file:return json.load(file)
    return json.loads(path.read_text())
def save(path,data):
    if path.suffix=='.gz':path.write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    else:path.write_text(json.dumps(data,indent=2)+'\n')
def components(graph):
    ids={};low={};stack=[];active=set();answer=[]
    def visit(i):
        ids[i]=low[i]=len(ids);stack.append(i);active.add(i)
        for j in sorted(graph[i]):
            if j not in ids:visit(j);low[i]=min(low[i],low[j])
            elif j in active:low[i]=min(low[i],ids[j])
        if low[i]==ids[i]:
            group=[]
            while True:
                j=stack.pop();active.remove(j);group.append(j)
                if j==i:break
            answer.append(sorted(group))
    for i in graph:
        if i not in ids:visit(i)
    return answer
for case,source,family_source in [
    ('dp','j_integral_bootstrap/reference/dp_global/connection_mod1000003.json.gz','j_integral_global/examples/dp_physical.json'),
    ('e6','j_integral_sixpoint/reference/connection_mod1000003.json','j_integral_sixpoint/inputs/family.json')]:
    dest=ROOT/case/'input';reference=load(ARCHIVE/source);physical=load(ARCHIVE/family_source);names=reference['masters'];index={m:i for i,m in enumerate(names)};p=reference['prime'];q=Fraction(reference['t0']);center=q.numerator*pow(q.denominator,-1,p)%p;gauge=3
    graph={i:{index[m] for m,v in reference['de'][name]['rhs'].items() if any(v['numerator_tau'])} for i,name in enumerate(names)}
    groups=components(graph);owner={i:k for k,group in enumerate(groups) for i in group}
    definition=family_definition(dict(family=dict(name=case+'_J',normalization='M=[t(t-1)]^3 J; tau=t-2/5',basis_identity=hashlib.sha256(json.dumps(names).encode()).hexdigest()),variable='tau',sectors=[dict(id='block'+str(k),masters=[names[i] for i in group],depends_on=sorted({'block'+str(owner[j]) for i in group for j in graph[i] if owner[j]!=k})) for k,group in enumerate(groups)],reductions={}))
    canonical=[m for sector in definition['sectors'] for m in sector['masters']];position={m:i for i,m in enumerate(canonical)}
    entries={(position[row],position[col]):RatFunc.make(v['numerator_tau'],v['denominator_tau'],p) for row,record in reference['de'].items() for col,v in record['rhs'].items()}
    shift=RatFunc.make([gauge],[center,1],p)+RatFunc.make([gauge],[(center-1)%p,1],p);zero=RatFunc.make([0],[1],p)
    for i in range(len(names)):entries[i,i]=entries.get((i,i),zero)+shift
    matrix=RatMat(entries,(len(names),len(names)),p)
    diagonals={}
    for sector in definition['sectors']:
        indices=[position[m] for m in sector['masters']];coords={j:i for i,j in enumerate(indices)}
        block=RatMat({(coords[i],coords[j]):v for (i,j),v in entries.items() if i in coords and j in coords},(len(indices),len(indices)),p)
        diagonals[sector['id']]=dict(masters=sector['masters'],matrix=matrix_record(block))
    E=len(physical['gram']);family=dict(name=case,loops=['k1','k2'],external_momenta=['p'+str(i+1) for i in range(E)],gram=physical['gram'],propagators=[{k:v for k,v in prop.items() if k!='branch'} for prop in physical['propagators']])
    parents=sorted({tuple(map(int,name[name.index('{')+1:-1].split(','))) for name in names})
    ibp=dict(schema_version=1,family=family,j=dict(coordinates='tz',branches=[prop['branch'] for prop in physical['propagators']],boundary=dict(z0='fbi_analytic',z1='regulated_zero')),point=dict(prime=p,dimension=physical['dimension'],t='17',parameters={}),request=dict(targets=parents,basis=names,derivatives='basis'),options=dict(completion=16,germ_order=16))
    save(dest/'family_request.json',ibp);save(dest/'basis.json',names);save(dest/'definition.json',definition);save(dest/'setting.json',dict(prime=p,assignments={}))
    save(dest/'anchor.json.gz',dict(prime=p,assignments={},masters=canonical,de=matrix_record(matrix),reductions={}))
    save(dest/'diagonals.json',dict(prime=p,assignments={},blocks=diagonals))
    save(dest/'normalization.json',dict(variable='tau',t0=reference['t0'],gauge_power=gauge,physical_basis=names))
    save(dest/'provenance.json',dict(reference_source=source,reference_sha256=hashlib.sha256((ARCHIVE/source).read_bytes()).hexdigest(),family_source=family_source,physical_master_count=len(names),block_count=len(groups),nonzero_physical_DE=sum(len(row['rhs']) for row in reference['de'].values()),FBI_branches=[sum(prop['branch']==b for prop in physical['propagators']) for b in range(3)],point=ibp['point'],normalization='uniform gauge only; no known off-diagonal coefficient enters numeric backend'))
    print(case,len(names),'masters',len(groups),'blocks',len(parents),'native J parents',flush=True)
