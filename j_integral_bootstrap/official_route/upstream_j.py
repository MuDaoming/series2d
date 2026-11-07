"""J-IBP capability adapters for the unmodified upstream production pipeline.

The anchor capability receives a complete rational connection once. Pointwise
production receives only numerical IBP and diagonal capabilities, not the anchor.
"""
import sys,json,subprocess,fcntl,hashlib
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tests'))
from reference_space import RatFunc,RatMat
from sector_de import components
from analytic_bootstrap.upstream.model import matrix_record,submatrix,family_definition
from analytic_bootstrap.storage import digest


def definition(answer):
    names=answer['masters'];index={m:i for i,m in enumerate(names)}
    graph={i:{index[n] for n,v in answer['de'][m]['rhs'].items() if any(v['numerator_tau'])} for i,m in enumerate(names)}
    groups=components(graph);owner={i:k for k,g in enumerate(groups) for i in g}
    sectors=[dict(id='block'+str(k),masters=[names[i] for i in g],depends_on=sorted({'block'+str(owner[j]) for i in g for j in graph[i] if owner[j]!=k})) for k,g in enumerate(groups)]
    return family_definition(dict(family={'name':'J-IBP','basis_identity':digest(names)},variable='tau',sectors=sectors,reductions={}))


def connection(answer,gauge):
    p=answer['prime'];names=answer['masters'];idx={m:i for i,m in enumerate(names)}
    entries={(idx[m],idx[n]):RatFunc.make(v['numerator_tau'],v['denominator_tau'],p) for m,row in answer['de'].items() for n,v in row['rhs'].items()}
    shift=normalization_shift(p,answer['t0'],gauge)
    for i in range(len(names)):entries[i,i]=entries.get((i,i),RatFunc.make([0],[1],p))+shift
    return RatMat(entries,(len(names),len(names)),p)


def normalization_shift(p,t0,gauge):
    q=Fraction(t0);c=q.numerator*pow(q.denominator,-1,p)%p
    return RatFunc.make([gauge],[c,1],p)+RatFunc.make([gauge],[(c-1)%p,1],p)


def physical_connection(matrix,names,t0,gauge):
    shift=normalization_shift(matrix.p,t0,gauge);zero=RatFunc.make([0],[1],matrix.p)
    de={m:{'rhs':{}} for m in names}
    for i,j in set(matrix.entries)|{(i,i) for i in range(len(names))}:
        rf=matrix.entries.get((i,j),zero)-(shift if i==j else zero)
        if not rf.is_zero():de[names[i]]['rhs'][names[j]]=dict(numerator_tau=rf.num,denominator_tau=rf.den)
    return de


class AnchorBackend:
    def __init__(self,answer,gauge=0):
        if answer.get('status')!='complete':raise ValueError('a complete anchor is required by learn_anchor')
        self.names=answer['masters'];self.matrix=connection(answer,gauge);self.prime=answer['prime']
    def anchor(self,req):
        if req['prime']!=self.prime or req['reductions']:raise ValueError('anchor setting or unsupported reductions')
        inds=[self.names.index(m) for m in req['masters']]
        return dict(request_sha256=digest(req),de=matrix_record(submatrix(self.matrix,inds,inds)),reductions={})


class DiagonalBackend:
    def __init__(self,anchor,definition):
        # Copy only diagonal matrices: no off-diagonal answer is retained.
        self.blocks={}
        for sector in definition['sectors']:
            inds=[anchor.names.index(m) for m in sector['masters']]
            self.blocks[sector['id']]=submatrix(anchor.matrix,inds,inds)
    def diagonal(self,req):
        return dict(request_sha256=digest(req),matrix=matrix_record(self.blocks[req['sector']]))


class NumericBackend:
    numeric_batch_size=1
    def __init__(self,family,basis,targets,binary,work,*,point_directories=(),gauge=0):
        self.cfg=json.loads(Path(family).read_text());self.names=json.loads(Path(basis).read_text())
        self.family=str(Path(family).resolve());self.basis=str(Path(basis).resolve());self.targets=str(Path(targets).resolve());self.binary=str(Path(binary).resolve())
        self.work=Path(work);self.work.mkdir(parents=True,exist_ok=True);self.gauge=gauge
        p=self.cfg['prime'];q=Fraction(self.cfg['t0']);self.center=q.numerator*pow(q.denominator,-1,p)%p;self.prime=p
        self.points={};self.calls=0
        for directory in point_directories:
            for path in Path(directory).glob('point*.json'):
                if '.tmp.' not in path.name:self.points[int(json.loads(path.read_text())['t_point'])]=path
        identity=dict(family_sha256=hashlib.sha256(Path(family).read_bytes()).hexdigest(),basis=self.names,targets_sha256=hashlib.sha256(Path(targets).read_bytes()).hexdigest(),binary_sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest(),gauge=gauge)
        ip=self.work/'identity.json'
        if ip.exists() and json.loads(ip.read_text())!=identity:raise ValueError('numeric adapter work identity changed')
        ip.write_text(json.dumps(identity,indent=2)+'\n')
    def point(self,x):
        t=(x+self.center)%self.prime;path=self.points.get(t,self.work/f'point{t}.json')
        with (self.work/f'point{t}.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if not path.exists():
                tmp=path.with_suffix('.tmp.json')
                with (self.work/f'point{t}.log').open('w') as log:
                    subprocess.run([self.binary,self.family,str(tmp),'--point',str(t),'--basis-file',self.basis,'--targets-file',self.targets,'--values','1','--compact','1','--integration','atomic'],stdout=log,stderr=subprocess.STDOUT,check=True)
                tmp.replace(path);self.calls+=1
        d=json.loads(path.read_text())
        if d['prime']!=self.prime or int(d['t_point'])!=t or d['masters_at_point']!=self.names or d['unclosed_derivatives'] or set(d['DE_values'])!=set(self.names):raise ValueError('numeric field/basis/closure/coverage mismatch')
        return d
    def prepare(self,req):
        if req['mode']!='de' or req['prime']!=self.prime or set(req['masters'])!=set(self.names):raise ValueError('unsupported numeric request')
        backend=self
        class Session:
            def evaluate(self,x):
                data=backend.point(x);t=(x+backend.center)%backend.prime
                if backend.gauge and (not t or t==1):
                    from analytic_bootstrap.upstream import BadSamplePoint
                    raise BadSamplePoint('normalization endpoint')
                shift=backend.gauge*(pow(t,-1,backend.prime)+pow((t-1)%backend.prime,-1,backend.prime))%backend.prime if backend.gauge else 0
                values=[]
                for target in req['targets']:
                    if target.get('kind')!='derivative' or target['variable']!='tau':raise ValueError('unsupported target kind')
                    m=target['master'];row=data['DE_values'][m]
                    values.append([(row.get(n,0)+(shift if n==m else 0))%backend.prime for n in req['masters']])
                return dict(request_sha256=digest(req),x=x,values=values)
            def close(self):pass
        return Session()
