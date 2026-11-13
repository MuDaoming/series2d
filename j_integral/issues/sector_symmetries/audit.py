"""Exact Feynman-parameter automorphisms of sectors with more J than FI masters."""
import json,gzip,itertools,collections
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[3]

def polynomials(f):
    props=f['family']['propagators'];n=len(props);x=s.symbols('x:'+str(n));gram=s.Matrix([[s.Rational(q) for q in row] for row in f['family']['gram']]);A=s.zeros(2);B=s.zeros(2,gram.rows);C=0
    for xi,p in zip(x,props):
        l=s.Matrix(p['loop']);q=s.Matrix(p['shift']);A+=xi*l*l.T;B+=xi*l*q.T;C+=xi*((q.T*gram*q)[0]-s.Rational(p['mass2']))
    U=s.expand(A.det());F=s.expand((A.adjugate()*B*gram*B.T).trace()-U*C)
    return x,U,F

def audit(case):
    base=ROOT/'j_integral';f=json.load(open(base/f'tests/{case}/input/family_request.json'));br=f['j']['branches'];n=len(br);x,U,F=polynomials(f)
    fi=json.load(open(base/f'archive/j_integral_bootstrap/reference/fi_master_counts/{case}.json'));de=json.loads(gzip.decompress((base/f'tests/{case}/reference/physical_de.json.gz').read_bytes()))
    idx=lambda name:tuple(map(int,name.split('{')[1][:-1].split(',')))
    Fm=collections.defaultdict(list);Jm=collections.defaultdict(list)
    for nu in fi['master_indices']:Fm[tuple(int(q>0) for q in nu[:n])].append(nu[:n])
    for name in de['masters']:
        if name.startswith('J'):nu=idx(name);Jm[tuple(int(q>0) for q in nu)].append(list(nu))
    rows=[]
    for mask in sorted(set(Fm)|set(Jm),key=lambda m:(sum(m),m)):
        if len(Jm[mask])<=len(Fm[mask]):continue
        active=[i for i,q in enumerate(mask) if q];groups=[[i for i in active if br[i]==b] for b in range(3)];zero={x[i]:0 for i in range(n) if not mask[i]};pf=s.Poly(F.subs(zero),*[x[i] for i in active]);terms=pf.terms();original=dict(terms);pos={i:k for k,i in enumerate(active)};autos=[]
        for bp in itertools.permutations(range(3)):
            if any(len(groups[b])!=len(groups[bp[b]]) for b in range(3)):continue
            for orders in itertools.product(*(itertools.permutations(groups[bp[b]]) for b in range(3))):
                mapping={i:j for b in range(3) for i,j in zip(groups[b],orders[b])};per=[pos[mapping[i]] for i in active];trans={}
                for ex,c in terms:
                    out=[0]*len(active)
                    for k,v in enumerate(ex):out[per[k]]=v
                    trans[tuple(out)]=c
                if trans==original:
                    autos.append(dict(branch_map=list(bp),propagator_map={str(i+1):mapping[i]+1 for i in active}))
        rows.append(dict(mask=''.join(map(str,mask)),branch_sizes=list(map(len,groups)),FI=len(Fm[mask]),J=len(Jm[mask]),FI_masters=Fm[mask],J_masters=Jm[mask],automorphisms=autos,branch_exchange=any(q['branch_map']!=[0,1,2] for q in autos),nontrivial_fixed_branch=any(q['branch_map']==[0,1,2] and any(int(i)!=j for i,j in q['propagator_map'].items()) for q in autos)))
        print(case,rows[-1]['mask'],rows[-1]['FI'],rows[-1]['J'],'autos',len(autos),'branch_exchange',rows[-1]['branch_exchange'],flush=True)
    return rows
if __name__=='__main__':
    result={case:audit(case) for case in ['e6','dp']}
    Path(__file__).with_name('audit.json').write_text(json.dumps(result,indent=2)+'\n')
