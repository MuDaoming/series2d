"""Compare complete rational connections by named masters, including zeros."""
import argparse,json
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'.deps'))
from flint import nmod_poly

def compare(left,right):
    if left.get('status')!='complete' or right.get('status')!='complete':raise ValueError('both connections must be complete')
    if left['prime']!=right['prime'] or Fraction(left['t0'])!=Fraction(right['t0']):raise ValueError('field/coordinate differs')
    names=left['masters'];p=left['prime']
    if len(set(names))!=len(names) or len(set(right['masters']))!=len(right['masters']) or set(names)!=set(right['masters']):raise ValueError('master identities differ')
    if set(left['de'])!=set(names) or set(right['de'])!=set(names):raise ValueError('row coverage differs')
    nonzero_checks=0;zero={'numerator_tau':[0],'denominator_tau':[1]}
    for m in names:
        a=left['de'][m]['rhs'];b=right['de'][m]['rhs']
        for n in set(a)|set(b):
            if n not in names:raise ValueError('unknown master column')
            x=a.get(n,zero);y=b.get(n,zero);nonzero_checks+=1
            if x==y:continue
            if nmod_poly(x['numerator_tau'],p)*nmod_poly(y['denominator_tau'],p)!=nmod_poly(y['numerator_tau'],p)*nmod_poly(x['denominator_tau'],p):raise ValueError('rational coefficient differs: '+m+' / '+n)
    return dict(all_pass=True,prime=p,master_count=len(names),matrix_entries_including_zeros=len(names)**2,explicit_union_entries=nonzero_checks)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('left');p.add_argument('right');p.add_argument('report');a=p.parse_args()
    report=compare(json.loads(Path(a.left).read_text()),json.loads(Path(a.right).read_text()));Path(a.report).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
