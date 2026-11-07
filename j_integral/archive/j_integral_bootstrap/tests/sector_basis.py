"""Choose a basis respecting sector filtration from existing numerical rows.

Every available lower-sector dot target precedes every higher-sector target.
This differs from sorting all targets solely by the sum of propagator powers.
"""
import argparse
from collections import Counter
import json
from pathlib import Path


def indices(name):
    return json.loads('['+name[name.index('{')+1:name.index('}')]+']')


def key(name):
    nu=indices(name)
    return (name[0]!='B',sum(v>0 for v in nu),sum(nu),name)


def insert(pivots,values,p):
    row={k:v%p for k,v in values.items() if v%p}
    while row:
        k=min(row)
        if k not in pivots:
            inv=pow(row[k],-1,p)
            pivots[k]={j:v*inv%p for j,v in row.items()}
            return True
        c=row[k]
        for j,v in pivots[k].items():
            value=(row.get(j,0)-c*v)%p
            if value:row[j]=value
            else:row.pop(j,None)
    return False


def main(args):
    first=json.loads(Path(args.points[0]).read_text());p=first['prime']
    old=first['masters_at_point'];positions={m:i for i,m in enumerate(old)}
    pivots={};chosen=[]
    for name in sorted(first['reduction_values'],key=key):
        values={positions[m]:v for m,v in first['reduction_values'][name].items()}
        if insert(pivots,values,p):chosen.append(name)
    if len(chosen)!=len(old):raise AssertionError('new basis loses rank')
    checks=[]
    for path in args.points[1:]:
        data=json.loads(Path(path).read_text());p=data['prime'];positions={m:i for i,m in enumerate(data['masters_at_point'])}
        pivots={}
        for name in chosen:
            values={positions[m]:v for m,v in data['reduction_values'][name].items()}
            if not insert(pivots,values,p):raise AssertionError('basis loses rank at independent point')
        checks.append(dict(path=path,prime=p,t_point=data['t_point'],rank=len(pivots)))
    Path(args.output).write_text(json.dumps(chosen,indent=2)+'\n')
    report=dict(source=args.points[0],master_count=len(chosen),
                changed=len(set(chosen)-set(old)),
                J_sector_count=len({tuple(v>0 for v in indices(m)) for m in chosen if m[0]=='J'}),
                by_active_lines=dict(Counter(sum(v>0 for v in indices(m)) for m in chosen if m[0]=='J')),
                independent_rank_checks=checks)
    Path(args.output+'.report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('output');ap.add_argument('points',nargs='+')
    main(ap.parse_args())
