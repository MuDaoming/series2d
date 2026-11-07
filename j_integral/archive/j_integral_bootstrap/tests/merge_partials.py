"""Merge rational partial answers, rejecting conflicts and checking numerical data.

Inputs and finite field are explicit. This program never chooses or acquires
verification points; pass only the intended learning/check datasets.
"""
import argparse,json
from fractions import Fraction
from pathlib import Path
from reference_space import RatFunc


def main(a):
    cfg=json.loads(Path(a.family).read_text());p=cfg['prime'];t0=Fraction(cfg['t0']);center=t0.numerator*pow(t0.denominator,-1,p)%p
    masters=json.loads(Path(a.basis).read_text());known={sec:{} for sec in a.section};points={};support={sec:{} for sec in a.section}
    def rf(v):return RatFunc.make(v['numerator_tau'],v['denominator_tau'],p)
    for path in a.input:
        data=json.loads(Path(path).read_text())
        if data['prime']!=p or Fraction(data['t0'])!=t0:raise ValueError('input field/coordinate differs: '+path)
        if 'masters' in data and data['masters']!=masters:raise ValueError('input basis differs: '+path)
        for sec in a.section:
            for m,row in data.get(sec,{}).items():
                dest=known[sec].setdefault(m,{})
                for n,v in (row['rhs'] if sec=='de' else row).items():
                    if n not in masters:raise ValueError('unknown output master: '+n)
                    value=rf(v)
                    if n in dest and not (dest[n]-value).is_zero():raise ValueError('rational sources conflict: '+m+' / '+n)
                    dest[n]=value
    for directory in a.points:
        for path in Path(directory).glob('point*.json'):
            if '.tmp.' in path.name:continue
            data=json.loads(path.read_text())
            if data['prime']!=p or data['masters_at_point']!=masters or data['unclosed_derivatives']:raise ValueError('numeric identity/closure differs: '+str(path))
            t=int(data['t_point']);point=points.setdefault(t,{sec:{} for sec in a.section})
            for sec in a.section:
                field='DE_values' if sec=='de' else 'reduction_values'
                for m,row in data[field].items():
                    if m in point[sec] and point[sec][m]!=row:raise ValueError('numeric providers conflict')
                    point[sec][m]=row;support[sec].setdefault(m,set()).update(n for n,v in row.items() if v)
    output=dict(prime=p,t0=cfg['t0'],masters=masters);report=dict(inputs=a.input,point_directories=a.points,point_values=sorted(points),sections={})
    for sec in a.section:
        if sec=='reductions':
            for m in masters:known[sec].setdefault(m,{})[m]=RatFunc.make([1],[1],p)
        checks=0
        for m,row in known[sec].items():
            for n,r in row.items():
                for t,point in points.items():
                    if m not in point[sec]:continue
                    checks+=1
                    if r.evaluate((t-center)%p)!=point[sec][m].get(n,0):raise ValueError('coefficient mismatch: '+m+' / '+n+' t='+str(t))
        missing={m:sorted(cols-set(known[sec].get(m,{}))) for m,cols in support[sec].items()}
        missing={m:cols for m,cols in missing.items() if cols}
        output[sec]={m:({'rhs':{n:dict(numerator_tau=r.num,denominator_tau=r.den) for n,r in row.items() if not r.is_zero()}} if sec=='de' else {n:dict(numerator_tau=r.num,denominator_tau=r.den) for n,r in row.items() if not r.is_zero()}) for m,row in known[sec].items()}
        report['sections'][sec]=dict(nonzero_known=sum(not v.is_zero() for row in known[sec].values() for v in row.values()),pointwise_checks=checks,missing_coefficients=sum(map(len,missing.values())),complete_rows=len(support[sec])-len(missing),observed_rows=len(support[sec]),missing=missing)
    if a.section==['de'] and set(support['de'])==set(masters) and not report['sections']['de']['missing']:
        output.update(status='complete',known_rows=masters)
        for m in masters:output['de'].setdefault(m,{'rhs':{}})
    else:output['status']='partial'
    for path,data in [(Path(a.output),output),(Path(a.output).with_suffix('.checks.json'),report)]:
        path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp.json');tmp.write_text(json.dumps(data,separators=(',',':'))+'\n');tmp.replace(path)
    print(json.dumps({sec:{k:v for k,v in val.items() if k!='missing'} for sec,val in report['sections'].items()}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',required=True);p.add_argument('--basis',required=True);p.add_argument('--input',action='append',required=True);p.add_argument('--points',action='append',required=True);p.add_argument('--section',action='append',choices=['de','reductions'],required=True);p.add_argument('--output',required=True)
    main(p.parse_args())
