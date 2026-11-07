#!/usr/bin/env python3
"""Check reconstructed coefficients against fresh coupled, uncut IBP points."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'j_integral_global'


def field(text,p):
    q=Fraction(str(text));return q.numerator*pow(q.denominator,-1,p)%p


def check(data,point):
    p=data['prime'];tau=(field(point['t_point'],p)-field(data['t0'],p))%p
    assert data['masters']==point['masters_at_point']
    def evaluate(r):
        def poly(v):
            out=0
            for x in reversed(v):out=(out*tau+x)%p
            return out
        return poly(r['numerator_tau'])*pow(poly(r['denominator_tau']),-1,p)%p
    checks=0;errors=[]
    for name,row in point['reduction_values'].items():
        expected={name:1} if name in data['masters'] else {k:evaluate(v) for k,v in data['reductions'].get(name,{}).items()}
        if any((row.get(k,0)-expected.get(k,0))%p for k in set(row)|set(expected)):errors.append(('reduction',name))
        checks+=1
    for name,row in point['DE_values'].items():
        expected={k:evaluate(v) for k,v in data['de'][name]['rhs'].items()}
        if any((row.get(k,0)-expected.get(k,0))%p for k in set(row)|set(expected)):errors.append(('DE',name))
        checks+=1
    return checks,errors


def main():
    p=argparse.ArgumentParser();p.add_argument('family',type=Path);p.add_argument('basis',type=Path);p.add_argument('targets',type=Path);p.add_argument('reduction',type=Path);p.add_argument('output',type=Path);p.add_argument('--workers',type=int,default=2);p.add_argument('--binary',type=Path,default=HERE/'build/global_point');args=p.parse_args();data=json.loads(args.reduction.read_text());args.output.mkdir(parents=True,exist_ok=True)
    points=['1/7','2/7','3/7','4/7','5/7','6/7','2/9','4/9']
    def worker(text):
        name=text.replace('/','over');out=args.output/(name+'.json')
        with (args.output/(name+'.log')).open('w') as log:
            subprocess.run([str(args.binary.resolve()),str(args.family.resolve()),str(out.resolve()),'--point',text,'--targets-file',str(args.targets.resolve()),'--basis-file',str(args.basis.resolve()),'--values','1','--compact','1'],stdout=log,stderr=subprocess.STDOUT,check=True)
        point=json.loads(out.read_text());assert not point['unclosed_derivatives'];checks,errors=check(data,point);assert not errors,(text,errors)
        return {'point':text,'checks':checks,'all_pass':True}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:results=list(pool.map(worker,points))
    # A concrete unequal coefficient must be rejected by the checker.
    bad=json.loads(json.dumps(data));changed=False
    for row in bad['reductions'].values():
        if row:
            coefficient=next(iter(row.values()));coefficient['numerator_tau'][0]=(coefficient['numerator_tau'][0]+1)%data['prime'];changed=True;break
    assert changed
    point=json.loads((args.output/'1over7.json').read_text());assert check(bad,point)[1]
    report={'all_pass':True,'prime':data['prime'],'points':results,'deliberate_unequal_coefficient_rejected':True};(args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
