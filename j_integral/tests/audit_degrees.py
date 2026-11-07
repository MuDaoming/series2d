"""Count ordinary rational-interpolation unknowns after exact Fp polynomial gcd."""
import argparse,gzip,json
from pathlib import Path
from flint import nmod_poly
p=argparse.ArgumentParser();p.add_argument('connection');p.add_argument('output');a=p.parse_args();source=Path(a.connection)
raw=gzip.decompress(source.read_bytes()) if source.suffix=='.gz' else source.read_bytes();data=json.loads(raw);prime=data['prime'];entries=[];cancelled=0
for row,record in data['de'].items():
    for col,rf in record['rhs'].items():
        n=nmod_poly(rf['numerator_tau'],prime);d=nmod_poly(rf['denominator_tau'],prime)
        if n.is_zero():continue
        g=n.gcd(d)
        if g.degree()>0:cancelled+=1;n=n//g;d=d//g
        m,k=n.degree(),d.degree();entries.append(dict(row=row,column=col,numerator_degree=m,denominator_degree=k,total_degree=m+k,interpolation_points=m+k+1))
maximum=max(v['total_degree'] for v in entries)
result=dict(prime=prime,t0=data['t0'],masters=len(data['masters']),nonzero_entries=len(entries),gcd_checked=True,entries_with_cancelled_common_factors=cancelled,maximum_numerator_degree=max(v['numerator_degree'] for v in entries),maximum_denominator_degree=max(v['denominator_degree'] for v in entries),maximum_total_degree=maximum,ordinary_rational_reconstruction_points=maximum+1,with_two_independent_check_points=maximum+3,hardest_entries=[v for v in entries if v['total_degree']==maximum])
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
