"""Remove the explicitly declared common gauge; no reconstruction or new samples."""
import argparse,json,gzip
from pathlib import Path
from fractions import Fraction
from analytic_bootstrap import RatFunc
from analytic_bootstrap.upstream import SolvedDE
p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--run',required=True);a=p.parse_args();case=Path(a.case).resolve();run=Path(a.run).resolve();norm=json.loads((case/'input/normalization.json').read_text());solved=SolvedDE.from_dict(json.loads((run/'solution.json').read_text()));prime=solved.context['prime'];q=Fraction(norm['t0']);center=q.numerator*pow(q.denominator,-1,prime)%prime;gauge=norm['gauge_power'];shift=RatFunc.make([gauge],[center,1],prime)+RatFunc.make([gauge],[(center-1)%prime,1],prime);zero=RatFunc.make([0],[1],prime);names=solved.context['masters'];de={m:{'rhs':{}} for m in norm['physical_basis']}
for i,j in set(solved.matrix.entries)|{(i,i) for i in range(len(names))}:
    value=solved.matrix.entries.get((i,j),zero)-(shift if i==j else zero)
    if not value.is_zero():de[names[i]]['rhs'][names[j]]=dict(numerator_tau=value.num,denominator_tau=value.den)
result=dict(status='complete',prime=prime,t0=norm['t0'],masters=norm['physical_basis'],de=de,source='original AB upstream CLI solve-de; common scalar normalization removed')
(run/'physical_de.json.gz').write_bytes(gzip.compress((json.dumps(result,separators=(',',':'))+'\n').encode(),mtime=0))
print('physical DE exported',len(names),'masters',sum(len(row['rhs']) for row in de.values()),'nonzero entries')
