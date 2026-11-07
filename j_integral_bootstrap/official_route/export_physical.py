"""Label an official solved matrix using its certified context master ordering."""
import argparse,json
from pathlib import Path
from upstream_j import physical_connection
from analytic_bootstrap.upstream import SolvedDE
p=argparse.ArgumentParser();p.add_argument('solved');p.add_argument('output');p.add_argument('--basis',required=True);p.add_argument('--t0',required=True);p.add_argument('--gauge',required=True,type=int);a=p.parse_args()
s=SolvedDE.from_dict(json.loads(Path(a.solved).read_text()));actual=s.context['masters'];basis=json.loads(Path(a.basis).read_text())
if len(set(basis))!=len(basis) or len(basis)!=len(actual) or set(basis)!=set(actual):raise ValueError('export basis is not a duplicate-free bijection')
de=physical_connection(s.matrix,actual,a.t0,a.gauge)
out=dict(status='complete',prime=s.context['prime'],t0=a.t0,masters=basis,known_rows=basis,de=de,checks=s.checks,global_check=s.global_check,normalization_power=a.gauge,method='unmodified upstream production; physical export uses certified context order')
Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
