"""Compare successful scalar-cut diagonal blocks against a full numerical DE."""
import argparse
import json
from pathlib import Path


def main(args):
    cut=json.loads(Path(args.cut).read_text());full=json.loads(Path(args.full).read_text())
    if cut['prime']!=full['prime'] or cut['t_point']!=full['t_point']:
        raise ValueError('field or point differs')
    checks=0;bad=[];failures=[]
    for block in cut['blocks']:
        if block['status']!='ok':
            failures.append(dict(sector=block['sector'],error=block['error']))
            continue
        for m in block['masters']:
            for n in block['masters']:
                checks+=1
                actual=cut['DE_values'][m].get(n,0);expected=full['DE_values'][m].get(n,0)
                if actual!=expected:bad.append(dict(row=m,column=n,actual=actual,expected=expected))
    result=dict(prime=cut['prime'],point=cut['t_point'],
                total_J_sectors=len(cut['blocks']),successful=len(cut['blocks'])-len(failures),
                unresolved=failures,scalar_checks=checks,mismatches=bad,
                cut_seconds=cut['total_seconds'],all_successful_blocks_match=not bad,
                complete=not bad and not failures)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('unresolved','mismatches')}))
    if bad:raise AssertionError('cut diagonal disagrees with full DE')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cut');p.add_argument('full');p.add_argument('output');main(p.parse_args())
