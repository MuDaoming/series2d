"""Obtain the first full rational J anchor from numerical IBP, without AB.

The upstream learner requires this one-time input. All later DE construction
is delegated to upstream. Independent point workers are dispatched by the
existing resumable collector, with an explicit worker budget.
"""
import argparse,json,subprocess,sys,fcntl
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main(a):
    work=Path(a.work).resolve();work.mkdir(parents=True,exist_ok=True)
    with (work/'anchor.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        for count in range(a.initial,a.max_points+1,a.batch):
            command=[sys.executable,str(ROOT/'tests/collect_points.py'),'full',str(work/'samples'),'--family',a.family,'--basis',a.basis,'--targets',a.targets,'--binary',a.binary,'--first',str(a.first),'--count',str(count),'--workers',str(a.workers)]
            subprocess.run(command,check=True)
            partial=work/'partial.json'
            proc=subprocess.run([a.reconstructor,str(work/'samples/manifest.json'),str(partial)])
            if proc.returncode not in (0,2):raise RuntimeError('rational reconstruction failed')
            data=json.loads(partial.read_text());masters=json.loads(Path(a.basis).read_text())
            status=dict(points=count,missing=len(data['unreconstructed']),reconstructed=data['reconstructed_count'])
            (work/'progress.json').write_text(json.dumps(status)+'\n');print(json.dumps(status),flush=True)
            if data['unreconstructed']:continue
            data.update(status='complete',masters=masters,known_rows=masters,method='first anchor from numeric J-IBP; univariate rational reconstruction',sample_count=count)
            for m in masters:data['de'].setdefault(m,{'rhs':{}})
            temp=work/'connection.tmp.json';temp.write_text(json.dumps(data,separators=(',',':'))+'\n');temp.replace(work/'connection.json')
            return 0
        raise RuntimeError('anchor point budget exhausted; partial answers retained, not complete')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['work','family','basis','targets','binary','reconstructor']:p.add_argument('--'+k,required=True)
    p.add_argument('--first',type=int,default=17);p.add_argument('--initial',type=int,default=64);p.add_argument('--batch',type=int,default=32);p.add_argument('--max-points',type=int,default=384);p.add_argument('--workers',type=int,default=1)
    a=p.parse_args()
    if a.workers<1 or a.initial<=8 or a.batch<1 or a.max_points<a.initial:p.error('invalid point/worker budget')
    raise SystemExit(main(a))
