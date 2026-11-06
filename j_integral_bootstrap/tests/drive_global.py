"""Advance full J bootstrap when new completed numerical points arrive."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent


def main(args):
    work=Path(args.work).resolve();work.mkdir(parents=True,exist_ok=True)
    pointdir=Path(args.full).resolve();last=-1
    while True:
        paths=sorted((p for p in pointdir.glob('point*.json') if '.tmp.' not in p.name),
                     key=lambda p:int(p.stem.replace('point','')))
        count=len(paths)
        if count<9 or count==last:
            time.sleep(5);continue
        if last>=0 and count<last+2 and count<args.max_points:
            time.sleep(5);continue
        manifest=work/f'pilot_{count}.json';partial=work/f'partial_{count}.json'
        manifest.write_text(json.dumps(dict(family=str(Path(args.family).resolve()),
                                           points=[str(p) for p in paths]),indent=2)+'\n')
        with (work/f'partial_{count}.log').open('w') as log:
            code=subprocess.run([args.reconstructor,str(manifest),str(partial)],stdout=log,stderr=subprocess.STDOUT).returncode
        if code not in (0,2):raise RuntimeError('partial learning failed')
        command=[sys.executable,str(HERE/'bootstrap_global.py'),str(work/'bootstrap'),
                 '--family',args.family,'--basis',args.basis,'--points',args.full,
                 '--points',args.extra_points,'--diagonals',args.cut_diagonal,
                 '--diagonals',args.extra_diagonal,'--partial',str(partial),'--partial',args.extra_partial]
        started=time.time()
        with (work/f'bootstrap_{count}.log').open('w') as log:
            code=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT).returncode
        checkpoint=json.loads((work/'bootstrap/connection.json').read_text())
        event=dict(point_count=count,returncode=code,status=checkpoint['status'],
                   known_rows=len(checkpoint['known_rows']),blocks=len(checkpoint['blocks']),
                   attempt_seconds=time.time()-started)
        with (work/'events.jsonl').open('a') as f:f.write(json.dumps(event)+'\n')
        print(json.dumps(event),flush=True)
        if code==0:return 0
        if code!=2:raise RuntimeError('bootstrap error: '+str(work/f'bootstrap_{count}.log'))
        if count>=args.max_points:return 2
        last=count


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('work')
    for name in ['family','basis','full','extra-points','cut-diagonal','extra-diagonal','extra-partial','reconstructor']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--max-points',type=int,default=33);raise SystemExit(main(p.parse_args()))
