"""Resumable local sampling, with immutable inputs and per-point receipts."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from time import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(args):
    work=Path(args.work).resolve();work.mkdir(parents=True,exist_ok=True)
    for directory in ['input','points','logs','worker','receipts']:(work/directory).mkdir(exist_ok=True)
    sources={'family':args.family,'basis':args.basis,'binary':args.binary}
    if args.targets:sources['targets']=args.targets
    if args.derivatives:sources['derivatives']=args.derivatives
    identity=dict(mode=args.mode,files={k:sha(v) for k,v in sources.items()})
    ip=work/'identity.json'
    if ip.exists() and json.loads(ip.read_text())!=identity:raise ValueError('run identity changed')
    ip.write_text(json.dumps(identity,indent=2)+'\n')
    frozen={}
    for role,path in sources.items():
        target=work/('worker' if role=='binary' else 'input')/(Path(path).name if role=='binary' else role+'.json')
        if not target.exists():shutil.copy2(path,target)
        if sha(target)!=identity['files'][role]:raise ValueError('frozen asset differs')
        frozen[role]=target
    cfg=json.loads(frozen['family'].read_text());basis=json.loads(frozen['basis'].read_text())
    def validate(path,point):
        d=json.loads(path.read_text())
        if d['prime']!=cfg['prime'] or int(d['t_point'])!=point:raise ValueError('point identity differs')
        if args.mode in ('full','subset'):
            if d['masters_at_point']!=basis or d['unclosed_derivatives']:raise ValueError('basis or closure differs')
            expected=basis if args.mode=='full' else json.loads(frozen['derivatives'].read_text())
            if set(d['DE_values'])!=set(expected):raise ValueError('DE coverage incomplete')
            expected_objects=args.expected_objects
            if expected_objects is None:
                expected_objects=len(json.loads(frozen['targets'].read_text()))
                if cfg.get('include_fbi_boundaries',False):expected_objects+=2*d['FBI_rank']
            if len(d['reduction_values'])!=expected_objects:raise ValueError('target coverage incomplete')
        return d
    for source in args.seed:
        data=json.loads(Path(source).read_text());point=int(data['t_point']);target=work/'points'/f'point{point}.json'
        validate(Path(source),point)
        if not target.exists():shutil.copy2(source,target)
        elif sha(target)!=sha(source):
            old=json.loads(target.read_text())
            for section in ['DE_values',*(['reduction_values','masters_at_point'] if args.mode=='full' else [])]:
                if old[section]!=data[section]:raise ValueError('seed conflicts with existing point')
        (work/'receipts'/f'point{point}.json').write_text(json.dumps(dict(
            point=point,seed_source=str(Path(source).resolve()),output_sha256=sha(target),validated=True),indent=2)+'\n')
    def worker(point):
        output=work/'points'/f'point{point}.json';receipt=work/'receipts'/f'point{point}.json'
        with (work/'receipts'/f'point{point}.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            if receipt.exists():
                saved=json.loads(receipt.read_text())
                if saved.get('validated') and output.exists() and saved['output_sha256']==sha(output):return output
            if output.exists():
                validate(output,point)
                receipt.write_text(json.dumps(dict(point=point,adopted_existing=True,validated=True,output_sha256=sha(output)),indent=2)+'\n')
                return output
            temp=output.with_suffix('.tmp.json');cpu=work/'logs'/f'point{point}.cpu.txt'
            if args.mode in ('full','subset'):
                command=[str(frozen['binary']),str(frozen['family']),str(temp),'--point',str(point),
                         '--targets-file',str(frozen['targets']),'--basis-file',str(frozen['basis']),
                         '--values','1','--compact','1','--integration','atomic']
                if args.mode=='subset':command+=['--derivative-names-file',str(frozen['derivatives'])]
            else:command=[str(frozen['binary']),str(frozen['family']),str(frozen['basis']),str(point),str(temp)]
            started=time()
            with (work/'logs'/f'point{point}.log').open('w') as log:
                process=subprocess.Popen(['/usr/bin/time','-f','%U %S %M','-o',str(cpu),*command],stdout=log,stderr=subprocess.STDOUT)
                code=process.wait()
            record=dict(point=point,returncode=code,start=started,end=time(),wrapper_pid=process.pid,validated=False)
            if code==0:
                data=validate(temp,point);temp.replace(output)
                record.update(validated=True,output_sha256=sha(output),timings=data.get('timings'))
                fields=cpu.read_text().split();record.update(cpu_user_seconds=float(fields[0]),cpu_system_seconds=float(fields[1]),max_rss_KiB=int(fields[2]))
            receipt.write_text(json.dumps(record,indent=2)+'\n')
            if code:raise RuntimeError('worker failed: '+str(receipt))
            print(json.dumps(dict(point=point,mode=args.mode,wall_seconds=record['end']-started,
                                  cpu_seconds=record['cpu_user_seconds']+record['cpu_system_seconds'])),flush=True)
            return output
    points=list(range(args.first,args.first+args.count))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(worker,p):p for p in points}
        for future in as_completed(futures):future.result()
    manifest=dict(family=str(frozen['family']),points=[str(work/'points'/f'point{p}.json') for p in points])
    (work/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(completed=len(points),manifest=str(work/'manifest.json'))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['full','cut','subset']);p.add_argument('work');p.add_argument('--family',required=True)
    p.add_argument('--basis',required=True);p.add_argument('--binary',required=True);p.add_argument('--targets')
    p.add_argument('--derivatives')
    p.add_argument('--first',type=int,default=17);p.add_argument('--count',type=int,default=17)
    p.add_argument('--workers',type=int,default=2);p.add_argument('--expected-objects',type=int)
    p.add_argument('--seed',action='append',default=[])
    args=p.parse_args()
    if args.mode in ('full','subset') and not args.targets:p.error('full/subset mode needs targets')
    if args.mode=='subset' and not args.derivatives:p.error('subset mode needs derivative names')
    if args.workers<1 or args.count<1:p.error('positive workers and count required')
    main(args)
