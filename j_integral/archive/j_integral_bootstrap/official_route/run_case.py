"""Stage launcher. Expensive numerical work stays in parameterized workers."""
import argparse,json,subprocess,sys,os,hashlib,socket,time
from pathlib import Path
PROJECT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent

def run(a):
    case=json.loads(Path(a.case).read_text())
    paths={k:str((PROJECT/case[k]).resolve()) for k in ['family','basis','targets','binary','reconstructor','work']}
    work=Path(paths['work']);work.mkdir(parents=True,exist_ok=True)
    if a.stage=='anchor':
        workers=a.workers or case['workers']
        if len(os.sched_getaffinity(0))<workers:raise RuntimeError('worker budget exceeds actual CPU affinity')
        cmd=[sys.executable,str(HERE/'cold_anchor.py'),'--work',str(work/'anchor'),'--workers',str(workers)]
        for k in ['family','basis','targets','binary','reconstructor']:cmd+=['--'+k,paths[k]]
        for k in ['initial','batch','max_points']:cmd+=['--'+k.replace('_','-'),str(case[k])]
    else:
        anchor=work/'anchor/connection.json'
        if not anchor.exists() or json.loads(anchor.read_text()).get('status')!='complete':raise RuntimeError('production requires complete first anchor')
        driver='run_restricted.py' if case.get('restricted_numeric',False) else 'run_upstream.py'
        cmd=[sys.executable,str(HERE/driver),'--work',str(work/'production'),'--anchor',str(anchor),'--points',str(work/'anchor/samples/points'),'--gauge',str(case['gauge'])]
        if case.get('restricted_numeric',False):cmd+=['--closure-from-bundle']
        for k in ['family','basis','targets','binary']:cmd+=['--'+k,paths[k]]
    output=work/('anchor/connection.json' if a.stage=='anchor' else 'production/connection.json')
    receipt=work/(a.stage+'-'+os.environ.get('SLURM_JOB_ID','local')+'.json')
    record=dict(stage=a.stage,case_sha256=hashlib.sha256(Path(a.case).read_bytes()).hexdigest(),
                input_sha256={k:hashlib.sha256(Path(paths[k]).read_bytes()).hexdigest() for k in ['family','basis','targets','binary']},
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),job_id=os.environ.get('SLURM_JOB_ID'),host=socket.gethostname(),pid=os.getpid(),affinity=sorted(os.sched_getaffinity(0)),start=time.time(),command=cmd)
    receipt.write_text(json.dumps(record,indent=2)+'\n')
    proc=subprocess.run(cmd,cwd=PROJECT)
    validation_code=0
    if proc.returncode==0 and a.stage=='production':
        validation=subprocess.run([sys.executable,str(HERE.parent/'tests/compare_connections.py'),str(work/'anchor/connection.json'),str(output),str(work/'production/anchor_exact_comparison.json')])
        validation_code=validation.returncode
    record.update(returncode=proc.returncode,validation_returncode=validation_code,end=time.time(),scientific_complete=proc.returncode==0 and validation_code==0 and output.exists() and json.loads(output.read_text()).get('status')=='complete')
    receipt.write_text(json.dumps(record,indent=2)+'\n')
    if not record['scientific_complete']:raise RuntimeError('stage incomplete; see '+str(receipt))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--stage',choices=['anchor','production'],required=True);p.add_argument('--workers',type=int)
    run(p.parse_args())
