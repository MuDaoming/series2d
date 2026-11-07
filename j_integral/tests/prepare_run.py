"""Write argv/config and provenance only; AB is executed by its original CLI."""
import argparse,json,sys,os,hashlib,time,socket
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--run',required=True);a=p.parse_args();case=Path(a.case).resolve();run=Path(a.run).resolve();run.mkdir(parents=True,exist_ok=False);tests=Path(__file__).resolve().parent;inp=case/'input';worker=tests/'backend_worker.py';binary=tests.parent/'src/IBP/build/j-ibp'
config=dict(anchor_command=[sys.executable,str(worker),'anchor','--input',str(inp/'anchor.json.gz')],cut_command=[sys.executable,str(worker),'cut','--input',str(inp/'diagonals.json')],numeric_command=[sys.executable,str(worker),'numeric','--input',str(inp/'family_request.json'),'--normalization',str(inp/'normalization.json'),'--binary',str(binary),'--cache',str(run/'numeric')],learn_profiles=True,sampling=dict(max_samples=512,holdouts=2,shared_de_samples=True,global_holdouts=2,fit_backend='stream'))
(run/'backend.json').write_text(json.dumps(config,indent=2)+'\n')
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
import analytic_bootstrap
pkg=Path(analytic_bootstrap.__file__).resolve().parent
files=[path for path in inp.iterdir() if path.is_file()]+[worker,binary,pkg/'upstream/__main__.py',pkg/'upstream/pipeline.py',pkg/'upstream/backends.py']
record=dict(start=time.time(),host=socket.gethostname(),pid=os.getpid(),case=case.name,python=sys.executable,AB_package=str(pkg),input_and_code_sha256={str(path):sha(path) for path in files},environment={k:v for k,v in os.environ.items() if k.startswith('ANALYTIC_BOOTSTRAP_') or k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')},entry='python -m analytic_bootstrap.upstream learn; python -m analytic_bootstrap.upstream solve-de',new_numeric_cache=True)
(run/'run.json').write_text(json.dumps(record,indent=2)+'\n')
print(run,flush=True)
