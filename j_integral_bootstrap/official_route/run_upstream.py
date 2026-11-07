"""Use upstream learn_anchor -> prepare_setting -> solve_de unchanged.

A complete single-prime anchor must first be obtained by the numeric J backend.
It is used only for upstream learning and diagonal input, never as a fitted C.
"""
import argparse,json,sys,fcntl,os,ctypes
from pathlib import Path
from time import perf_counter
from upstream_j import AnchorBackend,DiagonalBackend,NumericBackend,definition,physical_connection,ROOT
from analytic_bootstrap.upstream import learn_anchor,prepare_setting,LearnedBundle
from analytic_bootstrap.storage import digest
from analytic_bootstrap.upstream.model import matrix_record


def atomic(path,data):
    temp=path.with_suffix('.tmp.json');temp.write_text(json.dumps(data,separators=(',',':'))+'\n');temp.replace(path)


def run(a):
    os.environ.setdefault('ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR',str(ROOT/'build/upstream'))
    native=Path(os.environ['ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR'])
    for name in ['upstream_evaluate.so','upstream_fit.so']:ctypes.CDLL(str(native/name))
    work=Path(a.work);work.mkdir(parents=True,exist_ok=True)
    with (work/'run.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        answer=json.loads(Path(a.anchor).read_text());family=definition(answer)
        setting=dict(prime=answer['prime'],assignments={});start=perf_counter()
        identity=dict(anchor=digest(answer),family=digest(family),gauge=a.gauge,upstream='0.3.0')
        ip=work/'identity.json'
        if ip.exists() and json.loads(ip.read_text())!=identity:raise ValueError('official run identity changed')
        atomic(ip,identity)
        anchor=AnchorBackend(answer,a.gauge);cut=DiagonalBackend(anchor,family)
        bundle_path=work/'bundle.json'
        if bundle_path.exists():bundle=LearnedBundle(json.loads(bundle_path.read_text()))
        else:
            bundle=learn_anchor(family,setting,anchor,learn_profiles=True)
            atomic(bundle_path,bundle.to_dict())
        # The numeric object has no reference to anchor or any off-diagonal truth.
        numeric=NumericBackend(a.family,a.basis,a.targets,a.binary,work/'numeric',point_directories=a.points,gauge=a.gauge)
        with prepare_setting(bundle,setting,cut,numeric,max_samples=a.max_samples,holdouts=2,
                             shared_de_samples=True,global_holdouts=2,fit_backend='stream') as point:
            for k in range(1,len(family['sectors'])+1):
                solved=point.solve_de(_internal=True,_sector_limit=k)
                # Save official validation records after each successfully solved sector.
                atomic(work/'progress.json',dict(status='running',completed=k,total=len(family['sectors']),seconds=perf_counter()-start,checks=solved.checks,native_new_points=numeric.calls))
                print(json.dumps(dict(completed=k,total=len(family['sectors']),seconds=perf_counter()-start)),flush=True)
            atomic(work/'solved_upstream.json',point.solve_de().to_dict())
            names=point.solve_de().context['masters']
            matrix=point.solve_de().matrix;de=physical_connection(matrix,names,answer['t0'],a.gauge)
            # Return in the numerical worker's fixed basis order for the verifier.
            atomic(work/'connection.json',dict(status='complete',prime=answer['prime'],t0=answer['t0'],masters=answer['masters'],known_rows=answer['masters'],de=de,checks=point.checks,global_check=point.solve_de().global_check,method='unmodified upstream learning and bottom-up production',normalization_power=a.gauge,native_new_points=numeric.calls,seconds=perf_counter()-start))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['work','anchor','family','basis','targets','binary']:p.add_argument('--'+key,required=True)
    p.add_argument('--points',action='append',default=[]);p.add_argument('--gauge',type=int,default=3);p.add_argument('--max-samples',type=int,default=512)
    a=p.parse_args()
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','ANALYTIC_BOOTSTRAP_FIT_THREADS']:os.environ.setdefault(key,'1')
    run(a)
