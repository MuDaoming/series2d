"""Existing unmodified upstream lifecycle with the restricted numeric adapter."""
import argparse,os,json
from pathlib import Path
import run_upstream
from restricted_numeric import SectorRestrictedNumericBackend,closure_tasks

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['work','anchor','family','basis','targets','binary']:parser.add_argument('--'+key,required=True)
    parser.add_argument('--points',action='append',default=[])
    parser.add_argument('--gauge',type=int,default=3)
    parser.add_argument('--max-samples',type=int,default=512)
    parser.add_argument('--closure-from-bundle',action='store_true')
    args=parser.parse_args()
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','ANALYTIC_BOOTSTRAP_FIT_THREADS']:os.environ.setdefault(key,'1')
    def numeric(*pos,**kwargs):
        tasks=closure_tasks(json.loads((Path(args.work)/'bundle.json').read_text())) if args.closure_from_bundle else None
        return SectorRestrictedNumericBackend(*pos,closure_by_task=tasks,**kwargs)
    run_upstream.NumericBackend=numeric
    run_upstream.run(args)
