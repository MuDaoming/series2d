"""Isolated reproduction of the pending block, then original AB source-order tests."""
import sys,json,inspect,shutil,argparse,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
import bootstrap_global as bg
from analytic_bootstrap.pl import plan_pl_orders,build_pl_ansatz
from analytic_bootstrap.learning import learn_profiles
from native_fit import joint_fit
OUT=ROOT/'results/dp_official';OUT.mkdir(exist_ok=True)
class Captured(Exception):pass
capt={}
def capture(spec,**kw):
    capt.update(inspect.currentframe().f_back.f_locals);capt['spec']=spec
    raise Captured()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--extra',type=int,default=0);ap.add_argument('--profile',action='store_true');ap.add_argument('--official-fit',action='store_true');a=ap.parse_args()
    work=OUT/'capture';work.mkdir(exist_ok=True)
    shutil.copy(ROOT/'results/dp_complete_run/connection.json',work/'connection.json')
    bg.build_pl_ansatz=capture
    args=argparse.Namespace(work=str(work),family=str(ROOT.parent/'j_integral_global/examples/dp_physical.json'),basis=str(ROOT.parent/'j_integral_global/examples/dp496_basis.json'),points=[str(ROOT/'results'/s/'points') for s in ['dp_full','dp_extra_de']],diagonals=[str(ROOT/'results/dp_cuts/reconstruction.json'),str(ROOT/'results/dp_exceptional/diagonal33.json')],partial=[str(ROOT/'results/dp_driven/partial_33.json'),str(ROOT/'results/dp_extra_de/partial41.json')],backend='PL-structural',relax_all_poles=False,witness_budget=0,internal_holdouts=3,restart=False)
    try:bg.main(args)
    except Captured:pass
    c=capt;spec=c['spec'];start=time.perf_counter();report={'extra_order':a.extra,'block':[c['masters'][i] for i in c['block']],'ambient':sum(v+1 for row in c['caps'] for v in row),'plan':plan_pl_orders(spec,extra_order=a.extra)}
    def progress(d):print(json.dumps(d),flush=True)
    if a.profile:
        try:
            prof=learn_profiles(spec,progress=progress);report['profile_timing']=prof['timing'];(OUT/'profile.json').write_text(json.dumps(prof,indent=2))
        except Exception as e:report['profile_error']=repr(e)
    t=time.perf_counter();ansatz=build_pl_ansatz(spec,backend='PL-Hasse',laurent_mode='structural-refined-all',**({'order_profile':prof['order']} if a.profile and 'prof' in locals() else {'extra_order':a.extra}));report['construct_seconds']=time.perf_counter()-t;report['dimension']=ansatz.dimension
    fixed={(aa,b):c['learned'][i,j] for aa,i in enumerate(c['block']) for b,j in enumerate(c['sub']) if (i,j) in c['learned'] and (i,j) in c['active']}
    def sample(t):return dict(prime=c['p'],x=(t-c['center'])%c['p'],values=[[c['all_points'][t][c['masters'][i]].get(c['masters'][j],0) for j in c['sub']] for i in c['block']])
    fit=joint_fit(ansatz,[sample(t) for t in c['train']],fixed,[(t-c['center'])%c['p'] for t in c['times']]);report['fit']={k:v for k,v in fit.items() if k!='parameters'}
    (OUT/f'order{a.extra}.pre_fit.json').write_text(json.dumps(report,indent=2))
    if a.official_fit:
        from analytic_bootstrap.upstream.fitting import ProductionFitter
        fitter=ProductionFitter(ansatz,backend='native');t=time.perf_counter()
        try:
            for x in c['train']:
                result=fitter.add(sample(x))
                if result['parameters'] is not None:break
            report['official_fit']={k:v for k,v in result.items() if k!='parameters'}
        except Exception as e:report['official_fit_error']=repr(e)
        finally:fitter.close()
        report['official_fit_seconds']=time.perf_counter()-t
    report['seconds']=time.perf_counter()-start
    (OUT/f'order{a.extra}.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
if __name__=='__main__':main()
