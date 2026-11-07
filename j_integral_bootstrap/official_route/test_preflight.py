"""DP regression: reject the archived failed block; retain an actual solved block."""
import sys,json,argparse,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import probe_orders as probe
from reference_space import BootstrapSpec,RatFunc,RatMat,build_pl_ansatz
from preflight import preflight
OUT=ROOT/'results/dp_official';P=1000003;CENTER=200001

def main():
    cfg=ROOT.parent/'j_integral_global/examples'
    args=argparse.Namespace(work=str(OUT/'capture'),family=str(cfg/'dp_physical.json'),basis=str(cfg/'dp496_basis.json'),points=[str(ROOT/'results'/s/'points') for s in ['dp_full','dp_extra_de']],diagonals=[str(ROOT/'results/dp_cuts/reconstruction.json'),str(ROOT/'results/dp_exceptional/diagonal33.json')],partial=[str(ROOT/'results/dp_driven/partial_33.json'),str(ROOT/'results/dp_extra_de/partial41.json')],backend='PL-structural',relax_all_poles=False,witness_budget=0,internal_holdouts=3,restart=False)
    probe.bg.build_pl_ansatz=probe.capture
    try:probe.bg.main(args)
    except probe.Captured:pass
    c=probe.capt;failing=build_pl_ansatz(c['spec'],backend='PL-Hasse',laurent_mode='structural-refined-all')
    def samples(top,sub,times):
        return [dict(prime=P,x=(t-CENTER)%P,values=[[c['all_points'][t][m].get(j,0) for j in sub] for m in top]) for t in times if all(m in c['all_points'][t] for m in top)]
    lib=ROOT/'build/upstream';report={}
    report['failure_prefix']=preflight(failing,samples([c['masters'][i] for i in c['block']],[c['masters'][i] for i in c['sub']],c['train']),max_points=3,kernel_directory=lib)
    report['failure']=preflight(failing,samples([c['masters'][i] for i in c['block']],[c['masters'][i] for i in c['sub']],c['train']),max_points=38,kernel_directory=lib)
    checkpoint=json.loads((OUT/'capture/connection.json').read_text());de=checkpoint['de']
    previous=next(b for b in checkpoint['blocks'] if b.get('route')=='analyticity-bootstrap');top=previous['top'];deps=set()
    def follow(m):
        for j in de[m]['rhs']:
            if j not in top and j not in deps:deps.add(j);follow(j)
    for m in top:follow(m)
    sub=[m for m in checkpoint['known_rows'] if m in deps];zero=RatFunc.make([0],[1],P)
    def rf(r):return RatFunc.make(r['numerator_tau'],r['denominator_tau'],P)
    def matrix(names):return {(i,j):rf(de[m]['rhs'][n]) for i,m in enumerate(names) for j,n in enumerate(names) if n in de[m]['rhs']}
    A=matrix(top);B=matrix(sub);gauge=previous['gauge'];shift=RatFunc.make([400001*gauge%P,2*gauge],[200001*200000%P,400001,1],P)
    # endpoint polynomial is (tau+center)(tau+center-1).
    for M,n in [(A,len(top)),(B,len(sub))]:
        for i in range(n):M[i,i]=M.get((i,i),zero)+shift
    den=[];caps=[]
    for m in top:
        dr=[];cr=[]
        for j in sub:
            r=de[m]['rhs'].get(j);dr.append(r['denominator_tau'] if r else [1]);cr.append(len(r['numerator_tau'])-1 if r else -1)
        den.append(dr);caps.append(cr)
    spec=BootstrapSpec.from_inputs(family={'source':'actual-solved-DP-block'},variable='tau',assignments={},masters=sub+top,target=None,de={'A':RatMat(A,(len(top),len(top)),P),'B':RatMat(B,(len(sub),len(sub)),P)},denominators=den,numerator_caps=caps,prime=P,mode='de')
    ansatz=build_pl_ansatz(spec,backend='PL-Hasse',laurent_mode='structural-refined-all')
    report['success_block']=top;report['success_dimension']=ansatz.dimension
    report['success_prefix']=preflight(ansatz,samples(top,sub,c['train']),max_points=3,kernel_directory=lib)
    report['success_all']=preflight(ansatz,samples(top,sub,c['train']),max_points=len(c['train']),kernel_directory=lib)
    report['missing_native']=preflight(ansatz,samples(top,sub,c['train']),kernel_directory=OUT/'deliberately_missing')
    assert report['missing_native']['status']=='undetermined'
    assert report['failure']['status']=='incompatible',report
    assert report['success_prefix']['status']!='incompatible',report
    assert report['success_all']['status']=='compatible',report
    (OUT/'preflight_regression.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
if __name__=='__main__':main()
