"""Small correctness counterexamples; no production data or fresh points read."""
import sys,json,argparse,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
OUT=ROOT/'results/dp_diagnosis/review';OUT.mkdir(exist_ok=True)
from reference_space import BootstrapSpec,RatMat,RatFunc,build_pl_ansatz
from native_fit import joint_fit
import bootstrap_reductions as br
p=101;spec=BootstrapSpec.from_inputs(family={'name':'duplicate-identity'},variable='x',assignments={},masters=['m'],target=['f'],de=RatMat({},(1,1),p),denominators=[[[1]]],numerator_caps=[[1]],prime=p,mode='reduction');a=build_pl_ansatz(spec,backend='PL-Hasse');fit=joint_fit(a,[{'prime':p,'x':1,'values':[[1]]}],{(0,0):RatFunc.make([0],[1],p)},[0,p]);results={'duplicate_identity':{'compatible':fit['compatible'],'unique':fit['parameters'] is not None,'at_2':a.evaluate(2,fit['parameters']) if fit['parameters']is not None else None,'specified_fixed_function':0}}
# Cardinality equality must not stand in for target-set completion.
case=OUT/'cardinality';(case/'points').mkdir(parents=True,exist_ok=True)
m0,m1,target='J{1,0}','J{0,1}','J{1,1}'
connection={'prime':p,'t0':'0','masters':[m0,m1],'known_rows':[m0],'status':'running','de':{m0:{'rhs':{}},m1:{'rhs':{}}}}
point={'prime':p,'t_point':2,'masters_at_point':[m0,m1],'unclosed_derivatives':[],'DE_values':{},'reduction_values':{m0:{m0:1},target:{m1:1}}}
(case/'connection.json').write_text(json.dumps(connection));(case/'points/point2.json').write_text(json.dumps(point));br.main(argparse.Namespace(work=str(case/'work'),connection=str(case/'connection.json'),points=[str(case/'points')],partial=[],witness_budget=0,internal_holdouts=3));saved=json.loads((case/'work/reduction.json').read_text());results['false_complete_cardinality']={'status':saved['status'],'tested_objects':saved['tested_objects'],'known_targets':saved['known_targets'],'target_reduction_present':target in saved['reductions']}
# Empty verification manifests must not produce a success certificate.
manifest=OUT/'empty_manifest.json';manifest.write_text('{"points":[]}')
for label,binary,answer in [('DE','verify_connection',{'prime':p,'t0':'0','masters':[m0],'de':{}}),('R','verify_reductions',{'prime':p,'t0':'0','masters':[m0],'tested_objects':[],'reductions':{}})]:
 path=OUT/f'empty_{label}.json';path.write_text(json.dumps(answer));output=OUT/f'empty_{label}_result.json';run=subprocess.run([str(ROOT/'build'/binary),str(path),str(manifest),str(output)],capture_output=True,text=True);results['empty_verification_'+label]={'exit_code':run.returncode,'output':json.loads(output.read_text()) if output.exists() else run.stderr}
(OUT/'counterexamples.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
