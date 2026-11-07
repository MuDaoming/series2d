import json,sys,tempfile,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'j_integral_bootstrap/tests'));from rational_learning import make_learner
L=make_learner(1000003,200001);O=ROOT/'j_integral_bootstrap/build';D=L.one
for i in range(40):D*=L.poly([300+i,1])
with tempfile.TemporaryDirectory(prefix='witness_only_regression_') as temp:
 T=Path(temp)
 for zero_witness in [False,True]:
  N=L.poly([2,1]) if not zero_witness else L.poly([-(17-200001)%1000003,1]);H=L.poly([5,3,7]);E=L.poly([23,1]);true_num=H;true_den=D*E
  witness={'prime':1000003,'t0':'2/5','masters':['master'],'de':{'J_target':{'rhs':{'master':{'numerator_tau':[int(c) for c in N],'denominator_tau':[int(c) for c in D]}}}}};(T/'witness.json').write_text(json.dumps(witness));files=[]
  for t in range(17,58):
   x=(t-200001)%1000003;v=int(H(x))*pow(int(true_den(x)),-1,1000003)%1000003;path=T/f'point{t}.json';path.write_text(json.dumps({'prime':1000003,'masters_at_point':['master'],'t_point':str(t),'reduction_values':{'J_target':{'master':v}}}));files.append(str(path))
  manifest={'prime':1000003,'t0':'2/5','masters':['master'],'section':'reduction_values','points':files,'known':[],'holdout':8,'ratios':True,'ratio_budget':8,'witness_only':[{'path':str(T/'witness.json'),'section':'de'}]};(T/'manifest.json').write_text(json.dumps(manifest));subprocess.run([str(O/'preconditioned'),str(T/'manifest.json'),str(T/'out.json')],check=True,capture_output=True)
  out=json.load(open(T/'out.json'));assert out['new_count']==1 and out['witness_only_coefficients']==1
  c=out['reductions']['J_target']['master'];assert not (L.poly(c['numerator_tau'])*true_den-H*L.poly(c['denominator_tau']))
  print('same-label witness is not known R; witness-zero='+str(zero_witness),'PASS')
  manifest['witness_only']=[];(T/'no_witness.json').write_text(json.dumps(manifest));subprocess.run([str(O/'preconditioned'),str(T/'no_witness.json'),str(T/'without.json')],check=True,capture_output=True);assert not json.load(open(T/'without.json'))['new_count'];print('witness supplies candidates only, no witness gives no reconstruction PASS')
  witness['t0']='3/5';(T/'bad_witness.json').write_text(json.dumps(witness));manifest['witness_only']=[{'path':str(T/'bad_witness.json'),'section':'de'}];(T/'bad_manifest.json').write_text(json.dumps(manifest));p=subprocess.run([str(O/'preconditioned'),str(T/'bad_manifest.json'),str(T/'bad_out.json')],capture_output=True);assert p.returncode!=0;print('witness tau identity mismatch rejected PASS')
print('No blind verification data read.')
