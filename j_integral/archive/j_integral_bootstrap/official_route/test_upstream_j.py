import json,tempfile,unittest
from pathlib import Path
from upstream_j import *
class Adapter(unittest.TestCase):
 def test_normalization_and_missing_diagonal(self):
  p=101;shift=normalization_shift(p,'2/5',3)
  mat=RatMat({},(1,1),p)
  de=physical_connection(mat,['J'],'2/5',3)
  got=RatFunc.make(**dict(num=de['J']['rhs']['J']['numerator_tau'],den=de['J']['rhs']['J']['denominator_tau'],p=p))
  self.assertTrue((got+shift).is_zero())
 def test_upstream_canonical_sector_permutation(self):
  answer=dict(masters=['J0','J1','J2'],prime=101,t0='2/5',de={'J0':{'rhs':{}},'J1':{'rhs':{'J0':dict(numerator_tau=[1],denominator_tau=[1])}},'J2':{'rhs':{}}})
  family=definition(answer);actual=[m for sector in family['sectors'] for m in sector['masters']]
  self.assertEqual(actual,['J0','J2','J1'])
  matrix=RatMat({(2,0):RatFunc.make([7],[1],101)},(3,3),101)
  exported=physical_connection(matrix,actual,'2/5',0)
  self.assertEqual(exported['J1']['rhs']['J0']['numerator_tau'],[7]);self.assertEqual(exported['J2']['rhs'],{})
 def test_numeric_protocol_and_basis_order(self):
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);names=['J0','J1']
   for name,data in [('family',dict(prime=101,t0='2/5')),('basis',names),('targets',[[1],[2]])]:
    (d/(name+'.json')).write_text(json.dumps(data))
   (d/'cache').mkdir();(d/'cache/point7.json').write_text(json.dumps(dict(prime=101,t_point='7',masters_at_point=names,unclosed_derivatives=[],DE_values={'J0':{'J0':2,'J1':3},'J1':{'J1':4}})))
   b=NumericBackend(d/'family.json',d/'basis.json',d/'targets.json','/bin/true',d/'work',point_directories=[d/'cache'],gauge=0)
   req=dict(mode='de',prime=101,masters=list(reversed(names)),targets=[dict(kind='derivative',variable='tau',master='J0')])
   r=b.prepare(req).evaluate((7-b.center)%101)
   self.assertEqual(r['values'],[[3,2]]);self.assertEqual(r['request_sha256'],digest(req));self.assertEqual(b.calls,0)
   (d/'cache/point7.json').write_text(json.dumps(dict(prime=103,t_point='7',masters_at_point=names,unclosed_derivatives=[],DE_values={'J0':{},'J1':{}})))
   with self.assertRaisesRegex(ValueError,'field/basis'):b.point((7-b.center)%101)
if __name__=='__main__':unittest.main()
