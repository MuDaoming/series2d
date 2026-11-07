"""A unique native prefix solution must not be solved a second time."""
import unittest
from unittest.mock import patch
from reference_space import BootstrapSpec,RatMat,RatFunc,build_pl_ansatz
from native_fit import joint_fit
class Reuse(unittest.TestCase):
 def test_full_identity_checks_without_second_elimination(self):
  p=1000003;n=12
  spec=BootstrapSpec.from_inputs(family={'test':'polynomial-independent-entries'},variable='x',assignments={},prime=p,masters=['M'+str(i) for i in range(n)],target=['T'+str(i) for i in range(n)],mode='reduction',de=RatMat({},(n,n),p),denominators=[[[1] for j in range(n)] for i in range(n)],numerator_caps=[[4]*n for i in range(n)])
  a=build_pl_ansatz(spec);self.assertEqual(a.dimension,720)
  parameters=list(range(1,a.dimension+1));matrix=a.assemble(parameters)
  samples=[dict(prime=p,x=x,values=a.evaluate(x,parameters)) for x in range(10,18)]
  fixed={(0,j):matrix[0][j] for j in range(n)}
  with patch('native_fit._native_conditions_eliminate',side_effect=AssertionError('duplicate solve')):
   result=joint_fit(a,samples,fixed,list(range(10,18)))
   self.assertTrue(result['compatible']);self.assertEqual(result['route'],'reuse-upstream-unique-solution');self.assertEqual(result['parameters'],parameters);self.assertEqual(result['fixed_identity_value_checks'],60)
   bad=dict(fixed);bad[0,0]=bad[0,0]+RatFunc.make([1],[1],p)
   self.assertFalse(joint_fit(a,samples,bad,list(range(10,18)))['compatible'])
   broken=[dict(s) for s in samples];broken[-1]=dict(broken[-1],values=[row[:] for row in broken[-1]['values']]);broken[-1]['values'][0][0]+=1
   self.assertFalse(joint_fit(a,broken,fixed,list(range(10,18)))['compatible'])
if __name__=='__main__':unittest.main()
