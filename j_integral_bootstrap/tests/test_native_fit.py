import unittest
from reference_space import BootstrapSpec,RatFunc,build_reduction_ansatz,build_de_ansatz
from native_fit import joint_fit


class NativeFitTests(unittest.TestCase):
    def test_identity_nodes_are_distinct_in_the_field(self):
        p=101
        spec=BootstrapSpec.from_inputs(family={},variable='x',assignments={},target=['f'],
            masters=['M'],de=[['0']],denominators=[[[1]]],numerator_caps=[[1]],prime=p)
        ansatz=build_reduction_ansatz(spec)
        samples=[dict(prime=p,x=1,values=[[1]])]
        with self.assertRaisesRegex(ValueError,'insufficient points'):
            joint_fit(ansatz,samples,{(0,0):RatFunc.make([0],[1],p)},[0,p])
        known=RatFunc.make([0,1],[1],p)
        fitted=joint_fit(ansatz,samples,{(0,0):known},[0,p,1,p+1])
        self.assertTrue(fitted['compatible'])
        self.assertIsNotNone(fitted['parameters'])
        self.assertEqual(fitted['fixed_equation_count'],2)
        self.assertEqual(ansatz.evaluate(2,fitted['parameters']),[[2]])

    def test_shared_de_extension_with_one_fixed_entry(self):
        p=1000003
        spec=BootstrapSpec.from_inputs(family={},variable='x',assignments={},target=None,
            masters=['S','T1','T2'],de={'B':[['0']],'A':[['1/x','0'],['2/x','0']]},
            denominators=[['x'],['x']],numerator_caps=[[0],[0]],prime=p,mode='de')
        ansatz=build_de_ansatz(spec);self.assertEqual(ansatz.dimension,1)
        known=RatFunc.make([7],[0,1],p)
        result=joint_fit(ansatz,[],{(0,0):known},[2])
        self.assertIsNotNone(result['parameters'])
        fitted=ansatz.assemble(result['parameters'])
        self.assertTrue((fitted[0][0]-known).is_zero())
        self.assertTrue((fitted[1][0]-2*known).is_zero())

    def test_fixed_identity_and_unknown_row(self):
        p=1000003
        spec=BootstrapSpec.from_inputs(family={},variable='x',assignments={},target=['known','unknown'],
            masters=['M'],de=[['0']],denominators=[[[1]],[[1]]],numerator_caps=[[20],[1]],prime=p)
        ansatz=build_reduction_ansatz(spec);known=RatFunc.make([3]+[0]*19+[1],[1],p)
        samples=[dict(prime=p,x=x,values=[[known.evaluate(x)],[7+2*x]]) for x in [23,24]]
        result=joint_fit(ansatz,samples,{(0,0):known},list(range(1,22)))
        self.assertTrue(result['compatible']);self.assertIsNotNone(result['parameters'])
        rows=ansatz.assemble(result['parameters'])
        self.assertTrue((rows[0][0]-known).is_zero())
        self.assertTrue((rows[1][0]-RatFunc.make([7,2],[1],p)).is_zero())
        bad=dict(prime=p,x=25,values=[[known.evaluate(25)],[1]])
        self.assertFalse(joint_fit(ansatz,samples+[bad],{(0,0):known},list(range(1,22)))['compatible'])


if __name__=='__main__':unittest.main()
