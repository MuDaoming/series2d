import unittest
from reference_space import BootstrapSpec,RatFunc,build_reduction_ansatz
from fixed_fit import fixed_entry_fitter


class FixedFitTests(unittest.TestCase):
    def test_known_polynomial_does_not_consume_ibp_points(self):
        p=1000003
        spec=BootstrapSpec.from_inputs(family={},variable='x',assignments={},target=['known','unknown'],
            masters=['M'],de=[['0']],denominators=[[[1]],[[1]]],numerator_caps=[[20],[1]],prime=p)
        ansatz=build_reduction_ansatz(spec)
        known=RatFunc.make([3]+[0]*19+[1],[1],p)
        fitter,rank=fixed_entry_fitter(ansatz,{(0,0):known},list(range(1,22)))
        self.assertEqual(rank,21)
        self.assertEqual(len(fitter.samples),0)
        for x in [23,24]:fitter.add(dict(prime=p,x=x,values=[[known.evaluate(x)],[7+2*x]]))
        parameters=fitter.result()['parameters'];self.assertIsNotNone(parameters)
        fitted=ansatz.assemble(parameters)
        self.assertTrue((fitted[0][0]-known).is_zero())
        self.assertTrue((fitted[1][0]-RatFunc.make([7,2],[1],p)).is_zero())
        fitter.check_holdout(dict(prime=p,x=29,values=[[known.evaluate(29)],[65]]))


if __name__=='__main__':unittest.main()
