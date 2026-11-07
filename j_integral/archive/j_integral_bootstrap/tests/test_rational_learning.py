"""Exact rational witnesses for numerator learning and local pole-order changes."""
import unittest
from rational_learning import make_learner
from reference_space import RatFunc


class RationalLearningTests(unittest.TestCase):
    def setUp(self):
        self.p=1000003;self.L=make_learner(self.p,0);self.times=list(range(17,37))

    def values(self,numerator,denominator):
        return [int(numerator(t))*pow(int(denominator(t)),-1,self.p)%self.p for t in self.times]

    def same(self,got,numerator,denominator):
        self.assertIsNotNone(got)
        wanted=RatFunc.make([int(c) for c in numerator],[int(c) for c in denominator],self.p)
        self.assertTrue((got-wanted).is_zero())

    def test_low_numerator_over_large_known_denominator(self):
        n=self.L.poly([5,3]);d=self.L.poly([1,1])**14
        values=self.values(n,d);self.same(self.L.fit(d,self.times,values),n,d)
        values[-1]=(values[-1]+1)%self.p
        self.assertIsNone(self.L.fit(d,self.times,values))

    def test_remove_unneeded_factors(self):
        n=self.L.poly([2,1]);d=self.L.poly([5,3,1])**2
        over=d*self.L.poly([1,1,0,1])**6;values=self.values(n,d)
        self.assertIsNone(self.L.fit(over,self.times,values))
        self.same(self.L.prune(over,self.times,values),n,d)

    def test_raise_only_the_required_pole_order(self):
        n=self.L.poly([3,1]);factor=self.L.poly([-7,1]);d=factor**3
        self.same(self.L.boost(factor,self.times,self.values(n,d)),n,d)


if __name__=='__main__':unittest.main()
