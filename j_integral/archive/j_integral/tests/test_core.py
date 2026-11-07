import copy
import json
from pathlib import Path
import sys
import unittest
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from periods import FixedPeriodReducer
from family import TwoLoopFamily
from series_engine import Model, generate_series
from relations import search_relations


class CoreTests(unittest.TestCase):
    def test_singular_pole_map_retains_obstruction(self):
        r = FixedPeriodReducer(1000003, [['1/z']], endpoints={0:'regulated_zero',1:'regular'})
        coordinates, _ = r.reduce(['1/z**2'])
        self.assertEqual(coordinates.get(('period',2,0,0)), 1)

    def test_matrix_connection_and_boundary_direction(self):
        r = FixedPeriodReducer(1000003,
            [['1/(3*z)','1/(5*(z+1))'],['0','2/(7*z)']],
            endpoints={0:'regulated_zero',1:'regular'})
        g = [r.coerce('z/(z+1)**2'),r.coerce('1/(z+1)')]
        coordinates, _ = r.reduce(r.covariant(g))
        self.assertFalse(any(key[0] == 'period' for key in coordinates))
        self.assertEqual(coordinates[('endpoint',1,0)],pow(4,-1,r.p))
        self.assertEqual(coordinates[('endpoint',1,1)],pow(2,-1,r.p))

    def test_nonfuchsian_is_rejected(self):
        with self.assertRaises(ValueError):
            FixedPeriodReducer(1000003,[['1/z**2']],endpoints={0:'regulated_zero',1:'regular'})

    def test_nonregular_endpoint_is_not_silently_zero(self):
        with self.assertRaises(ValueError):
            FixedPeriodReducer(1000003,[['1/z']],endpoints={0:'regular',1:'regular'})

    def test_vac_fbi_normalization_and_independent_t_coefficients(self):
        cfg=json.loads((ROOT/'examples/vac.json').read_text())
        fam=TwoLoopFamily(cfg)
        model=Model(fam.connection_model())
        self.assertTrue(model.flatness_check())
        result=generate_series(model,4)
        t,z,ta=s.symbols('t z tau')
        d=s.Rational(cfg['dimension']);t0=s.Rational(cfg['t0']);h=d/2
        U=(1-z)*(z+(1-z)*t*(1-t))
        # Independent generalized-binomial construction, not the connection recurrence.
        A=s.cancel(s.diff(U,t).subs(t,t0)/U.subs(t,t0))
        B=s.cancel(s.diff(U,t,2).subs(t,t0)/(2*U.subs(t,t0)))
        phi=[]
        for n in range(5):
            phi.append(sum(s.binomial(-h,j)*s.binomial(j,n-j)*A**(2*j-n)*B**(n-j)
                           for j in range((n+1)//2,n+1)))
        for nu in cfg['targets']:
            dot=sum(nu)-3
            factor=s.prod(d-3-j for j in range(dot))/s.prod(s.factorial(x-1) for x in nu)
            expected=factor*((1-z)*t)**(nu[0]-1)*((1-z)*(1-t))**(nu[1]-1)*z**(nu[2]-1)
            label='J{'+','.join(map(str,nu))+'}'
            self.assertFalse(model.targets[label][0]-model.rational(expected))
            rp=s.Poly(s.expand(expected.subs(t,t0+ta)),ta)
            for n in range(5):
                coeff=sum(rp.nth(j)*phi[n-j] for j in range(min(n,rp.degree())+1))
                for point in (2,3,7):
                    numerical=s.cancel(coeff.subs(z,point))
                    expected_mod=int(numerical.p)*pow(int(numerical.q),-1,model.p)%model.p
                    self.assertEqual(model.reducer.evaluate(result['raw_rows'][label][n][0],point),expected_mod)

    def test_false_train_relation_fails_holdout(self):
        key=('period',1,0,0)
        data={'a':[{key:1},{key:1},{key:1},{key:1}],
              'b':[{key:1},{key:1},{key:2},{key:1}]}
        found=search_relations(data,1000003,0,1)
        self.assertEqual(found['candidate_count'],1)
        self.assertEqual(found['rejected_count'],1)
        self.assertEqual(found['relations'],[])


if __name__=='__main__':
    unittest.main()
