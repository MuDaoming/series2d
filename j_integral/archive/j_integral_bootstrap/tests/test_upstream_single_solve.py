"""Exercise real upstream learning/production: partial solves add prefix checks."""
from collections import Counter
import unittest
from unittest.mock import patch
from reference_space import RatMat
from analytic_bootstrap.expressions import connection_matrix
from analytic_bootstrap.upstream import learn_anchor,prepare_setting
from analytic_bootstrap.upstream.model import matrix_record,submatrix
from analytic_bootstrap.storage import digest

FAMILY=dict(family={'name':'single-solve-audit'},variable='x',sectors=[
    dict(id='base',masters=['I0'],depends_on=[]),
    dict(id='top',masters=['I1'],depends_on=['base'])],reductions={})
SETTING=dict(prime=1000003,assignments={'a':'1'})


def truth(req):
    return connection_matrix([['0','0'],['2*x/(x^2+a)','2*x/(x^2+a)']],
                             'x',req['assignments'],req['prime'])


class Anchor:
    def __init__(self):self.calls=0
    def anchor(self,req):
        self.calls+=1
        return dict(request_sha256=digest(req),de=matrix_record(truth(req)),reductions={})


class Cuts:
    def diagonal(self,req):
        i=0 if req['sector']=='base' else 1
        return dict(request_sha256=digest(req),matrix=matrix_record(submatrix(truth(req),[i],[i])))


class Numeric:
    """A pointwise exact toy oracle, caching one full evaluation per x."""
    def __init__(self):self.events=[];self.cache={}
    def prepare(self,req):
        owner=self
        class Session:
            def evaluate(self,x):
                owner.events.append(dict(task=req['task'],x=x))
                if x not in owner.cache:
                    matrix=truth(req)
                    owner.cache[x]=[[matrix.entries[i,j].evaluate(x) if (i,j) in matrix.entries else 0
                                     for j in range(2)] for i in range(2)]
                ids=['I0','I1'];columns=[ids.index(m) for m in req['masters']]
                values=[[owner.cache[x][ids.index(t['master'])][j] for j in columns]
                        for t in req['targets']]
                return dict(request_sha256=digest(req),x=x,values=values)
            def close(self):pass
        return Session()


def solve_toy(bundle,sector_loop):
    numeric=Numeric();prefix=[]
    with prepare_setting(bundle,SETTING,Cuts(),numeric,shared_de_samples=True,holdouts=2,
                         global_holdouts=2,fit_backend='python',
                         builder_options={'laurent_mode':'classic'}) as point:
        with patch.object(point,'_check_complete_de',wraps=point._check_complete_de) as check:
            if sector_loop:
                for limit in range(1,len(FAMILY['sectors'])+1):
                    part=point.solve_de(_internal=True,_sector_limit=limit)
                    if limit<len(FAMILY['sectors']):prefix.append(part.global_check)
            else:point.solve_de(_internal=True)
            solved=point.solve_de(_internal=True)
            calls=check.call_count;before=len(numeric.events)
            point.solve_de();point.solve_de(_internal=True)
            if len(numeric.events)!=before or check.call_count!=calls:
                raise AssertionError('cached completed solves performed numerical checks')
        exact=matrix_record(solved.matrix)==matrix_record(truth(SETTING))
        return dict(exact=exact,checks=solved.checks,global_check=solved.global_check,
                    prefix_checks=prefix,complete_check_calls=calls,
                    numeric_unique_points=len(numeric.cache),
                    numeric_tasks=dict(Counter(e['task'] for e in numeric.events)))


class UpstreamSingleSolve(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.anchor=Anchor()
        # This is the actual package learning path, including source profiles.
        cls.bundle=learn_anchor(FAMILY,SETTING,cls.anchor,learn_profiles=True)
        cls.single=solve_toy(cls.bundle,False)
        cls.loop=solve_toy(cls.bundle,True)

    def test_one_solve_performs_only_final_global_validation(self):
        result=self.single
        self.assertEqual(self.anchor.calls,1)
        self.assertTrue(result['exact'])
        self.assertEqual(result['complete_check_calls'],1)
        self.assertFalse(result['prefix_checks'])
        self.assertFalse(any(k.startswith('__prefix_de_holdouts__') for k in result['numeric_tasks']))
        self.assertEqual(result['numeric_tasks']['__complete_de_holdouts__'],2)
        self.assertEqual(result['global_check']['status'],'fresh-full-de-holdouts')

    def test_per_sector_calls_add_prefix_validation(self):
        result=self.loop
        self.assertTrue(result['exact'])
        self.assertEqual(result['complete_check_calls'],2)
        self.assertEqual(result['numeric_tasks']['__prefix_de_holdouts__1'],2)
        self.assertEqual(result['numeric_tasks']['__complete_de_holdouts__'],2)
        self.assertEqual(result['prefix_checks'][0]['status'],'fresh-prefix-de-holdouts')
        self.assertEqual(result['numeric_unique_points'],self.single['numeric_unique_points']+2)
        self.assertEqual(result['checks'],self.single['checks'])

    def test_removing_prefix_calls_changes_final_fresh_points(self):
        old=self.loop;new=self.single
        self.assertEqual(new['global_check']['holdout_points'],old['prefix_checks'][0]['holdout_points'])
        self.assertNotEqual(new['global_check']['holdout_points'],old['global_check']['holdout_points'])
        used={x for c in new['checks'] for key in ['training_points','holdout_points'] for x in c[key]}
        self.assertFalse(used&set(new['global_check']['holdout_points']))


if __name__=='__main__':unittest.main()
