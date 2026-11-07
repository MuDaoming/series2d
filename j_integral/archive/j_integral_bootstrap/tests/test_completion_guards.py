"""Reduction completion uses target coverage, independently of extra masters."""
import argparse
import json
from pathlib import Path
import tempfile
import unittest
from bootstrap_reductions import main


class CompletionGuards(unittest.TestCase):
    def run_case(self, root, with_target):
        m0,m1,target='J{1,0}','J{0,1}','J{1,1}'
        (root/'points').mkdir()
        connection=dict(prime=101,t0='0',masters=[m0,m1],known_rows=[m0],status='running',
                        de={m0:dict(rhs={}),m1:dict(rhs={})})
        point=dict(prime=101,t_point=2,masters_at_point=[m0,m1],unclosed_derivatives=[],
                   DE_values={},reduction_values={m0:{m0:1},target:{m1:1}})
        (root/'connection.json').write_text(json.dumps(connection))
        (root/'points/point2.json').write_text(json.dumps(point))
        partial=[]
        if with_target:
            (root/'partial.json').write_text(json.dumps(dict(prime=101,t0='0',
                reductions={target:{m1:dict(numerator_tau=[1],denominator_tau=[1])}})))
            partial=[str(root/'partial.json')]
        code=main(argparse.Namespace(work=str(root/'work'),connection=str(root/'connection.json'),
            points=[str(root/'points')],partial=partial,witness_budget=0,internal_holdouts=3))
        return code,json.loads((root/'work/reduction.json').read_text()),target

    def test_equal_cardinality_does_not_hide_a_missing_target(self):
        with tempfile.TemporaryDirectory() as td:
            code,result,target=self.run_case(Path(td),False)
            self.assertEqual(code,2)
            self.assertEqual(result['status'],'waiting_for_connection')
            self.assertIn(target,result['pending_targets'])
            self.assertNotIn(target,result['reductions'])

    def test_extra_masters_do_not_prevent_complete_target_coverage(self):
        with tempfile.TemporaryDirectory() as td:
            code,result,target=self.run_case(Path(td),True)
            self.assertEqual(code,0)
            self.assertEqual(result['status'],'complete')
            self.assertTrue(set(result['tested_objects'])<=set(result['known_targets']))
            self.assertIn(target,result['reductions'])


if __name__=='__main__':unittest.main()
