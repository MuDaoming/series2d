"""Rational tau coefficients cannot be reused with a different expansion center."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest
import bootstrap_global as bg
import bootstrap_reductions as br


class CoordinateGuards(unittest.TestCase):
    p=101
    master='J{1}'

    def write(self,path,data):
        path.write_text(json.dumps(data))
        return str(path)

    def global_case(self,root,diagonal_center=None,checkpoint_center=None,partial_center=None,
                    omit_diagonal_center=False,omit_checkpoint_center=False):
        m=self.master;work=root/'work';work.mkdir();points=root/'points';points.mkdir()
        field_center=61 # 2/5 mod 101
        for t in range(2,11):
            self.write(points/f'point{t}.json',dict(prime=self.p,t_point=t,masters_at_point=[m],
                unclosed_derivatives=[],DE_values={m:{m:pow((t-field_center)%self.p,-1,self.p)}}))
        coefficient=dict(numerator_tau=[1],denominator_tau=[0,1])
        row={m:dict(rhs={m:coefficient})}
        diagonals=[];partial=[]
        if diagonal_center is not None or omit_diagonal_center:
            value=dict(prime=self.p,t0=diagonal_center,de=row)
            if omit_diagonal_center:del value['t0']
            diagonals=[self.write(root/'diagonal.json',value)]
        if partial_center is not None:
            partial=[self.write(root/'partial.json',dict(prime=self.p,t0=partial_center,de={}))]
        checkpoint=work/'connection.json'
        if checkpoint_center is not None or omit_checkpoint_center:
            value=dict(prime=self.p,t0=checkpoint_center,masters=[m],known_rows=[m],blocks=[],de=row)
            if omit_checkpoint_center:del value['t0']
            self.write(checkpoint,value)
        args=argparse.Namespace(work=str(work),family=self.write(root/'family.json',dict(prime=self.p,t0='2/5')),
            basis=self.write(root/'basis.json',[m]),points=[str(points)],diagonals=diagonals,partial=partial,
            restart=False,witness_budget=0,internal_holdouts=3,backend='PL-structural',relax_all_poles=False,
            extra_order=0)
        return args,checkpoint

    def test_global_rejects_wrong_or_missing_diagonal_center(self):
        for options in [dict(diagonal_center='0'),dict(omit_diagonal_center=True)]:
            with self.subTest(options=options),tempfile.TemporaryDirectory() as td:
                args,checkpoint=self.global_case(Path(td),**options)
                with self.assertRaisesRegex(ValueError,'diagonal field/coordinate'):
                    bg.main(args)
                self.assertFalse(checkpoint.exists())

    def test_global_rejects_wrong_or_missing_checkpoint_center_without_rewriting(self):
        for options in [dict(checkpoint_center='0'),dict(omit_checkpoint_center=True)]:
            with self.subTest(options=options),tempfile.TemporaryDirectory() as td:
                args,checkpoint=self.global_case(Path(td),**options);before=checkpoint.read_bytes()
                with self.assertRaisesRegex(ValueError,'checkpoint field/basis/coordinate'):
                    bg.main(args)
                self.assertEqual(checkpoint.read_bytes(),before)

    def test_global_accepts_equivalent_rational_centers(self):
        for options in [dict(diagonal_center='4/10'),dict(checkpoint_center='4/10')]:
            with self.subTest(options=options),tempfile.TemporaryDirectory() as td:
                args,checkpoint=self.global_case(Path(td),**options)
                self.assertEqual(bg.main(args),0)
                result=json.loads(checkpoint.read_text())
                self.assertEqual(result['status'],'complete')
                self.assertEqual(Fraction(result['t0']),Fraction(2,5))

    def reduction_case(self,root,checkpoint_center):
        m=self.master;work=root/'work';work.mkdir();points=root/'points';points.mkdir()
        connection=dict(prime=self.p,t0='2/5',masters=[m],known_rows=[m],status='complete',de={m:dict(rhs={})})
        self.write(points/'point2.json',dict(prime=self.p,t_point=2,masters_at_point=[m],
            unclosed_derivatives=[],DE_values={m:{}},reduction_values={m:{m:1}}))
        checkpoint=work/'reduction.json'
        saved=dict(prime=self.p,t0=checkpoint_center,masters=[m],known_targets=[m],groups=[],reductions={})
        if checkpoint_center is None:del saved['t0']
        self.write(checkpoint,saved)
        return argparse.Namespace(work=str(work),connection=self.write(root/'connection.json',connection),
            points=[str(points)],partial=[],witness_budget=0,internal_holdouts=3),checkpoint

    def test_reduction_rejects_wrong_or_missing_checkpoint_center_without_rewriting(self):
        for center in ['0',None]:
            with self.subTest(center=center),tempfile.TemporaryDirectory() as td:
                args,checkpoint=self.reduction_case(Path(td),center);before=checkpoint.read_bytes()
                with self.assertRaisesRegex(ValueError,'checkpoint field/basis/coordinate'):
                    br.main(args)
                self.assertEqual(checkpoint.read_bytes(),before)

    def test_reduction_accepts_equivalent_rational_center(self):
        with tempfile.TemporaryDirectory() as td:
            args,checkpoint=self.reduction_case(Path(td),'4/10')
            self.assertEqual(br.main(args),0)
            self.assertEqual(json.loads(checkpoint.read_text())['status'],'complete')

    def test_existing_partial_coordinate_guards_remain_active(self):
        with tempfile.TemporaryDirectory() as td:
            args,checkpoint=self.global_case(Path(td),partial_center='0')
            with self.assertRaisesRegex(ValueError,'partial reconstruction field/coordinate'):
                bg.main(args)
            self.assertFalse(checkpoint.exists())
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);args,checkpoint=self.reduction_case(root,'2/5')
            args.partial=[self.write(root/'partial.json',dict(prime=self.p,t0='0',reductions={}))]
            before=checkpoint.read_bytes()
            with self.assertRaisesRegex(ValueError,'partial field/coordinate'):
                br.main(args)
            self.assertEqual(checkpoint.read_bytes(),before)


if __name__=='__main__':unittest.main()
