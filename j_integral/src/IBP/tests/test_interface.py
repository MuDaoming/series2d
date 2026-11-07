"""Exercise numerical J-IBP and its explicit input contract."""
import json,subprocess,tempfile,unittest,sys,gzip
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BINARY=ROOT/'build/j-ibp'
REFERENCE='--reference' in sys.argv
if REFERENCE:sys.argv.remove('--reference')

class Interface(unittest.TestCase):
    def run_request(self,request,ok=True):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'request.json';dest=Path(directory)/'result.json'
            source.write_text(json.dumps(request))
            result=subprocess.run([str(BINARY),'--input',str(source),'--output',str(dest)],text=True,capture_output=True)
            if ok:
                self.assertEqual(result.returncode,0,result.stderr)
                return json.loads(dest.read_text())
            self.assertNotEqual(result.returncode,0,result.stderr)
            self.assertFalse(dest.exists(),'an invalid request must not replace a result')
    def vacuum(self):return json.loads((ROOT/'examples/vacuum.json').read_text())
    def test_numeric_reduction_and_fixed_basis_higher_dots(self):
        req=self.vacuum();out=self.run_request(req)
        self.assertEqual(out['status'],'complete');self.assertFalse(out['unclosed_derivatives'])
        basis=out['basis'];req['request']['basis']=basis;req['request']['generate_targets']['max_dot']=3
        higher=self.run_request(req)
        self.assertEqual(higher['basis'],basis)
        self.assertEqual(higher['derivatives'],out['derivatives'])
        for name,row in out['reductions'].items():self.assertEqual(higher['reductions'][name],row)
        # An independent numerical t must also close in this fixed basis.
        req['point']['t']='23';other=self.run_request(req);self.assertEqual(other['basis'],basis)
    def test_exact_symbolic_assignments_and_inferred_branches(self):
        req=self.vacuum();base=self.run_request(req)
        del req['j']['branches']
        for prop in req['family']['propagators']:prop['mass2']='1'
        req['point']['parameters']={}
        inferred=self.run_request(req)
        self.assertEqual(base['reductions'],inferred['reductions'])
        self.assertEqual(base['derivatives'],inferred['derivatives'])
    def test_double_box_is_only_a_different_input(self):
        req=json.loads((ROOT/'examples/double_box.json').read_text())
        out=self.run_request(req);self.assertEqual(out['status'],'complete')
        self.assertTrue(out['basis']);self.assertFalse(out['unclosed_derivatives'])
    def test_input_errors_are_rejected(self):
        changes=[lambda r:r['point'].update(prime=1000000),
                 lambda r:r['family']['propagators'][0].update(mass2=0.1),
                 lambda r:r['point'].update(parameters={}),
                 lambda r:r['request'].update(targets=[[1,1,1]]),
                 lambda r:r['request'].update(basis=['J{1,1,1}','J{1,1,1}']),
                 lambda r:r['j'].update(cut_sector=[1,1,1])]
        for change in changes:
            with self.subTest(change=change):
                req=self.vacuum();change(req);self.run_request(req,False)
        req=self.vacuum();del req['request']['generate_targets'];req['request']['targets']=[[-1,1,1]];self.run_request(req,False)
        req=self.vacuum();del req['request']['generate_targets'];req['request']['targets']=[[1,0,1]];self.run_request(req,False)
    def test_r_only_and_requested_derivative(self):
        req=self.vacuum();req['request']['derivatives']=[]
        out=self.run_request(req);self.assertEqual(out['derivatives'],{})
        req['request']['basis']=out['basis'];req['request']['derivatives']=[out['basis'][0]]
        one=self.run_request(req);self.assertEqual(set(one['derivatives']),{out['basis'][0]})
    @unittest.skipUnless(REFERENCE,'use --reference to run the DP independent-reference check')
    def test_dp_full_connection_against_frozen_rational_reference(self):
        req=json.loads((ROOT/'examples/dp.json').read_text());req['request']['basis']=json.loads((ROOT/'examples/dp_basis.json').read_text());req['point']['t']='23'
        out=self.run_request(req)
        archive=ROOT.parents[1]/'archive'
        with gzip.open(archive/'j_integral_bootstrap/reference/dp_global/connection_mod1000003.json.gz','rt') as file:ref=json.load(file)
        self.assertEqual(out['basis'],ref['masters']);p=ref['prime'];center=Fraction(ref['t0']);x=(23-center.numerator*pow(center.denominator,-1,p))%p
        def ev(poly):
            result=0
            for c in reversed(poly):result=(result*x+c)%p
            return result
        count=0
        for row in ref['masters']:
            actual=out['derivatives'][row]
            for col in ref['masters']:
                rf=ref['de'][row]['rhs'].get(col)
                expected=ev(rf['numerator_tau'])*pow(ev(rf['denominator_tau']),-1,p)%p if rf else 0
                self.assertEqual(actual.get(col,0),expected,(row,col));count+=1
        self.assertEqual(count,496*496)
        print('DP: all246016 entries including zeros agree with the independently reconstructed rational connection.')

if __name__=='__main__':unittest.main()
