"""Prove t separation when the two retained branch intrinsic polynomials vanish."""
import json
from pathlib import Path
import sympy as s
base=Path(__file__).resolve().parent;root=base.parents[3]
t,z,w,h,A,B,C=s.symbols('t z w h A B C');f,k=s.symbols('f k');v=A+B+C;a=t*(1-t);L=a+(1-a)*z;D=1-(1-a)*w;Z=a*w/D
U=(1-z)*L;F=z*(1-z)*(L*f-a*(1-z)*k)
assert s.cancel(Z/(a+(1-a)*Z)-w)==0
assert s.factor(s.diff(Z,w)-a/D**2)==0
assert s.cancel(F.subs(z,Z)-a**2*w*(1-w)/D**3*(f-(1-w)*k))==0
# Before transformation: (1-z)^(A+C-h-1) z^(2h-A-C-1)
# L^(v-3h) [L*f-a*(1-z)*k]^(2h-v) dz.
exponents=[A+C-h-1,2*h-A-C-1,v-3*h,2*h-v,2]
assert s.simplify(sum(exponents))==0
assert s.simplify((2*h-A-C-1)+(v-3*h)+(2*h-v)+1-(h-A-C))==0
# t^(A-1)(1-t)^(C-1)*a^(h-A-C).
assert s.simplify(A-1+h-A-C-(h-C-1))==0
assert s.simplify(C-1+h-A-C-(h-A-1))==0
# Old projection retains t inside the kernel.
Fold=(1-z)**2*(1-t)*(L*f-t*z*k)
assert s.cancel(Fold.subs(z,Z)-a*(1-t)*(1-w)**2/D**3*(f-t*w*k))==0
assert s.simplify((h-C-1)+(C-1)+(v-3*h)+(2*h-v)+2)==0
assert s.simplify(C-1+v-3*h+2*h-v+1-(C-h))==0
# Intrinsic branch polynomials are zero iff all masses vanish and all
# pairwise shift differences inside the retained branch are null.
fam=json.load(open(root/'j_integral/tests/e6/input/family_request.json'));G=s.Matrix([[s.Rational(x) for x in row] for row in fam['family']['gram']]);props=fam['family']['propagators'];checks=[]
for group in [[0,1],[5,6]]:
 for i in group:assert s.Rational(props[i]['mass2'])==0
 for i in group:
  for j in group:
   q=s.Matrix(props[i]['shift'])-s.Matrix(props[j]['shift']);assert (q.T*G*q)[0]==0
 checks.append([i+1 for i in group])
# Explicit duplicate in the original Blade list: both are the same massless sunrise.
scales=[]
for indices in [[0,2,5],[1,4,6]]:
    for i in indices:assert s.Rational(props[i]['mass2'])==0
    kvec=s.Matrix(props[indices[2]]['shift'])-s.Matrix(props[indices[0]]['shift'])-s.Matrix(props[indices[1]]['shift'])
    scales.append(s.factor((kvec.T*G*kvec)[0]))
assert scales==[s.Rational(1),s.Rational(1)]
fi=json.load(open(root/'j_integral/archive/j_integral_bootstrap/reference/fi_master_counts/e6.json'))
expected_basis={'J{'+','.join(map(str,nu[:7]))+'}' for nu in fi['master_indices'] if nu[:7]!=[1,0,1,0,0,1,0]}
assert len(fi['master_indices'])==62 and len(expected_basis)==61
# Every saved new-projection diagonal DE must equal the exact t power.
rows=[]
for name in ['j_result.json','j_dot2.json','j_blade_basis.json','j_blade_basis_second.json','j_blade_basis_t23.json']:
 path=base/'e6_full'/name
 if not path.exists():continue
 native=json.load(open(path));
 if name.startswith('j_blade_basis'):
  assert set(native['basis'])==expected_basis
  assert native['basis_status']=='provided_basis_independent_and_spans_requested_targets'
 p=native['point']['prime'];tt=s.Rational(native['point']['t']);hh=s.Rational(native['point']['dimension'])/2
 for label,deriv in native['derivatives'].items():
  nu=list(map(int,label.split('{')[1][:-1].split(',')));na=sum(nu[:2]);nc=sum(nu[5:]);value=s.cancel((hh-nc-1)/tt-(hh-na-1)/(1-tt));num,den=s.fraction(value);expected=int(num)%p*pow(int(den)%p,-1,p)%p
  assert deriv==({label:expected} if expected else {}),(name,label,deriv,expected)
 rows.append({'output':name,'master_count':len(native['basis']),'every_diagonal_exponent_matches':True,'DE_closed':not native['unclosed_derivatives']})
report={'generic_parameter_transformation_checked':True,'retained_zero_intrinsic_groups':checks,'Blade_duplicate_sunrise_invariants':['1','1'],'Blade_raw_count':62,'Blade_after_explicit_duplicate_merge':61,'t_factor':'t^(h-C-1)*(1-t)^(h-A-1)','normalized_identity':'J_nu(t)=t^(h-C-1)*(1-t)^(h-A-1)/Beta(h-C,h-A) * I_nu','diagonal_checks':rows}
(base/'theorem_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
