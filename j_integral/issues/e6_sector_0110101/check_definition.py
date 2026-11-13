"""Verify the independent four-propagator input and its inherited numerical DE."""
import json,sys
from pathlib import Path
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
load=lambda p:json.loads(Path(p).read_text())
request=load(HERE/'input/family_request.json');mapping=load(HERE/'input/index_map.json');original=load(ROOT/'j_integral/tests/e6/input/family_request.json');G=[[Q(v) for v in row] for row in original['family']['gram']]
P=[1,0,-1,-1,0];R=[1,1,-1,-1,0]
def dot(a,b):return sum(a[i]*G[i][j]*b[j] for i in range(5) for j in range(5))
assert Q(request['point']['parameters']['sP'])==dot(P,P)
assert Q(request['point']['parameters']['sQ'])==dot(R,R)
assert Q(mapping['invariants']['P_dot_Q'])==dot(P,R)
diff=[b-a for a,b in zip(P,R)];assert Q(request['point']['parameters']['sR'])==dot(diff,diff)
# Undo ell1=k1+p3 in each new propagator and check the original linear momenta.
for prop,old_index in zip(request['family']['propagators'],mapping['original_propagator_numbers']):
    old=original['family']['propagators'][old_index-1];a,b=map(Q,prop['loop']);u,v=map(Q,prop['shift']);shift=[u*P[i]+v*R[i]+(a if i==2 else 0) for i in range(5)]
    assert [a,b]==list(map(Q,old['loop'])) and shift==list(map(Q,old['shift']))
    assert Q(prop['mass2'])==Q(old['mass2'])
result=dict(momentum_routing_and_kinematics_equal=True)
if len(sys.argv)>1:
    actual=load(sys.argv[1]);reference=load(HERE/'reference/j_de.json');p=reference['prime'];c=Q(reference['t0']);center=c.numerator*pow(c.denominator,-1,p)%p;t=int(actual['point']['t']);x=(t-center)%p
    assert actual['basis']==reference['masters'] and actual['status']=='complete' and not actual['unclosed_derivatives']
    def ev(poly):
        total=0
        for v in reversed(poly):total=(total*x+v)%p
        return total
    count=0
    for row in reference['masters']:
        for col in reference['masters']:
            value=reference['de'][row]['rhs'].get(col);expected=ev(value['numerator_tau'])*pow(ev(value['denominator_tau']),-1,p)%p if value else 0
            assert actual['derivatives'][row].get(col,0)==expected,(row,col);count+=1
    result.update(point=t,DE_scalar_checks_including_zeros=count,all_DE_coefficients_equal=True)
(HERE/'reference/definition_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
