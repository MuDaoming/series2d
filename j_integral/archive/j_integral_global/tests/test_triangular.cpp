#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"
#include "../native/triangular.hpp"
#include <cassert>
int main(){prime=1000003;nmod_init(&mod,prime);P z;z.set(1,1);ZR one(1),minus_one(prime-1);
 // A Laurent-log fundamental solution: f1=(c1+c2 log z)/z, f2=c2.
 std::vector<std::vector<ZR>>A(2,std::vector<ZR>(2));A[0][0]=ZR(P(prime-1),z);A[0][1]=ZR(P(1),z*z);EndpointGerm germ(A,6);
 assert(germ.finite_part(0,zp(z))==C({{"endpoint:0:0",1}}));assert(germ.finite_part(0,one).empty());assert(germ.finite_part(1,one)==C({{"endpoint:0:1",1}}));
 A[0][0]=ZR(P(prime-2),z);A[0][1]=ZR(P(1),z*z*z);EndpointGerm order_two(A,6);assert(order_two.finite_part(0,zp(z*z))==C({{"endpoint:0:0",1}}));A[0][0]=ZR(P(prime-1),z);A[0][1]=ZR(P(1),z*z);
 // Local Hermite identity, resonant pole and ordinary pole in the same input.
 Factors factors;ScalarHermite reducer(ZR(P(1),z),factors);ZR target=ZR(P(1),z*z)+ZR(P(7),(z+P(2))*(z+P(2))*(z+P(2)));auto [h,g]=reducer.reduce(target);assert((target-h-g.diff()-g*reducer.A).zero());
 // An off-diagonal total derivative must propagate into lower components.
 TriangularIntegrator tri(A,8);tri.complete(8);std::vector<ZR>v(2);ZR gp=zp(z);v[0]=gp.diff()+gp*A[0][0];v[1]=gp*A[0][1];C c=tri.integrate(v);axpy(c,germ.finite_part(0,gp),1);assert(tri.completion.reduce(c).empty());
 std::cout<<"Laurent-log endpoint, resonant local reduction, coupled triangular total derivative: PASS\n";
}
