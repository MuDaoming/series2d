#define DIRECT_J_NO_MAIN
#include "../native/direct_j.cpp"
#undef DIRECT_J_NO_MAIN
#include <cassert>
int main(){
 prime=1000003;nmod_init(&mod,prime);
 Span s;assert(s.insert(C{{"a",1},{"b",2}},0));assert(s.insert(C{{"b",1}},1));
 assert(s.solve(C{{"a",3},{"b",10}},2)==std::vector<U>({3,4}));
 std::vector<U>xs,ys;for(U x=2;x<43;++x){xs.push_back(x);ys.push_back(nmod_div(add(mul(x,x),7),add(mul(x,x),3),mod));}
 Interpolator interpolator(std::vector<U>(xs.begin(),xs.begin()+25));Fraction f;
 assert(interpolator.reconstruct(ys,xs,f));assert(f.n.deg()==2&&f.d.deg()==2);
 ys.back()=add(ys.back(),1);assert(!interpolator.reconstruct(ys,xs,f));
 // Compare dense ansatz and Hermite on a two-component connection, including a nonzero endpoint.
 P z;z.set(1,1);P q=(z+P(1))*(z+P(2));Mat B(2,Row(2));B[0][0]=scale(z+P(2),nmod_inv(3,mod));B[0][1]=z+P(1);B[1][1]=scale(z+P(1),nmod_div(2,7,mod));
 Row a={z,P(1)},n(2);for(int j=0;j<2;++j){n[j]=q*derivative(a[j])-derivative(q)*a[j];for(int i=0;i<2;++i)n[j]+=a[i]*B[i][j];}
 RationalRow target;for(auto&x:n)target.push_back({x,q*q});std::vector<RationalRow>rows={target};
 for(auto&method:{std::string("hermite"),std::string("ansatz")}){PointIntegrator p(q,B,{"regular","regular"},method,rows,0,6);C result=p.integrate(target);RationalRow at;for(auto&x:a)at.push_back({x,q});C right=p.endpoint(at,1),lower=p.endpoint(at,0);axpy(right,lower,prime-1);assert(result==right);}
 std::cout<<"basis solve, rational reconstruction/heldout rejection, matrix z-IBP and both nonzero endpoints: PASS\n";
}
