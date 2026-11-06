#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"
#include "../native/triangular.hpp"
#include <cassert>
ZR rz(const RF&x){return ZR(univariate(x.n),univariate(x.d));}
int main(int argc,char**argv){assert(argc==2);json cfg;{std::ifstream f(argv[1]);f>>cfg;}prime=cfg["prime"];nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);cfg["specialize_t"]="17";NativeFamily family(cfg);std::vector<Nu>targets;for(auto&[s,info]:family.sectors){Nu nu=s;std::function<void(int,int)>visit=[&](int i,int k){if(i==family.n){targets.push_back(nu);return;}if(!s[i]){visit(i+1,k);return;}for(int j=0;j<=k;++j){nu[i]=1+j;visit(i+1,k-j);}nu[i]=1;};visit(0,2);}family.cfg["targets"]=targets;RMat At,Az;std::map<std::string,RRow>raw;family.model(At,Az,raw);int r=family.rank;
 std::vector<std::vector<ZR>>A(r);std::vector<RationalRow>all;std::map<std::string,RationalRow>rows,derivatives;P q(1);auto take=[&](const RF&x){auto v=specialize(x,17);P g;nmod_poly_gcd(g.v,q.v,v.second.v);q=quotient(q,g)*v.second;return v;};
 for(int i=0;i<r;++i)for(auto&x:Az[i]){A[i].push_back(rz(x));take(x);}for(auto&[name,row]:raw){for(auto&x:row)rows[name].push_back(take(x));all.push_back(rows[name]);}for(auto&nu:targets){std::string name="J{";for(int i=0;i<family.n;++i)name+=(i?",":"")+std::to_string(nu[i]);name+="}";for(auto&x:family.target_derivative(nu,0))derivatives[name].push_back(take(x));all.push_back(derivatives[name]);}
 Mat B(r,Row(r));for(int i=0;i<r;++i)for(int j=0;j<r;++j){auto x=specialize(Az[i][j],17);B[i][j]=quotient(q*x.first,x.second);}PointIntegrator dense(q,B,{"regular","regulated_zero"},"ansatz",all,1,16);TriangularIntegrator tri(A,16);tri.complete(16);
 auto convert=[](const RationalRow&row){std::vector<ZR>v;for(auto&x:row)v.emplace_back(x.first,x.second);return v;};
 std::map<std::string,C>tc,dc;for(auto&[name,row]:rows){tc[name]=tri.integrate(convert(row));dc[name]=dense.integrate(row);}
 Span basis;std::vector<std::string>masters;for(auto&[name,c]:tc)if(basis.insert(c,masters.size()))masters.push_back(name);Span db;for(int i=0;i<int(masters.size());++i)assert(db.insert(dc.at(masters[i]),i));
 int checks=0;for(auto&[name,c]:tc){assert(basis.solve(c,masters.size())==db.solve(dc.at(name),masters.size()));++checks;}for(auto&name:masters){C t=tri.integrate(convert(derivatives.at(name))),d=dense.integrate(derivatives.at(name));assert(basis.solve(t,masters.size())==db.solve(d,masters.size()));++checks;}
 std::cout<<"Independent dense coupled IBP vs factorwise triangular IBP: r="<<r<<" targets="<<targets.size()<<" basis="<<masters.size()<<" checks="<<checks<<" PASS\n";
}
