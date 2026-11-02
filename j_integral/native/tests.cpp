#define main j_series_entry_for_tests
#include "series.cpp"
#undef main
#include <cassert>
int main(){
 prime=1000003;nmod_init(&mod,prime);
 // q-adic conversion must be exact, including the Newton-division path.
 P q;q.set(0,3);q.set(1,2);q.set(3,1);Qadic qa(q);
 for(int d:{0,1,2,3,9,57,200,600}){P p;for(int j=0;j<=d;++j)p.set(j,(j*97+113)%prime);auto v=qa.digits(p);assert((qa.compose(v,0,v.size())-p).zero());}
 // Matrix orientation and nonzero regular endpoint; input is exactly nabla(g).
 P z;z.set(1,1);P q2=z*(z+P(1));Mat B(2,Row(2));
 B[0][0]=scale(z+P(1),nmod_inv(3,mod));B[0][1]=scale(z,nmod_inv(5,mod));B[1][1]=scale(z+P(1),nmod_div(2,7,mod));
 Row G={z,P(1)},input(2);for(int j=0;j<2;++j){input[j]=q2*derivative(G[j])-derivative(q2)*G[j];for(int i=0;i<2;++i)input[j]+=G[i]*B[i][j];}
 Reducer reducer(q2,B,{"regulated_zero","regular"});auto c=reducer.reduce(input,2,true);assert(c.size()==2);assert(c.at("endpoint:1:0")==nmod_inv(2,mod));assert(c.at("endpoint:1:1")==nmod_inv(2,mod));
 // A resonant pole-map obstruction is a period, not zero.
 Reducer resonant(z,Mat{{P(1)}},{"regulated_zero","regular"});auto obs=resonant.reduce(Row{P(1)},2,true);assert(obs.at("period:2:0:0")==1);
 // A null direction at infinity supplies an additional exact period relation.
 P qr=(z+P(1))*(z+P(2))*(z+P(3));P br=scale((z+P(2))*(z+P(3)),nmod_div(prime-1,3,mod))-scale((z+P(1))*(z+P(3)),nmod_div(2,3,mod))-(z+P(1))*(z+P(2));
 Reducer completed(qr,Mat{{br}},{"regular","regular"});Row gg={z*z};Row dg={qr*derivative(gg[0])+gg[0]*br};auto before=completed.reduce(dg,1,true);assert(before.size()>1);completed.complete_relations(8);assert(!completed.relations.empty());auto after=completed.reduce(dg,1,true);assert(after.size()==1&&after.at("endpoint:1:0")==1);
 // Independent all-order vac certificate for the two-function closed basis.
 nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);
 RF t=Expression("t").parse(),zz=Expression("z").parse(),h=RF(7)/RF(26),a=t*(RF(1)-t),L=a+(RF(1)-a)*zz;
 RF At=-h*(RF(1)-zz)*a.diff(0)/L;
 RF Az=-(RF(1)-h)/(RF(1)-zz)-h*(RF(1)-a)/L;
 RF g=a.diff(0)*(RF(1)-zz)/(RF(1)-a);
 assert((At-a.diff(0)*(RF(2)-h)/(RF(1)-a)-g.diff(1)-g*Az).zero());
 std::cout<<"q-adic conversion, matrix direction/endpoints, resonant obstruction, exact vac DE certificate: PASS\n";
}
