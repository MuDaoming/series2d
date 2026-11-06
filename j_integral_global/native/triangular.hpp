#pragma once
#include <flint/nmod_poly_factor.h>
// Univariate rational arithmetic after t specialization.
struct ZR {
 P n,d;
 ZR():d(1){} explicit ZR(U c):n(c),d(1){} ZR(P n,P d):n(std::move(n)),d(std::move(d)){normal();}
 bool zero()const{return n.zero();}
 void normal(){if(d.zero())throw std::runtime_error("zero z denominator");if(n.zero()){d=P(1);return;}P g;nmod_poly_gcd(g.v,n.v,d.v);n=quotient(n,g);d=quotient(d,g);U c=nmod_inv(d.at(d.deg()),mod);n=scale(n,c);d=scale(d,c);}
 ZR diff()const{return ZR(derivative(n)*d-n*derivative(d),d*d);}
};
ZR operator+(const ZR&a,const ZR&b){if(a.zero())return b;if(b.zero())return a;P g;nmod_poly_gcd(g.v,a.d.v,b.d.v);P ad=quotient(a.d,g),bd=quotient(b.d,g);return ZR(a.n*bd+b.n*ad,ad*b.d);}
ZR operator-(const ZR&a,const ZR&b){return a+ZR(scale(b.n,prime-1),b.d);}
ZR operator*(const ZR&a,const ZR&b){if(a.zero()||b.zero())return ZR();return ZR(a.n*b.n,a.d*b.d);}
ZR zscale(const ZR&a,U c){return ZR(scale(a.n,c),a.d);}
ZR zp(const P&a){return ZR(a,P(1));}
std::string pkey(const P&p){return polyjson(p).dump();}
int valuation(const P&p){if(p.zero())return 1000000;int n=0;while(!p.at(n))++n;return n;}
int valuation(const ZR&p){return p.zero()?1000000:valuation(p.n)-valuation(p.d);}
struct Factors {
 std::map<std::string,std::vector<std::pair<P,int>>>cache;
 const std::vector<std::pair<P,int>>&get(const P&den){std::string key=pkey(den);auto it=cache.find(key);if(it!=cache.end())return it->second;std::vector<std::pair<P,int>>v;nmod_poly_factor_t f;nmod_poly_factor_init(f);nmod_poly_factor(f,den.v);for(slong i=0;i<f->num;++i){P p;nmod_poly_set(p.v,f->p+i);nmod_poly_make_monic(p.v,p.v);v.push_back({std::move(p),int(f->exp[i])});}nmod_poly_factor_clear(f);return cache.emplace(key,std::move(v)).first->second;}
};
P modp(const P&a,const P&f){P p;nmod_poly_rem(p.v,a.v,f.v);return p;}
P invert_mod(const P&a,const P&f){P p;if(!nmod_poly_invmod(p.v,a.v,f.v))throw std::runtime_error("noninvertible local residue");return p;}
struct Laurent {
 int low=0;std::vector<U>a;
 U at(int k)const{return k<low||k>=low+int(a.size())?0:a[k-low];}
};
Laurent laurent(const ZR&f,int high){Laurent out;if(f.zero())return out;int vn=valuation(f.n),vd=valuation(f.d);out.low=vn-vd;if(high<out.low)return out;int count=high-out.low+1;P n,d,di,p;nmod_poly_shift_right(n.v,f.n.v,vn);nmod_poly_shift_right(d.v,f.d.v,vd);nmod_poly_inv_series(di.v,d.v,count);nmod_poly_mullow(p.v,n.v,di.v,count);for(int i=0;i<count;++i)out.a.push_back(p.at(i));return out;}
using LogJet=std::vector<C>;
void log_axpy(LogJet&v,const LogJet&w,U c){v.resize(std::max(v.size(),w.size()));for(size_t i=0;i<w.size();++i)axpy(v[i],w[i],c);while(!v.empty()&&v.back().empty())v.pop_back();}
// Regularized z=0 boundary functional: coefficient of z^0 log(z)^0.
// Integer-power Laurent-log sector of the triangular FBI connection.
struct EndpointGerm {
 int r;std::vector<int>low,high;std::vector<std::map<int,LogJet>>jets;
 EndpointGerm(const std::vector<std::vector<ZR>>&A,int order):r(A.size()),low(r,0),high(r,order),jets(r){
  for(int i=0;i<r;++i)for(int j=i+1;j<r;++j)if(!A[i][j].zero())high[j]=std::max(high[j],high[i]-1-valuation(A[i][j]));
  for(int i=r-1;i>=0;--i){if(valuation(A[i][i])<-1)throw std::runtime_error("irregular diagonal at z=0");U lambda=laurent(A[i][i],-1).at(-1);int integer_bound=order+r+16;low[i]=lambda>=prime-U(integer_bound)?int(lambda)-int(prime):0;for(int j=i+1;j<r;++j)if(!A[i][j].zero())low[i]=std::min(low[i],low[j]+1+valuation(A[i][j]));
   std::vector<Laurent>a(r);for(int j=i;j<r;++j)if(!A[i][j].zero())a[j]=laurent(A[i][j],high[i]-1-low[j]);
   for(int n=low[i];n<=high[i];++n){LogJet rhs;for(int j=i;j<r;++j){if(A[i][j].zero())continue;for(int k=a[j].low;k<int(a[j].low+a[j].a.size());++k){if(i==j&&k<0)continue;U coeff=a[j].at(k);int index=n-1-k;if(!coeff)continue;auto it=jets[j].find(index);if(it!=jets[j].end())log_axpy(rhs,it->second,coeff);}}
    U factor=sub(n<0?(prime-U(-n)%prime)%prime:U(n)%prime,lambda);LogJet value;
    if(!factor){value.resize(rhs.size()+1);value[0]={{"endpoint:0:"+std::to_string(i),1}};for(size_t l=0;l<rhs.size();++l)value[l+1]=scaled(rhs[l],nmod_inv(l+1,mod));}
    else{value=rhs;U inv=nmod_inv(factor,mod);for(int l=int(value.size())-1;l>=0;--l){if(l+1<int(value.size()))axpy(value[l],value[l+1],prime-U(l+1));value[l]=scaled(value[l],inv);}}
    while(!value.empty()&&value.back().empty())value.pop_back();if(!value.empty())jets[i].emplace(n,std::move(value));
   }
  }
 }
 C finite_part(int i,const ZR&g)const{C out;if(g.zero())return out;int v=valuation(g);if(-v>high[i])throw std::runtime_error("endpoint germ order insufficient");auto a=laurent(g,-low[i]);for(auto&[n,v]:jets[i])if(!v.empty()){U coeff=a.at(-n);if(coeff)axpy(out,v[0],coeff);}return out;}
};
struct ScalarHermite {
 ZR A;Factors* factors;long long checks=0;
 ScalarHermite(ZR A,Factors&f):A(std::move(A)),factors(&f){}
 P residue(const P&f){P g;nmod_poly_gcd(g.v,A.d.v,f.v);if(g.deg()<1)return P();P rest=quotient(A.d,f);if(modp(rest,f).zero())throw std::runtime_error("non-Fuchsian diagonal");return modp(A.n*invert_mod(modp(rest,f),f),f);}
 std::pair<ZR,ZR>reduce(ZR v){ZR original=v,primitive,held;
  auto initial=factors->get(v.d);for(auto&[f,power]:initial){P df=derivative(f),res=residue(f);for(int k=power;k>=2;--k){P fk=ppow(f,k),rest;if(!nmod_poly_divides(rest.v,v.d.v,fk.v))continue;P part=modp(v.n*invert_mod(modp(rest,f),f),f);P map=modp(res-scale(df,k-1),f);if(map.zero()){ZR obstruction(part,fk);held=held+obstruction;v=v-obstruction;continue;}P a=modp(part*invert_mod(map,f),f);ZR g(a,ppow(f,k-1));v=v-g.diff()-g*A;primitive=primitive+g;}}
  if(A.n.deg()>=A.d.deg())throw std::runtime_error("polynomial diagonal at infinity unsupported");U lambda=A.n.deg()+1==A.d.deg()?nmod_div(A.n.at(A.n.deg()),A.d.at(A.d.deg()),mod):0;
  // Leave resonant monomials as cohomology representatives.
  P polynomial,rem;nmod_poly_divrem(polynomial.v,rem.v,v.n.v,v.d.v);for(int d=polynomial.deg();d>=0;--d){U leading=polynomial.at(d),den=add(U(d+1)%prime,lambda);if(!leading||!den)continue;P p;p.set(d+1,nmod_div(leading,den,mod));ZR g=zp(p);v=v-g.diff()-g*A;primitive=primitive+g;nmod_poly_divrem(polynomial.v,rem.v,v.n.v,v.d.v);}
  // Infinity also has the z^(-1) layer: a constant primitive removes it when lambda != 0.
  nmod_poly_divrem(polynomial.v,rem.v,v.n.v,v.d.v);
  if(lambda&&rem.deg()+1==v.d.deg()){
   U c=nmod_div(rem.at(rem.deg()),v.d.at(v.d.deg()),mod);ZR g(nmod_div(c,lambda,mod));v=v-g*A;primitive=primitive+g;
  }
  v=v+held;
  if(!(original-v-primitive.diff()-primitive*A).zero())throw std::runtime_error("scalar Hermite rational certificate failed");++checks;return {std::move(v),std::move(primitive)};
 }
 void coordinates(const ZR&v,int row,C&out){P poly,proper;nmod_poly_divrem(poly.v,proper.v,v.n.v,v.d.v);std::string prefix="a:"+std::to_string(row)+":";for(int j=0;j<=poly.deg();++j)if(poly.at(j))out[prefix+"poly:"+std::to_string(j)]=poly.at(j);
  for(auto&[f,k]:factors->get(v.d)){P fk=ppow(f,k),rest=quotient(v.d,fk);P digit=modp(proper*invert_mod(modp(rest,fk),fk),fk);for(int power=k;power>=1;--power){P hi,lo;nmod_poly_divrem(hi.v,lo.v,digit.v,f.v);for(int j=0;j<=lo.deg();++j)if(lo.at(j))out[prefix+pkey(f)+":"+std::to_string(power)+":"+std::to_string(j)]=lo.at(j);digit=std::move(hi);}}
 }
};
struct TriangularIntegrator {
 int r;std::vector<std::vector<ZR>>A;Factors factors;std::vector<std::unique_ptr<ScalarHermite>>scalar;EndpointGerm endpoint;Quotient completion;long long steps=0;int completion_count=0;
 TriangularIntegrator(std::vector<std::vector<ZR>>A,int germ_order):r(A.size()),A(std::move(A)),endpoint(this->A,germ_order){for(int i=0;i<r;++i){for(int j=0;j<i;++j)if(!this->A[i][j].zero())throw std::runtime_error("FBI connection not sector triangular");scalar.push_back(std::make_unique<ScalarHermite>(this->A[i][i],factors));}}
 C raw(std::vector<ZR>v){C out;for(int i=0;i<r;++i)if(!v[i].zero()){auto [h,g]=scalar[i]->reduce(std::move(v[i]));scalar[i]->coordinates(h,i,out);axpy(out,endpoint.finite_part(i,g),prime-1);if(!g.zero())for(int j=i+1;j<r;++j)if(!A[i][j].zero())v[j]=v[j]-g*A[i][j];++steps;}return out;}
 C integrate(std::vector<ZR>v){return completion.reduce(raw(std::move(v)));}
 C boundary(const std::vector<ZR>&v,int end){if(end==1)return {};C out;for(int i=0;i<r;++i)if(!v[i].zero())axpy(out,endpoint.finite_part(i,v[i]),1);return completion.reduce(std::move(out));}
 void impose_analytic_boundary(const std::vector<int>&vanishing_order){
  if(int(vanishing_order.size())!=r)throw std::runtime_error("boundary order size mismatch");
  for(int i=0;i<r;++i)for(auto&[n,logs]:endpoint.jets[i])for(size_t l=0;l<logs.size();++l)
   if(l||n<vanishing_order[i])completion.insert(logs[l]);
 }
 void complete(int limit){for(int i=0;i<r;++i){auto submit=[&](const ZR&g){std::vector<ZR>v(r);v[i]=g.diff();for(int j=i;j<r;++j)v[j]=v[j]+g*A[i][j];C c=raw(std::move(v));axpy(c,endpoint.finite_part(i,g),1);if(!c.empty()){completion.insert(std::move(c));++completion_count;}};
   auto&a=*scalar[i];for(auto&[f,k]:factors.get(a.A.d)){P res=a.residue(f),df=derivative(f);for(int n=2;n<=limit+1;++n)if(modp(res-scale(df,n-1),f).zero())for(int j=0;j<f.deg();++j){P p;p.set(j,1);submit(ZR(p,ppow(f,n-1)));}}
   U lambda=a.A.n.deg()+1==a.A.d.deg()?nmod_div(a.A.n.at(a.A.n.deg()),a.A.d.at(a.A.d.deg()),mod):0;for(int n=0;n<=limit;++n)if(add(U(n),lambda)==0){P p;p.set(n,1);submit(zp(p));}
  }}
};
