// Polynomial LRR and q-adic coefficient integration. No J reduction/DE input.
#include <flint/nmod_poly.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using json=nlohmann::json;
using U=mp_limb_t;
static U prime;
static nmod_t mod;
struct P {
 nmod_poly_t v;
 P(){nmod_poly_init(v,prime);}
 explicit P(U c):P(){nmod_poly_set_coeff_ui(v,0,c);}
 P(const P& b):P(){nmod_poly_set(v,b.v);}
 P(P&& b) noexcept:P(){nmod_poly_swap(v,b.v);}
 P& operator=(const P& b){nmod_poly_set(v,b.v);return *this;}
 P& operator=(P&& b) noexcept {nmod_poly_swap(v,b.v);return *this;}
 ~P(){nmod_poly_clear(v);}
 long deg()const{return nmod_poly_degree(v);}
 U at(long i)const{return i<0?0:nmod_poly_get_coeff_ui(v,i);}
 void set(long i,U x){nmod_poly_set_coeff_ui(v,i,x);}
 bool zero()const{return nmod_poly_is_zero(v);}
 P& operator+=(const P& b){nmod_poly_add(v,v,b.v);return *this;}
 P& operator-=(const P& b){nmod_poly_sub(v,v,b.v);return *this;}
};
P operator+(P a,const P& b){return a+=b;}
P operator-(P a,const P& b){return a-=b;}
P operator*(const P& a,const P& b){P c;nmod_poly_mul(c.v,a.v,b.v);return c;}
P scale(const P& a,U x){P c;nmod_poly_scalar_mul_nmod(c.v,a.v,x);return c;}
P derivative(const P& a){P c;nmod_poly_derivative(c.v,a.v);return c;}
P quotient(const P& a,const P& b){P q,r;nmod_poly_divrem(q.v,r.v,a.v,b.v);if(!r.zero())throw std::runtime_error("nonexact polynomial division");return q;}
P readp(const json& a){P p;for(size_t i=0;i<a.size();++i)p.set(i,a[i].get<U>()%prime);return p;}
using Row=std::vector<P>;
using Mat=std::vector<Row>;
using Jet=std::vector<P>;
Jet readj(const json& a){Jet j;for(auto& x:a)j.push_back(readp(x));return j;}
U add(U a,U b){return nmod_add(a,b,mod);} U sub(U a,U b){return nmod_sub(a,b,mod);} U mul(U a,U b){return nmod_mul(a,b,mod);}
struct PartialInverse {
 std::vector<std::vector<U>> m;
 std::vector<std::vector<U>> kernel;
 size_t rank=0;
 explicit PartialInverse(std::vector<std::vector<U>> a){
  size_t n=a.size(),row=0;m.assign(n,std::vector<U>(n));
  for(size_t i=0;i<n;++i){a[i].resize(2*n);a[i][n+i]=1;}
  std::vector<size_t> pivots;
  for(size_t col=0;col<n && row<n;++col){size_t v=row;while(v<n&&!a[v][col])++v;if(v==n)continue;
   std::swap(a[row],a[v]);U inv=nmod_inv(a[row][col],mod);
   for(U& x:a[row])x=mul(x,inv);
   for(size_t i=0;i<n;++i)if(i!=row&&a[i][col]){U c=a[i][col];for(size_t j=0;j<2*n;++j)a[i][j]=sub(a[i][j],mul(c,a[row][j]));}
   pivots.push_back(col);++row;
  }
  for(size_t i=0;i<pivots.size();++i)for(size_t j=0;j<n;++j)m[pivots[i]][j]=a[i][n+j];
  rank=pivots.size();
  for(size_t col=0;col<n;++col)if(std::find(pivots.begin(),pivots.end(),col)==pivots.end()){std::vector<U>v(n);v[col]=1;for(size_t i=0;i<pivots.size();++i)v[pivots[i]]=a[i][col]?prime-a[i][col]:0;kernel.push_back(std::move(v));}
 }
 std::vector<U> solve(const std::vector<U>& b)const{
  std::vector<U>x(m.size());for(size_t i=0;i<m.size();++i)for(size_t j=0;j<m.size();++j)x[i]=add(x[i],mul(m[i][j],b[j]));return x;
 }
};
struct PoleData {
 PartialInverse inverse;
 std::vector<std::vector<U>> transfer,residual;
 explicit PoleData(PartialInverse x):inverse(std::move(x)){}
};
struct Qadic {
 P q;std::map<int,P> powers;
 explicit Qadic(P q):q(std::move(q)){powers.emplace(0,P(1));powers.emplace(1,this->q);}
 const P& pow(int n){auto it=powers.find(n);if(it!=powers.end())return it->second;P a=pow(n/2);P b=a*a;if(n%2)b=b*q;return powers.emplace(n,std::move(b)).first->second;}
 void split(const P& p,int n,std::vector<P>& out,int start){
  if(p.zero())return;
  if(n==1){out[start]=p;return;}
  int h=n/2;P hi,lo;const P& divisor=pow(h);
  if(divisor.deg()>32&&p.deg()>divisor.deg()+32)nmod_poly_divrem_newton(hi.v,lo.v,p.v,divisor.v);
  else nmod_poly_divrem(hi.v,lo.v,p.v,divisor.v);
  split(lo,h,out,start);split(hi,n-h,out,start+h);
 }
 std::vector<P> digits(const P& p,int minimum=0){int n=std::max(minimum,int(p.deg()/q.deg()+1));n=std::max(1,n);std::vector<P>d(n);split(p,n,d,0);return d;}
 P compose(const std::vector<P>& d,int a,int b){if(a>=b)return P();if(a+1==b)return d[a];int h=(a+b)/2;return compose(d,a,h)+pow(h-a)*compose(d,h,b);}
};
struct Reducer {
 int r,s;P q,dq;Mat B;Qadic qa;std::map<int,PoleData> poles;std::map<int,PartialInverse> infinity;
 std::vector<std::string> modes;
 long long pole_steps=0,infinity_steps=0,checks=0;
 double digit_seconds=0,pole_seconds=0,certificate_seconds=0;
 bool completing=false;
 using Classes=std::map<std::string,U>;
 std::vector<std::pair<std::string,Classes>> relations;
 json completion_certificates=json::array();
 Reducer(P q,Mat B,std::vector<std::string> modes):r(B.size()),s(q.deg()),q(q),dq(derivative(q)),B(std::move(B)),qa(q),modes(std::move(modes)){}
 const PoleData& pole(int k){auto it=poles.find(k);if(it!=poles.end())return it->second;
  std::vector<std::vector<U>>a(r*s,std::vector<U>(r*s));
  for(int i=0;i<r;++i)for(int j=0;j<s;++j)for(int t=0;t<r;++t){P x=B[i][t];if(i==t)x-=scale(dq,(k-1)%prime);P y;nmod_poly_shift_left(y.v,x.v,j);nmod_poly_rem(y.v,y.v,q.v);for(int d=0;d<s;++d)a[t*s+d][i*s+j]=y.at(d);}
  PoleData data(PartialInverse(std::move(a)));
  int size=r*s;data.transfer.assign(size,std::vector<U>(size));
  if(data.inverse.rank!=size_t(size))data.residual.assign(size,std::vector<U>(size));
  // Precompose the entire local pole elimination, including the quotient update.
  for(int col=0;col<size;++col){Row aa(r);for(int i=0;i<r;++i)for(int j=0;j<s;++j)aa[i].set(j,data.inverse.m[i*s+j][col]);
   for(int j=0;j<r;++j){P num=q*derivative(aa[j])-scale(dq*aa[j],U(k-1)%prime);for(int i=0;i<r;++i)num+=aa[i]*B[i][j];P hi,lo;nmod_poly_divrem(hi.v,lo.v,num.v,q.v);
    for(int d=0;d<s;++d){data.transfer[j*s+d][col]=hi.at(d);if(!data.residual.empty())data.residual[j*s+d][col]=sub(j*s+d==col?1:0,lo.at(d));}
   }
  }
  return poles.emplace(k,std::move(data)).first->second;
 }
 const PartialInverse& inf(int n){auto it=infinity.find(n);if(it!=infinity.end())return it->second;
  std::vector<std::vector<U>>a(r,std::vector<U>(r));for(int i=0;i<r;++i)for(int j=0;j<r;++j)a[i][j]=add(B[j][i].at(s-1),i==j?U(n)%prime:0);
  return infinity.emplace(n,PartialInverse(std::move(a))).first->second;
 }
 static bool preferred(const std::string&a,const std::string&b){bool pa=a.rfind("period:",0)==0,pb=b.rfind("period:",0)==0;if(pa!=pb)return pa;return a>b;}
 Classes canonicalize(Classes c)const{for(auto&[pivot,row]:relations){auto it=c.find(pivot);if(it==c.end()||!it->second)continue;U x=it->second;for(auto&[key,v]:row)c[key]=sub(c[key],mul(x,v));}for(auto it=c.begin();it!=c.end();)if(!it->second)it=c.erase(it);else ++it;return c;}
 void add_relation(Classes c){c=canonicalize(std::move(c));if(c.empty())return;std::string pivot=c.begin()->first;for(auto&[key,v]:c)if(preferred(key,pivot))pivot=key;U inv=nmod_inv(c.at(pivot),mod);for(auto&[key,v]:c)v=mul(v,inv);for(auto&entry:relations){U x=entry.second.count(pivot)?entry.second.at(pivot):0;if(x)for(auto&[key,v]:c)entry.second[key]=sub(entry.second[key],mul(x,v));}relations.push_back({pivot,std::move(c)});std::sort(relations.begin(),relations.end(),[](const auto&a,const auto&b){return preferred(a.first,b.first);});}
 void complete_relations(int limit){completing=true;
  auto submit=[&](const Row&a,int denominator_power,const std::string&kind,int degree){Row input(r);for(int j=0;j<r;++j){input[j]=q*derivative(a[j])-scale(dq*a[j],U(denominator_power)%prime);for(int i=0;i<r;++i)input[j]+=a[i]*B[i][j];}Classes c=reduce(input,denominator_power+1,true);
   for(int point=0;point<2;++point)if(modes[point]=="regular"){U qe=nmod_poly_evaluate_nmod(q.v,point);if(!qe)throw std::runtime_error("completion requires regular common denominator at explicit endpoint");U den=nmod_pow_ui(qe,denominator_power,mod);for(int j=0;j<r;++j){U v=nmod_div(nmod_poly_evaluate_nmod(a[j].v,point),den,mod);if(point==0&&v)v=prime-v;std::string key="endpoint:"+std::to_string(point)+":"+std::to_string(j);c[key]=sub(c[key],v);}}
   for(auto it=c.begin();it!=c.end();)if(!it->second)it=c.erase(it);else++it;
   if(!c.empty()){json numer=json::array();for(auto&p:a){json v=json::array();for(int j=0;j<=p.deg();++j)v.push_back(p.at(j));numer.push_back(v);}completion_certificates.push_back({{"kind",kind},{"degree",degree},{"g_numerators",numer},{"g_q_denominator_power",denominator_power},{"raw_relation",c}});add_relation(std::move(c));}
  };
  for(int n=0;n<=limit;++n){auto kernel=inf(n).kernel;for(auto&v:kernel){Row a(r);for(int j=0;j<r;++j)a[j].set(n,v[j]);submit(a,0,"infinity_kernel",n);}}
  for(int k=2;k<=limit+1;++k){auto kernel=pole(k).inverse.kernel;for(auto&v:kernel){Row a(r);for(int i=0;i<r;++i)for(int j=0;j<s;++j)a[i].set(j,v[i*s+j]);submit(a,k-1,"finite_pole_kernel",k);}}
  completing=false;pole_steps=infinity_steps=checks=0;digit_seconds=pole_seconds=certificate_seconds=0;
 }
 void record(Classes& c,int k,int i,int j,U x){if(!x)return;std::string key="period:"+std::to_string(k)+":"+std::to_string(i)+":"+std::to_string(j);c[key]=add(c[key],x);}
 Classes reduce(const Row& input,int K,bool verify){
  auto stamp=std::chrono::steady_clock::now();
  std::vector<std::vector<P>> digits;for(auto& x:input)digits.push_back(qa.digits(x,K));
  digit_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-stamp).count();stamp=std::chrono::steady_clock::now();
  std::vector<std::vector<P>> gd,hd;if(verify){gd.assign(r,std::vector<P>(K));hd.assign(r,std::vector<P>(K));}
  Row gpoly(r),hpoly(r);std::vector<U> endpoint_num(2*r);Classes classes;
  std::vector<U>rhs(r*s);bool boundary=modes[0]=="regular"||modes[1]=="regular";
  for(int k=K;k>=2;--k){int idx=K-k;for(int i=0;i<r;++i)for(int j=0;j<s;++j)rhs[i*s+j]=digits[i][idx].at(j);
   const auto& map=pole(k);std::vector<U>x;if(verify||boundary)x=map.inverse.solve(rhs);
   for(int j=0;j<r;++j){P a,obs;
    for(int d=0;d<s;++d){int row=j*s+d;U change=0;for(int c=0;c<r*s;++c)change=add(change,mul(map.transfer[row][c],rhs[c]));digits[j][idx+1].set(d,sub(digits[j][idx+1].at(d),change));
     if(!map.residual.empty()){U left=0;for(int c=0;c<r*s;++c)left=add(left,mul(map.residual[row][c],rhs[c]));record(classes,k,j,d,left);if(verify)obs.set(d,left);}
     if(verify||boundary)a.set(d,x[row]);
    }
    if(verify){gd[j][idx]=a;hd[j][idx]=obs;}
    for(int e=0;e<2;++e)if(modes[e]=="regular"){
     U qe=nmod_poly_evaluate_nmod(q.v,e);if(!qe)throw std::runtime_error("regular endpoint lies in common q; use a different gauge or endpoint backend");
     endpoint_num[e*r+j]=add(endpoint_num[e*r+j],mul(nmod_pow_ui(qe,idx,mod),nmod_poly_evaluate_nmod(a.v,e)));
    }
   }++pole_steps;
  }
  pole_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-stamp).count();
  Row p(r);for(int i=0;i<r;++i)p[i]=qa.compose(digits[i],K-1,digits[i].size());
  while(true){int d=-1;for(auto& x:p)d=std::max(d,int(x.deg()));if(d<s-1)break;int n=d-s+1;
   std::vector<U>rhs(r);for(int i=0;i<r;++i)rhs[i]=p[i].at(d);auto a=inf(n).solve(rhs);
   // Sparse shifted update: work per step depends on q and r, not the current degree.
   for(int j=0;j<r;++j){for(int h=0;h<s;++h){U c=0;for(int i=0;i<r;++i)c=add(c,mul(a[i],B[i][j].at(h)));p[j].set(n+h,sub(p[j].at(n+h),c));}
    if(n)for(int h=0;h<=s;++h){int ix=n-1+h;U c=mul(mul(a[j],U(n)%prime),q.at(h));p[j].set(ix,sub(p[j].at(ix),c));}
    U obs=p[j].at(d);record(classes,1,j,d,obs);if(verify)hpoly[j].set(d,obs);p[j].set(d,0);
    if(verify)gpoly[j].set(n,a[j]);
    for(int e=0;e<2;++e)if(modes[e]=="regular"){
     U term=(e==0&&n)?0:a[j]; // 0^0=1; 1^n=1
     U qe=nmod_poly_evaluate_nmod(q.v,e);endpoint_num[e*r+j]=add(endpoint_num[e*r+j],mul(term,nmod_pow_ui(qe,K-1,mod)));
    }
   }++infinity_steps;
  }
  for(int i=0;i<r;++i)for(int j=0;j<=p[i].deg();++j)record(classes,1,i,j,p[i].at(j));
  stamp=std::chrono::steady_clock::now();
  if(verify){Row G(r),H(r);for(int i=0;i<r;++i){G[i]=qa.compose(gd[i],0,K)+gpoly[i]*qa.pow(K-1);H[i]=qa.compose(hd[i],0,K)+(p[i]+hpoly[i])*qa.pow(K-1);}
   for(int j=0;j<r;++j){P deriv=q*derivative(G[j])-scale(dq*G[j],U(K-1)%prime);for(int i=0;i<r;++i)deriv+=G[i]*B[i][j];if(!(input[j]-H[j]-deriv).zero())throw std::runtime_error("global polynomial integration certificate failed");}++checks;
  }
  certificate_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-stamp).count();
  for(int e=0;e<2;++e)if(modes[e]=="regular")for(int i=0;i<r;++i){U den=nmod_poly_evaluate_nmod(qa.pow(K-1).v,e);U v=nmod_div(endpoint_num[e*r+i],den,mod);if(e==0&&v)v=prime-v;if(v)classes["endpoint:"+std::to_string(e)+":"+std::to_string(i)]=v;}
  for(auto it=classes.begin();it!=classes.end();)if(!it->second)it=classes.erase(it);else++it;
  return completing?classes:canonicalize(std::move(classes));
 }
};
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
#include "fbi.hpp"
int main(int argc,char**argv){try{
 if(argc<3)throw std::runtime_error("usage: j_series INPUT OUTPUT [--order N] [--verify-all]");
 json in;{std::ifstream f(argv[1]);f>>in;}
 for(int i=3;i<argc;++i){std::string arg=argv[i];if(arg=="--verify-all")in["verify_all"]=true;else if(arg=="--boundaries")in["include_z_boundaries"]=true;else if(arg=="--max-dot"&&i+1<argc)in["max_dot"]=std::stoi(argv[++i]);else if(arg=="--prime"&&i+1<argc)in["prime"]=std::stoul(argv[++i]);else if(arg=="--t0"&&i+1<argc)in["t0"]=argv[++i];else if(arg=="--order"&&i+1<argc)in["order"]=std::stoi(argv[++i]);else throw std::runtime_error("unknown/incomplete option: "+arg);}
 auto process_start=Clock::now();
 prime=in.at("prime").get<U>();nmod_init(&mod,prime);if(in.contains("kind"))in=prepare(in);int r=in.at("rank"),order=in.at("order");if(order<0||U(order+1)>=prime)throw std::runtime_error("invalid order/prime");
 P q=readp(in.at("q"));Mat B(r,Row(r));for(int i=0;i<r;++i)for(int j=0;j<r;++j)B[i][j]=readp(in["B"][i][j]);
 Reducer reducer(q,B,in.at("endpoints").get<std::vector<std::string>>());
 reducer.complete_relations(in.value("completion_degree_limit",16));
 Jet qt=readj(in.at("qt"));std::vector<std::vector<Jet>> nt(r,std::vector<Jet>(r));for(int i=0;i<r;++i)for(int j=0;j<r;++j)nt[i][j]=readj(in["nt"][i][j]);
 std::vector<Jet> den;std::vector<std::vector<Jet>> numer;std::vector<std::string> labels;
 for(auto& t:in.at("targets")){labels.push_back(t.at("label"));den.push_back(readj(t.at("den")));std::vector<Jet> v;for(auto& x:t.at("numer"))v.push_back(readj(x));numer.push_back(std::move(v));}
 int e=1;P H=q;while(true){bool ok=true;P rem;for(auto& d:den){nmod_poly_rem(rem.v,H.v,d[0].v);ok&=rem.zero();}nmod_poly_rem(rem.v,H.v,qt[0].v);ok&=rem.zero();if(ok)break;if(++e>100)throw std::runtime_error("denominator outside common q");H=H*q;}
 int maxlag=qt.size();for(auto& row:nt)for(auto& j:row)maxlag=std::max(maxlag,int(j.size()+1));for(auto& x:numer)for(auto& j:x)maxlag=std::max(maxlag,int(j.size()+1));for(auto& x:den)maxlag=std::max(maxlag,int(x.size()));
 std::vector<P> hp(maxlag+2);hp[0]=P(1);for(size_t i=1;i<hp.size();++i)hp[i]=hp[i-1]*H;
 P hq=quotient(H,qt[0]);for(int i=0;i<r;++i)for(int j=0;j<r;++j)for(size_t d=0;d<nt[i][j].size();++d)nt[i][j][d]=nt[i][j][d]*hp[d]*hq;
 for(size_t d=1;d<qt.size();++d)qt[d]=qt[d]*hp[d-1]*hq;
 int T=labels.size();std::vector<std::vector<Row>> vbuf(T);std::vector<int> vsize(T);
 for(int a=0;a<T;++a){P hd=quotient(H,den[a][0]);for(int i=0;i<r;++i)for(size_t d=0;d<numer[a][i].size();++d)numer[a][i][d]=numer[a][i][d]*hp[d]*hd;for(size_t d=1;d<den[a].size();++d)den[a][d]=den[a][d]*hp[d-1]*hd;
  vsize[a]=std::max(1,int(den[a].size()));vbuf[a].assign(vsize[a],Row(r));}
 std::vector<Mat> fbuf(maxlag,Mat(r,Row(r)));for(int i=0;i<r;++i)fbuf[0][i][i]=P(1);
 json output;output["coordinates"]=json::object();for(auto& l:labels)output["coordinates"][l]=json::array();
 bool boundaries=in.value("include_z_boundaries",false);
 if(boundaries)for(auto&l:labels)for(int point=0;point<2;++point)output["coordinates"]["Bz"+std::to_string(point)+l.substr(1)]=json::array();
 double ftime=0,vtime=0,itime=0;long long maxdegree=0;auto total=Clock::now();
 for(int n=0;n<=order;++n){auto t=Clock::now();if(n){Mat f(r,Row(r));for(int i=0;i<r;++i)for(int j=0;j<r;++j){P x;for(int a=0;a<r;++a)for(size_t d=0;d<nt[i][a].size()&&int(d)<n;++d)if(!nt[i][a][d].zero())x+=nt[i][a][d]*fbuf[(n-1-d)%maxlag][a][j];for(size_t d=1;d<qt.size()&&int(d)<=n;++d)if(n!=int(d))x-=scale(qt[d]*fbuf[(n-d)%maxlag][i][j],U(n-d)%prime);f[i][j]=scale(x,nmod_inv(n,mod));}fbuf[n%maxlag]=std::move(f);}ftime+=seconds(t);
  for(int a=0;a<T;++a){t=Clock::now();Row row(r);for(int j=0;j<r;++j){P x;for(int i=0;i<r;++i)for(size_t d=0;d<numer[a][i].size()&&int(d)<=n;++d)if(!numer[a][i][d].zero())x+=numer[a][i][d]*fbuf[(n-d)%maxlag][i][j];for(size_t d=1;d<den[a].size()&&int(d)<=n;++d)x-=den[a][d]*vbuf[a][(n-d)%vsize[a]][j];row[j]=std::move(x);maxdegree=std::max(maxdegree,(long long)row[j].deg());}vbuf[a][n%vsize[a]]=row;vtime+=seconds(t);
   t=Clock::now();bool verify=in.value("verify_all",false)||n<5||n==order;auto c=reducer.reduce(row,e*(n+1),verify);output["coordinates"][labels[a]].push_back(c);itime+=seconds(t);
   if(boundaries)for(int point=0;point<2;++point){Reducer::Classes b;
    if(reducer.modes[point]=="regular"){U hp=nmod_poly_evaluate_nmod(H.v,point);if(!hp)throw std::runtime_error("boundary evaluation requires nonzero common denominator");U inv=nmod_inv(nmod_pow_ui(hp,n+1,mod),mod);for(int j=0;j<r;++j){U value=mul(nmod_poly_evaluate_nmod(row[j].v,point),inv);if(value)b["endpoint:"+std::to_string(point)+":"+std::to_string(j)]=value;}}
    output["coordinates"]["Bz"+std::to_string(point)+labels[a].substr(1)].push_back(reducer.canonicalize(std::move(b)));
   }
  }
  if(n%100==0||n==order)std::cerr<<"coefficient "<<n<<"/"<<order<<" elapsed="<<seconds(total)<<"s\n";
 }
 output["timings"]={{"fbi_jet_seconds",ftime},{"target_jet_seconds",vtime},{"integration_seconds",itime},{"native_total_seconds",seconds(total)},{"preparation_seconds",in.value("preparation_seconds",0.)},{"qadic_seconds",reducer.digit_seconds},{"pole_elimination_seconds",reducer.pole_seconds},{"certificate_seconds",reducer.certificate_seconds},{"process_compute_seconds",seconds(process_start)}};
 output["prime"]=prime;output["order"]=order;output["fbi_rank"]=r;output["metadata"]=in.value("metadata",json::object());output["t0"]=in.value("t0",json());output["endpoints"]=in.at("endpoints");output["include_z_boundaries"]=boundaries;
 output["operations"]={{"pole_steps",reducer.pole_steps},{"infinity_steps",reducer.infinity_steps},{"certificate_checks",reducer.checks},{"coefficient_integrations",T*(order+1)},{"max_numerator_degree",maxdegree},{"denominator_stride",e}};
 output["period_completion"]={{"degree_limit",in.value("completion_degree_limit",16)},{"independent_relations",reducer.relations.size()},{"certificates",reducer.completion_certificates}};
 {std::filesystem::path path(argv[2]);if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());std::ofstream f(path);if(!f)throw std::runtime_error("cannot open coefficient output");f<<output.dump();if(!f)throw std::runtime_error("cannot write coefficient output");}return 0;
 }catch(const std::exception&e){std::cerr<<"j_integral native: "<<e.what()<<"\n";return 1;}}
