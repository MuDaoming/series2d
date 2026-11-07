// Direct z-IBP at sampled t values, followed by rational reconstruction.
// No Taylor coefficients, J reduction, or J differential equation are inputs.
#define main shared_series_unused_main
#include "polynomial.cpp"
#undef main
#include <memory>
#include <sstream>
using C=Reducer::Classes;
void axpy(C&v,const C&w,U c){if(!c)return;for(auto&[k,x]:w){U y=add(v[k],mul(c,x));if(y)v[k]=y;else v.erase(k);}}
C scaled(C v,U c){for(auto&[k,x]:v)x=mul(x,c);return v;}
struct Span {
 struct Entry {C row,weights;};std::map<std::string,Entry> rows;
 bool insert(C v,int index){C w={{std::to_string(index),1}};for(auto&[key,e]:rows){auto it=v.find(key);if(it!=v.end()){U a=prime-it->second;axpy(v,e.row,a);axpy(w,e.weights,a);}}
  if(v.empty())return false;
  U inv=nmod_inv(v.begin()->second,mod);std::string key=v.begin()->first;rows.emplace(key,Entry{scaled(v,inv),scaled(w,inv)});return true;}
 std::vector<U> solve(C v,int size)const{C answer;for(auto&[key,e]:rows){auto it=v.find(key);if(it!=v.end()){U c=it->second;axpy(v,e.row,prime-c);axpy(answer,e.weights,c);}}if(!v.empty())throw std::runtime_error("chosen J basis does not span target/derivative");std::vector<U>a(size);for(auto&[k,x]:answer)a.at(std::stoi(k))=x;return a;}
};
struct Quotient {
 std::map<std::string,C> rows;
 C reduce(C v)const{for(auto&[key,row]:rows){auto it=v.find(key);if(it!=v.end())axpy(v,row,prime-it->second);}return v;}
 void insert(C v){v=reduce(std::move(v));if(v.empty())return;U inv=nmod_inv(v.begin()->second,mod);std::string key=v.begin()->first;rows.emplace(key,scaled(v,inv));}
};
std::pair<P,P> specialize(const RF&a,U point){MP n,d;nmod_mpoly_evaluate_one_ui(n.v,a.n.v,0,point,mc);nmod_mpoly_evaluate_one_ui(d.v,a.d.v,0,point,mc);P nn=univariate(n),dd=univariate(d),g;if(dd.zero())throw std::runtime_error("singular specialization");nmod_poly_gcd(g.v,nn.v,dd.v);return {quotient(nn,g),quotient(dd,g)};}
using RationalRow=std::vector<std::pair<P,P>>;
U ev(const std::pair<P,P>&a,U z){U den=nmod_poly_evaluate_nmod(a.second.v,z);if(!den)throw std::runtime_error("singular endpoint");return nmod_div(nmod_poly_evaluate_nmod(a.first.v,z),den,mod);}
P ppow(P x,int n){P y(1);while(n){if(n&1)y=y*x;n/=2;if(n)x=x*x;}return y;}
int required_power(const P&q,const P&d){P x(1),rem;for(int k=0;k<=100;++k){nmod_poly_rem(rem.v,x.v,d.v);if(rem.zero())return k;x=x*q;}throw std::runtime_error("denominator not covered by q");}
struct PointIntegrator {
 int r;P q;Mat B;std::vector<std::string>modes;std::unique_ptr<Reducer> hermite;Quotient dense;
 int K=0,N=0,ansatz_generators=0;long long integrations=0;
 PointIntegrator(P q,Mat B,std::vector<std::string> modes,const std::string&method,const std::vector<RationalRow>&rows,int extra,int degree):r(B.size()),q(q),B(B),modes(modes){
  if(method=="hermite"){hermite=std::make_unique<Reducer>(q,B,modes);hermite->complete_relations(degree);return;}
  for(auto&row:rows)for(auto&x:row)K=std::max(K,required_power(q,x.second));
  K=std::max(1,K)+extra;P D=ppow(q,K);
  N=K*q.deg()+degree-1;for(auto&row:rows)for(auto&x:row)N=std::max(N,int(quotient(D*x.first,x.second).deg()));
  P dq=derivative(q);for(int i=0;i<r;++i)for(int j=0;j<=N-int(q.deg())+1;++j){P a;a.set(j,1);C relation;
   for(int c=0;c<r;++c){P v=a*B[i][c];if(i==c)v+=q*derivative(a)-scale(dq*a,U(K-1)%prime);for(int d=0;d<=v.deg();++d)if(v.at(d))relation[period(c,d)]=v.at(d);}
   for(int e=0;e<2;++e)if(modes[e]=="regular"){U den=nmod_pow_ui(nmod_poly_evaluate_nmod(q.v,e),K-1,mod);if(!den)throw std::runtime_error("ansatz regular endpoint singular");U aev=j==0?1:U(e);U value=nmod_div(aev,den,mod);if(e&&value)value=prime-value;if(value)relation["endpoint:"+std::to_string(e)+":"+std::to_string(i)]=value;}
   dense.insert(std::move(relation));++ansatz_generators;
  }
 }
 std::string period(int i,int d)const{std::string x=std::to_string(N-d);return "a:"+std::string(8-x.size(),'0')+x+":"+std::to_string(i);}
 C integrate(const RationalRow&row){++integrations;
  if(hermite){int k=1;for(auto&x:row)k=std::max(k,required_power(q,x.second));P D=ppow(q,k);Row nums;for(auto&x:row)nums.push_back(quotient(D*x.first,x.second));return hermite->reduce(nums,k,true);}
  P D=ppow(q,K);C v;for(int i=0;i<r;++i){P p=quotient(D*row[i].first,row[i].second);if(p.deg()>N)throw std::runtime_error("ansatz degree too small");for(int d=0;d<=p.deg();++d)if(p.at(d))v[period(i,d)]=p.at(d);}return dense.reduce(std::move(v));
 }
 C endpoint(const RationalRow&row,int e){C v;if(modes[e]=="regular")for(int i=0;i<r;++i)if(U x=ev(row[i],e))v["endpoint:"+std::to_string(e)+":"+std::to_string(i)]=x;return hermite?hermite->canonicalize(std::move(v)):dense.reduce(std::move(v));}
};
struct Fraction {
 P n,d;
 Fraction():d(1){}
 void normal(){if(n.zero()){d=P(1);return;}P g;nmod_poly_gcd(g.v,n.v,d.v);n=quotient(n,g);d=quotient(d,g);U inv=nmod_inv(d.at(d.deg()),mod);n=scale(n,inv);d=scale(d,inv);}
 json value()const{return {{"numerator_tau",polyjson(n)},{"denominator_tau",polyjson(d)}};}
};
struct Interpolator {
 P modulus;std::vector<P>basis;
 explicit Interpolator(const std::vector<U>&xs){modulus=P(1);for(U x:xs){P f;f.set(0,x?prime-x:0);f.set(1,1);modulus=modulus*f;}for(U x:xs){P f;f.set(0,x?prime-x:0);f.set(1,1);P a=quotient(modulus,f);basis.push_back(scale(a,nmod_inv(nmod_poly_evaluate_nmod(a.v,x),mod)));}}
 bool reconstruct(const std::vector<U>&values,const std::vector<U>&xs,Fraction&out)const{
  P a;for(size_t i=0;i<basis.size();++i)a+=scale(basis[i],values[i]);int bound=(basis.size()-1)/2;
  if(a.zero()){out=Fraction();}else{P r0=modulus,r1=a,s0,s1(1);while(r1.deg()>bound){P quo,rem;nmod_poly_divrem(quo.v,rem.v,r0.v,r1.v);P next=s0-quo*s1;r0=std::move(r1);r1=std::move(rem);s0=std::move(s1);s1=std::move(next);if(r1.zero())return false;}if(s1.zero()||s1.deg()>bound)return false;out.n=r1;out.d=s1;out.normal();}
  for(size_t i=0;i<xs.size();++i){U d=nmod_poly_evaluate_nmod(out.d.v,xs[i]);if(!d||nmod_poly_evaluate_nmod(out.n.v,xs[i])!=mul(values[i],d))return false;}return true;
 }
};
int weight(const std::string&s){int sum=0,n=0;for(char c:s){if(std::isdigit(c))n=n*10+c-'0';else{sum+=n;n=0;}}return sum+n;}
std::string parent_j(const std::string&s){return "J"+s.substr(s.find('{'));}
json all_targets(json cfg,int dot){int n=cfg.at("propagators").size();cfg["targets"]=json::array();Nu nu(n,1);std::function<void(int,int)>visit=[&](int i,int left){if(i==n){cfg["targets"].push_back(nu);return;}for(int k=0;k<=left;++k){nu[i]=1+k;visit(i+1,left-k);}nu[i]=1;};visit(0,dot);return cfg;}
int row_degree(const std::vector<Fraction>&v){P common(1);for(auto&x:v){P g;nmod_poly_gcd(g.v,common.v,x.d.v);common=quotient(common,g)*x.d;}int b=common.deg();for(auto&x:v)b=std::max(b,int(quotient(common*x.n,x.d).deg()));return b;}
#ifndef DIRECT_J_NO_MAIN
int main(int argc,char**argv){try{
 if(argc<3)throw std::runtime_error("usage: direct_j FAMILY OUTPUT [--method hermite|ansatz --max-dot 2 --max-points 193 --degree-limit 16 --extra-poles 0]");
 json cfg;{std::ifstream f(argv[1]);f>>cfg;}std::string method="hermite";int dot=2,maxpoints=193,limit=16,extra=0;
 for(int i=3;i<argc;++i){std::string a=argv[i];if(i+1==argc)throw std::runtime_error("missing option value");std::string v=argv[++i];if(a=="--method")method=v;else if(a=="--max-dot")dot=std::stoi(v);else if(a=="--max-points")maxpoints=std::stoi(v);else if(a=="--degree-limit")limit=std::stoi(v);else if(a=="--extra-poles")extra=std::stoi(v);else throw std::runtime_error("unknown option");}
 if(method!="hermite"&&method!="ansatz")throw std::runtime_error("unknown direct backend");
 if(dot<0||maxpoints<17||limit<0||extra<0)throw std::runtime_error("invalid bounds");
 auto start=Clock::now();prime=cfg.at("prime");nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);cfg=all_targets(cfg,dot);NativeFamily family(cfg);RMat At,Az;std::map<std::string,RRow>raw;family.model(At,Az,raw);int r=At.size();
 std::map<std::string,RRow>derivs;for(auto&[name,row]:raw){RRow dr(r);for(int i=0;i<r;++i){dr[i]=row[i].diff(0);for(int j=0;j<r;++j)dr[i]+=row[j]*At[j][i];}derivs[name]=std::move(dr);}
 std::vector<std::string>labels;for(auto&[name,row]:raw){labels.push_back(name);for(int e=0;e<2;++e)labels.push_back("Bz"+std::to_string(e)+name.substr(1));}
 std::sort(labels.begin(),labels.end(),[](auto&a,auto&b){bool ba=a[0]=='B',bb=b[0]=='B';if(ba!=bb)return ba;int ca=weight(a),cb=weight(b);return ca!=cb?ca<cb:a<b;});
 U t0=constant(readr(cfg.at("t0")));std::vector<std::string>modes={cfg["endpoints"]["0"],cfg["endpoints"]["1"]};
 std::vector<std::string>masters,nonzero,zero;std::vector<U>xs,absolute_points;std::vector<std::vector<U>>samples;int skipped=0,attempted=0;long long integration_count=0,certificates=0;double probe_seconds=0,reconstruct_seconds=0;json diagnostics;
 auto prep_seconds=seconds(start);
 std::vector<Fraction>fits;bool done=false;
 for(int a=0;int(xs.size())<maxpoints&&!done;++a){U point=U(17+a)%prime;if(point==0||point==1)continue;if(a>maxpoints+100)throw std::runtime_error("too many singular sample points");++attempted;auto stamp=Clock::now();
  try{
   std::map<std::string,RationalRow>v,dv;P q(1);auto include=[&](const std::pair<P,P>&x){P g;nmod_poly_gcd(g.v,q.v,x.second.v);q=quotient(q,g)*x.second;};
   for(auto&[name,row]:raw){for(auto&x:row){auto s=specialize(x,point);include(s);v[name].push_back(std::move(s));}for(auto&x:derivs.at(name)){auto s=specialize(x,point);include(s);dv[name].push_back(std::move(s));}}
   std::vector<RationalRow> all;for(auto&[name,row]:v)all.push_back(row);for(auto&[name,row]:dv)all.push_back(row);
   std::vector<RationalRow>az(r);for(int i=0;i<r;++i)for(auto&x:Az[i]){auto s=specialize(x,point);include(s);az[i].push_back(std::move(s));}
   P g;nmod_poly_gcd(g.v,q.v,derivative(q).v);q=quotient(q,g);nmod_poly_make_monic(q.v,q.v);if(q.deg()<1)throw std::runtime_error("polynomial-only connection unsupported");Mat B(r,Row(r));for(int i=0;i<r;++i)for(int j=0;j<r;++j){B[i][j]=quotient(q*az[i][j].first,az[i][j].second);if(B[i][j].deg()>=q.deg())throw std::runtime_error("non-Fuchsian connection");}
   PointIntegrator integrator(q,B,modes,method,all,extra,limit);std::map<std::string,C>coords;
   for(auto&[name,row]:v){coords[name]=integrator.integrate(row);for(int e=0;e<2;++e)coords["Bz"+std::to_string(e)+name.substr(1)]=integrator.endpoint(row,e);}
   if(masters.empty()){Span choose;for(auto&name:labels){if(coords.at(name).empty()){zero.push_back(name);continue;}nonzero.push_back(name);if(choose.insert(coords.at(name),masters.size()))masters.push_back(name);}std::cerr<<"direct "<<method<<" selected "<<masters.size()<<" masters\n";}
   Span span;for(int i=0;i<int(masters.size());++i)if(!span.insert(coords.at(masters[i]),i))throw std::runtime_error("rank drop at probe point");
   std::vector<U>sample;for(auto&name:nonzero){auto c=span.solve(coords.at(name),masters.size());sample.insert(sample.end(),c.begin(),c.end());}for(auto&name:zero)if(!coords.at(name).empty())throw std::runtime_error("zero-object classification changed");
   for(auto&name:masters){std::string j=parent_j(name);C derivative=name[0]=='J'?integrator.integrate(dv.at(j)):integrator.endpoint(dv.at(j),name[2]-'0');auto c=span.solve(std::move(derivative),masters.size());sample.insert(sample.end(),c.begin(),c.end());}
   samples.push_back(std::move(sample));xs.push_back(sub(point,t0));absolute_points.push_back(point);integration_count+=integrator.integrations;if(integrator.hermite)certificates+=integrator.hermite->checks;
   diagnostics={{"q_degree",q.deg()},{"FBI_rank",r},{"ansatz_denominator_power",integrator.K},{"ansatz_numerator_degree",integrator.N},{"ansatz_generators",integrator.ansatz_generators},{"ansatz_relation_rank",integrator.dense.rows.size()}};
  }catch(const std::exception&e){++skipped;if(skipped>8)throw;std::cerr<<"skip probe "<<point<<": "<<e.what()<<"\n";}
  probe_seconds+=seconds(stamp);int count=xs.size();if(count<17||(count!=maxpoints&&(count-1)%16!=0))continue;
  stamp=Clock::now();int train=count-8;std::vector<U>training(xs.begin(),xs.begin()+train);Interpolator interp(training);fits.assign(samples[0].size(),Fraction());bool ok=true;
  for(size_t col=0;col<fits.size();++col){std::vector<U>values;for(auto&row:samples)values.push_back(row[col]);if(!interp.reconstruct(values,xs,fits[col])){ok=false;break;}}
  reconstruct_seconds+=seconds(stamp);std::cerr<<"probes="<<count<<" train="<<train<<" reconstruction="<<(ok?"PASS":"needs more points")<<" elapsed="<<seconds(start)<<"s\n";done=ok;
 }
 if(!done)throw std::runtime_error("rational reconstruction did not close within max-points");
 json out={{"method",method=="hermite"?"pointwise Hermite z-IBP and rational reconstruction":"pointwise bounded total-derivative ansatz and rational reconstruction"},{"prime",prime},{"t0",cfg.at("t0")},{"masters",masters},{"master_count",masters.size()},{"tested_objects",nonzero},{"zero_objects",zero},{"reductions",json::object()},{"de",json::object()},{"input",argv[1]},{"max_dot",dot},{"Taylor_order",0},{"sample_points",absolute_points},{"training_points",xs.size()-8},{"heldout_points",8},{"skipped_points",skipped},{"diagnostics",diagnostics},{"degree_limit",limit},{"extra_poles",extra}};
 size_t col=0;for(auto&name:nonzero){json row=json::object();for(auto&m:masters){if(!fits[col].n.zero())row[m]=fits[col].value();++col;}if(std::find(masters.begin(),masters.end(),name)==masters.end())out["reductions"][name]=row;}
 for(auto&name:masters){json row=json::object();std::vector<Fraction>vf;for(auto&m:masters){vf.push_back(fits[col]);if(!fits[col].n.zero())row[m]=fits[col].value();++col;}out["de"][name]={{"rhs",row},{"degree",row_degree(vf)},{"verification_points",8}};}
 out["timings"]={{"preparation_seconds",prep_seconds},{"probe_seconds",probe_seconds},{"reconstruction_seconds",reconstruct_seconds},{"total_compute_seconds",seconds(start)}};out["operations"]={{"coefficient_integrations",integration_count},{"exact_integration_certificate_checks",certificates},{"sample_attempts",attempted}};
 out["caveat"]="Finite-field rational reconstruction with heldout points; independently verify at new points. Bounds on total-derivative ansatz/kernel completion are explicit and not a universal completeness proof.";
 std::filesystem::path output(argv[2]);if(!output.parent_path().empty())std::filesystem::create_directories(output.parent_path());{std::ofstream f(output);if(!f)throw std::runtime_error("cannot open output");f<<out.dump(2);}
 std::cout<<json({{"masters",masters},{"points",xs.size()},{"timings",out["timings"]}}).dump(2)<<"\n";return 0;
 }catch(const std::exception&e){std::cerr<<"direct_j: "<<e.what()<<"\n";return 1;}}

#endif
