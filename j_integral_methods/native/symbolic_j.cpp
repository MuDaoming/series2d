// Full symbolic-t total-derivative ansatz over F_p(t), without sampling or Taylor jets.
#define DIRECT_J_NO_MAIN
#include "direct_j.cpp"
using SV=std::vector<RF>;
int zdegree(const MP&p){int d=-1;for(slong i=0;i<nmod_mpoly_length(p.v,mc);++i){ulong e[2];nmod_mpoly_get_term_exp_ui(e,p.v,i,mc);d=std::max(d,int(e[1]));}return d;}
SV zcoeff(const RF&p){if(zdegree(p.d)>0)throw std::runtime_error("coefficient denominator still depends on z");std::vector<MP>ns(zdegree(p.n)+1);for(slong i=0;i<nmod_mpoly_length(p.n.v,mc);++i){ulong e[2];nmod_mpoly_get_term_exp_ui(e,p.n.v,i,mc);ulong j=e[1];e[1]=0;U c=nmod_mpoly_get_term_coeff_ui(p.n.v,i,mc);nmod_mpoly_set_coeff_ui_ui(ns[j].v,c,e,mc);}SV out;for(auto&n:ns)out.emplace_back(n,p.d);return out;}
RF at_z(const RF&p,int z){MP n,d;nmod_mpoly_evaluate_one_ui(n.v,p.n.v,1,z,mc);nmod_mpoly_evaluate_one_ui(d.v,p.d.v,1,z,mc);return RF(std::move(n),std::move(d));}
struct SymbolicSpan {
 bool track; explicit SymbolicSpan(bool track=true):track(track){}
 struct Entry{SV row,weights;};std::map<int,Entry>rows;
 static void ax(SV&v,const SV&w,const RF&c){if(c.zero())return;v.resize(std::max(v.size(),w.size()));for(size_t i=0;i<w.size();++i)if(!w[i].zero())v[i]+=c*w[i];}
 bool insert(SV v,int index){SV weights;if(track){weights.resize(index+1);weights[index]=RF(1);}for(auto&[k,row]:rows)if(k<int(v.size())&&!v[k].zero()){RF c=-v[k];ax(v,row.row,c);if(track)ax(weights,row.weights,c);}int k=0;while(k<int(v.size())&&v[k].zero())++k;if(k==int(v.size()))return false;RF pivot=v[k];for(auto&x:v)x=x/pivot;for(auto&x:weights)x=x/pivot;rows.emplace(k,Entry{std::move(v),std::move(weights)});return true;}
 SV reduce(SV v)const{for(auto&[k,row]:rows)if(k<int(v.size())&&!v[k].zero())ax(v,row.row,-v[k]);return v;}
 SV solve(SV v,int size)const{SV weights(size);for(auto&[k,row]:rows)if(k<int(v.size())&&!v[k].zero()){RF c=v[k];ax(v,row.row,-c);if(track)ax(weights,row.weights,c);}for(auto&x:v)if(!x.zero())throw std::runtime_error("symbolic J basis does not span derivative");return weights;}
};
Fraction as_fraction(const RF&x,U t0){if(zdegree(x.n)>0||zdegree(x.d)>0)throw std::runtime_error("result still depends on z");Fraction f;auto n=shiftjet(x.n,t0),d=shiftjet(x.d,t0);for(size_t i=0;i<n.size();++i)f.n.set(i,n[i].at(0));f.d=P();for(size_t i=0;i<d.size();++i)f.d.set(i,d[i].at(0));f.normal();return f;}
int main(int argc,char**argv){try{
 if(argc<3)throw std::runtime_error("usage: symbolic_j FAMILY OUTPUT [--max-dot 2 --degree-limit 4]");
 json cfg;{std::ifstream f(argv[1]);f>>cfg;}int dot=2,limit=4;for(int i=3;i<argc;++i){std::string a=argv[i];if(i+1==argc)throw std::runtime_error("missing value");int v=std::stoi(argv[++i]);if(a=="--max-dot")dot=v;else if(a=="--degree-limit")limit=v;else throw std::runtime_error("unknown option");}if(dot<0||limit<2)throw std::runtime_error("invalid bounds");
 auto start=Clock::now();prime=cfg.at("prime");nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);cfg=all_targets(cfg,dot);NativeFamily family(cfg);RMat At,Az;std::map<std::string,RRow>raw,dv;family.model(At,Az,raw);int r=At.size();
 for(auto&[name,row]:raw){RRow d(r);for(int i=0;i<r;++i){d[i]=row[i].diff(0);for(int j=0;j<r;++j)d[i]+=row[j]*At[j][i];}dv[name]=std::move(d);}
 MP qn(1);auto include=[&](const RF&x){qn=mlcm(qn,x.d);};for(auto&row:Az)for(auto&x:row)include(x);for(auto&[l,row]:raw)for(auto&x:row)include(x);for(auto&[l,row]:dv)for(auto&x:row)include(x);
 MP qz;nmod_mpoly_derivative(qz.v,qn.v,1,mc);qn=mquot(qn,mgcd(qn,qz));RF q(qn);SV qc=zcoeff(q);if(qc.size()<2)throw std::runtime_error("no z poles");q=q/qc.back();qc=zcoeff(q);int s=qc.size()-1;
 std::vector<SV>B(r*r);for(int i=0;i<r;++i)for(int j=0;j<r;++j)B[i*r+j]=zcoeff(q*Az[i][j]);
 int K=1;RF D;std::map<std::string,std::vector<SV>>nums,dnums;for(;K<=12;++K){D=rpow(q,K);bool ok=true;for(auto&[name,row]:raw)for(auto&x:row)ok&=zdegree((D*x).d)<=0;for(auto&[name,row]:dv)for(auto&x:row)ok&=zdegree((D*x).d)<=0;if(ok)break;}if(K>12)throw std::runtime_error("could not clear z denominators");
 int N=K*s+limit-1;for(auto&[name,row]:raw)for(auto&x:row){auto a=zcoeff(D*x);N=std::max(N,int(a.size())-1);nums[name].push_back(std::move(a));}for(auto&[name,row]:dv)for(auto&x:row){auto a=zcoeff(D*x);N=std::max(N,int(a.size())-1);dnums[name].push_back(std::move(a));}
 int cols=r*(N+1)+2*r;auto slot=[&](int c,int d){return c*(N+1)+N-d;};auto epi=[&](int e,int c){return r*(N+1)+e*r+c;};std::vector<std::string>modes={cfg["endpoints"]["0"],cfg["endpoints"]["1"]};
 auto clear=[&](const std::vector<SV>&row){SV out(cols);for(int c=0;c<r;++c)for(int d=0;d<int(row[c].size());++d)out[slot(c,d)]=row[c][d];return out;};
 auto endpoint=[&](const RRow&row,int e){SV out(cols);if(modes[e]=="regular")for(int i=0;i<r;++i)out[epi(e,i)]=at_z(row[i],e);return out;};
 double prep=seconds(start);SymbolicSpan quotient(false);int generators=0;std::cerr<<"symbolic q_degree="<<s<<" K="<<K<<" N="<<N<<" columns="<<cols<<"\n";
 for(int i=0;i<r;++i)for(int j=0;j<=N-s+1;++j){SV row(cols);for(int c=0;c<r;++c)for(int d=0;d<int(B[i*r+c].size());++d)row[slot(c,j+d)]+=B[i*r+c][d];
  for(int d=0;d<=s;++d){int power=j+d-1;if(power>=0)row[slot(i,power)]+=RF(j-(K-1)*d)*qc[d];}
  for(int e=0;e<2;++e)if(modes[e]=="regular"){RF value=(j==0||e==1)?RF(1)/rpow(at_z(q,e),K-1):RF(0);row[epi(e,i)]+=e?-value:value;}
  quotient.insert(std::move(row),generators++);if(j%8==0)std::cerr<<"symbolic generator="<<j<<" elapsed="<<seconds(start)<<"s\n";
 }
 std::map<std::string,SV>coords;std::vector<std::string>labels;for(auto&[name,row]:raw){coords[name]=quotient.reduce(clear(nums.at(name)));labels.push_back(name);for(int e=0;e<2;++e){std::string l="Bz"+std::to_string(e)+name.substr(1);coords[l]=quotient.reduce(endpoint(row,e));labels.push_back(l);}}
 std::sort(labels.begin(),labels.end(),[](auto&a,auto&b){bool ba=a[0]=='B',bb=b[0]=='B';if(ba!=bb)return ba;int ca=weight(a),cb=weight(b);return ca!=cb?ca<cb:a<b;});SymbolicSpan basis;std::vector<std::string>masters,nonzero,zero;auto iszero=[](const SV&v){for(auto&x:v)if(!x.zero())return false;return true;};for(auto&l:labels)if(iszero(coords[l]))zero.push_back(l);else{nonzero.push_back(l);if(basis.insert(coords[l],masters.size()))masters.push_back(l);}
 json out={{"method","symbolic F_p(t) z-IBP total-derivative ansatz"},{"prime",prime},{"t0",cfg.at("t0")},{"masters",masters},{"master_count",masters.size()},{"tested_objects",nonzero},{"zero_objects",zero},{"reductions",json::object()},{"de",json::object()},{"Taylor_order",0},{"sample_points",json::array()},{"input",argv[1]},{"max_dot",dot},{"diagnostics",{{"q_degree",s},{"denominator_power",K},{"numerator_degree",N},{"IBP_generators",generators},{"IBP_rank",quotient.rows.size()}}}};U t0=constant(readr(cfg.at("t0")));
 auto format=[&](const SV&v){json row=json::object();for(int i=0;i<int(masters.size());++i)if(!v[i].zero())row[masters[i]]=as_fraction(v[i],t0).value();return row;};
 for(auto&l:nonzero)if(std::find(masters.begin(),masters.end(),l)==masters.end())out["reductions"][l]=format(basis.solve(coords[l],masters.size()));
 for(auto&l:masters){std::string j=parent_j(l);SV derivative=l[0]=='J'?quotient.reduce(clear(dnums.at(j))):quotient.reduce(endpoint(dv.at(j),l[2]-'0'));auto solution=basis.solve(std::move(derivative),masters.size());std::vector<Fraction>fractions;for(auto&x:solution)fractions.push_back(as_fraction(x,t0));out["de"][l]={{"rhs",format(solution)},{"degree",row_degree(fractions)}};}
 out["timings"]={{"preparation_seconds",prep},{"total_compute_seconds",seconds(start)}};out["caveat"]="Exact symbolic-t ansatz identities at the supplied finite-field kinematics; finite ansatz bounds do not prove completeness for arbitrary families.";std::filesystem::path path(argv[2]);if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());{std::ofstream f(path);f<<out.dump(2);}std::cout<<"symbolic masters="<<masters.size()<<" elapsed="<<seconds(start)<<"s\n";return 0;
 }catch(const std::exception&e){std::cerr<<"symbolic_j: "<<e.what()<<"\n";return 1;}}
