// Numerical first-integral connection in MFBI periods with proper rational rows.
// The t derivative includes the derivative of the rational row prefactor.
#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"
#include "../../j_integral_global/native/triangular.hpp"

ZR univar(const RF&x){return ZR(univariate(x.n),univariate(x.d));}
std::string nu_label(const Nu&nu){std::string s="{";for(size_t j=0;j<nu.size();++j)s+=(j?",":"")+std::to_string(nu[j]);return s+"}";}
RF cut_pole_polynomial(const json&cfg,const Nu&nu){
 json local=cfg;local.erase("specialize_t");local.erase("fbi_corner_aliases");
 local["canonical_fbi_symmetries"]=false;local["z0_fbi_analytic"]=false;
 local["propagators"]=json::array();for(int j=0;j<int(nu.size());++j)if(nu[j])local["propagators"].push_back(cfg["propagators"][j]);
 local["cut_sector"]=Nu(local["propagators"].size(),1);local["targets"]=json::array();
 NativeFamily cut(local);RMat At,Az;std::map<std::string,RRow>raw;cut.model(At,Az,raw);
 if(cut.rank!=1)throw std::runtime_error("non-scalar canonical cut FBI");
 MP den=Az[0][0].d,deriv;nmod_mpoly_derivative(deriv.v,den.v,1,mc);
 return RF(mquot(den,mgcd(den,deriv)));
}

int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("period_point FAMILY POINT OUTPUT");
 json cfg;{std::ifstream f(argv[1]);f>>cfg;}
 prime=cfg.at("prime");nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);
 cfg["specialize_t"]=argv[2];cfg["targets"]=json::array();auto start=Clock::now();
 NativeFamily family(cfg);RMat At,Az;std::map<std::string,RRow> unused;family.model(At,Az,unused);
 int r=family.rank;std::vector<std::vector<ZR>>at(r),az(r);
 for(int i=0;i<r;++i){for(auto&x:At[i])at[i].push_back(univar(x));for(auto&x:Az[i])az[i].push_back(univar(x));}
 double prep=seconds(start);std::cerr<<"FBI="<<r<<" model="<<prep<<"s\n";
 TriangularIntegrator integrator(az,24);integrator.atomic_mode=cfg.value("atomic_period",false);integrator.complete(24);
 if(cfg.value("z0_fbi_analytic",false)){
  std::vector<int> orders;for(auto&nu:family.masters){int n=-1;for(int j=0;j<family.n;++j)if(family.branch[j]==2)n+=nu[j];orders.push_back(n);}
  integrator.impose_analytic_boundary(orders);
  for(auto&[i,j]:family.endpoint_corner_aliases()){
   C relation=integrator.endpoint.finite_part(i,ZR(1));axpy(relation,integrator.endpoint.finite_part(j,ZR(1)),prime-1);integrator.completion.insert(std::move(relation));
  }
 }
 struct Master{std::string name;int row,power;bool boundary;ZR weight,dweight;};
 std::vector<Master> masters;Span span;int candidate_count=0;
 auto offer=[&](std::string name,int i,int k,bool boundary,const ZR&weight,const ZR&dweight,C c){++candidate_count;
  if(span.insert(std::move(c),masters.size()))masters.push_back({name,i,k,boundary,weight,dweight});};
 for(int i=r-1;i>=0;--i){std::vector<ZR>v(r);v[i]=ZR(1);offer("Bz0"+nu_label(family.masters[i]),i,0,true,ZR(1),ZR(0),integrator.boundary(v,0));}
 int boundaries=masters.size();
 for(int i=r-1;i>=0;--i){
  RF q=cut_pole_polynomial(cfg,family.masters[i]);
  int bound=nmod_mpoly_degree_si(q.n.v,1,mc);
  U point=constant(readr(std::string(argv[2])));
  auto qs=specialize(q,point);P rem;nmod_poly_rem(rem.v,qs.first.v,az[i][i].d.v);
  std::cerr<<"pole row="<<i<<" global="<<az[i][i].d.deg()<<" cut="<<qs.first.deg()<<" coverage="<<rem.zero()<<"\n";
  for(int k=-1;k<bound;++k){RF weight(1);if(k>=0)weight=rpow(Expression("z").parse(),k)/q;
   auto w=specialize(weight,point),dw=specialize(weight.diff(0),point);
   ZR v0(w.first,w.second),v1(dw.first,dw.second);std::vector<ZR>v(r);v[i]=v0;
   offer("P"+nu_label(family.masters[i])+":"+std::to_string(k),i,k,false,v0,v1,integrator.integrate(std::move(v)));}
  if(i%10==0)std::cerr<<"basis row="<<i<<" dimension="<<masters.size()<<" elapsed="<<seconds(start)<<"s\n";
 }
 double basis_time=seconds(start)-prep;
 json out={{"prime",prime},{"t_point",argv[2]},{"FBI_rank",r},{"FBI_masters",family.masters},
           {"masters_at_point",json::array()},{"DE_values",json::object()},{"unclosed_derivatives",json::array()},
           {"boundary_rank",boundaries},{"candidate_count",candidate_count},{"master_descriptors",json::array()}};
 for(auto&m:masters){out["masters_at_point"].push_back(m.name);out["master_descriptors"].push_back({{"name",m.name},{"FBI_row",m.row},{"power",m.power},{"boundary",m.boundary}});}
 for(int a=0;a<int(masters.size());++a){auto&m=masters[a];std::vector<ZR>v(r);
  for(int j=m.row;j<r;++j)v[j]=at[m.row][j]*m.weight;
  v[m.row]=v[m.row]+m.dweight;
  C d=m.boundary?integrator.boundary(v,0):integrator.integrate(std::move(v));
  try{auto row=span.solve(std::move(d),masters.size());json values=json::object();for(int j=0;j<int(masters.size());++j)if(row[j])values[masters[j].name]=row[j];out["DE_values"][m.name]=std::move(values);}
  catch(const std::exception&){out["unclosed_derivatives"].push_back(m.name);}
  if(a%50==0)std::cerr<<"derivative="<<a<<" elapsed="<<seconds(start)<<"s\n";
 }
 out["master_count_at_point"]=masters.size();out["timings"]={{"FBI_model_seconds",prep},{"basis_seconds",basis_time},{"total_seconds",seconds(start)}};
 out["atomic_cache"]={{"hits",integrator.atomic_hits},{"misses",integrator.atomic_misses}};
 std::filesystem::path output(argv[3]);std::filesystem::create_directories(output.parent_path());std::ofstream f(output);f<<out.dump(2)<<"\n";
 std::cerr<<"done dimension="<<masters.size()<<" unclosed="<<out["unclosed_derivatives"].size()<<" total="<<seconds(start)<<"s\n";
 return out["unclosed_derivatives"].empty()?0:2;
}catch(const std::exception&e){std::cerr<<"period_point: "<<e.what()<<"\n";return 1;}}
