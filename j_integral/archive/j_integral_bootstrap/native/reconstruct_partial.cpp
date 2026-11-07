// Reconstruct all currently determined columns, retaining explicit failures.
#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"

bool reconstruct_asymmetric(const Interpolator& interp,const std::vector<U>&values,
                            const std::vector<U>&xs,Fraction&out){
 P a;for(size_t i=0;i<interp.basis.size();++i)a+=scale(interp.basis[i],values[i]);
 if(a.zero()){
  if(std::any_of(values.begin(),values.end(),[](U v){return v!=0;}))return false;
  out=Fraction();return true;
 }
 P r0=interp.modulus,r1=a,s0,s1(1);std::vector<Fraction> candidates;
 while(!r1.zero()){
  if(!s1.zero()&&r1.deg()+s1.deg()<int(interp.basis.size())){
   Fraction f;f.n=r1;f.d=s1;f.normal();candidates.push_back(std::move(f));
  }
  P q,rem;nmod_poly_divrem(q.v,rem.v,r0.v,r1.v);P next=s0-q*s1;
  r0=std::move(r1);r1=std::move(rem);s0=std::move(s1);s1=std::move(next);
 }
 std::sort(candidates.begin(),candidates.end(),[](const Fraction&a,const Fraction&b){
  return std::make_pair(a.n.deg()+a.d.deg(),std::max(a.n.deg(),a.d.deg()))<
         std::make_pair(b.n.deg()+b.d.deg(),std::max(b.n.deg(),b.d.deg()));
 });
 for(auto&f:candidates){bool valid=true;
  // Reject wrong candidates at a held-out point before checking training nodes.
  for(size_t k=0;k<xs.size();++k){size_t i=(k+interp.basis.size())%xs.size();
   U d=nmod_poly_evaluate_nmod(f.d.v,xs[i]);
   if(!d||nmod_poly_evaluate_nmod(f.n.v,xs[i])!=mul(values[i],d)){valid=false;break;}
  }
  if(valid){out=std::move(f);return true;}
 }
 return false;
}

int main(int argc,char**argv){try{
 if(argc<3||argc>4)throw std::runtime_error("reconstruct_partial MANIFEST OUTPUT [SELECTION]");
 json manifest,cfg;{std::ifstream f(argv[1]);f>>manifest;}{std::ifstream f(manifest.at("family").get<std::string>());f>>cfg;}
 prime=cfg.at("prime");nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);
 std::string section=manifest.value("section",std::string("DE_values"));
 if(section!="DE_values"&&section!="reduction_values")throw std::runtime_error("unknown section");
 U t0=constant(readr(cfg.at("t0")));auto start=Clock::now();
 using Key=std::pair<std::string,std::string>;std::set<Key> selection;
 bool selected=argc==4;
 std::set<std::string> row_selection;if(manifest.contains("rows"))for(auto&row:manifest.at("rows"))row_selection.insert(row.get<std::string>());
 if(selected){json s;{std::ifstream f(argv[3]);f>>s;}for(auto&v:s.at("pairs"))selection.emplace(v[0],v[1]);}
 std::map<Key,std::vector<U>>columns;std::vector<U>xs;std::set<std::string>all_rows;
 for(auto&path:manifest.at("points")){
  json point;{std::ifstream f(path.get<std::string>());f>>point;}
  if(point.at("prime").get<U>()!=prime)throw std::runtime_error("field differs");
  if(point.contains("unclosed_derivatives")&&!point.at("unclosed_derivatives").empty())throw std::runtime_error("unclosed point");
  U x=sub(constant(readr(point.at("t_point"))),t0);if(std::find(xs.begin(),xs.end(),x)!=xs.end())throw std::runtime_error("duplicate point");
  int k=xs.size();xs.push_back(x);std::set<std::string> present;
  for(auto row=point.at(section).begin();row!=point.at(section).end();++row)
   if(row_selection.empty()||row_selection.count(row.key()))present.insert(row.key());
  if(!row_selection.empty()&&present!=row_selection)throw std::runtime_error("requested DE row missing");
  if(k&&present!=all_rows)throw std::runtime_error("DE row coverage changes between points");
  for(auto row=point.at(section).begin();row!=point.at(section).end();++row){
   if(!row_selection.empty()&&!row_selection.count(row.key()))continue;
   all_rows.insert(row.key());
   for(auto col=row.value().begin();col!=row.value().end();++col){
    Key key={row.key(),col.key()};if(selected&&!selection.count(key))continue;
    auto&v=columns[key];v.resize(k+1);v[k]=col.value().get<U>();
   }
  }
 }
 int holdout=8,train=int(xs.size())-holdout;
 if(train<1)throw std::runtime_error("need >8 points");
 Interpolator interpolator(std::vector<U>(xs.begin(),xs.begin()+train));
 json out={{"prime",prime},{"t0",cfg.at("t0")},{"training_points",train},{"heldout_points",holdout},
           {"de",json::object()},{"reductions",json::object()},{"section",section},
           {"method","Euclidean rational reconstruction allowing asymmetric degrees"},
           {"unreconstructed",json::array()},{"reconstructed_count",0}};
 int count=0,passed=0,maxdegree=0;
 for(auto&[key,values]:columns){values.resize(xs.size());Fraction fit;
  if(reconstruct_asymmetric(interpolator,values,xs,fit)){
   ++passed;maxdegree=std::max({maxdegree,int(fit.n.deg()),int(fit.d.deg())});
   if(section=="DE_values")out["de"][key.first]["rhs"][key.second]=fit.value();
   else out["reductions"][key.first][key.second]=fit.value();
  }else out["unreconstructed"].push_back({key.first,key.second});
  if(++count%2000==0)std::cerr<<"columns="<<count<<" reconstructed="<<passed<<" elapsed="<<seconds(start)<<"s\n";
 }
 out["reconstructed_count"]=passed;out["column_count"]=columns.size();out["max_coefficient_degree"]=maxdegree;
 out["rows_present"]=all_rows;out["total_seconds"]=seconds(start);
 std::ofstream f(argv[2]);f<<out.dump(2)<<"\n";
 std::cout<<"columns="<<columns.size()<<" reconstructed="<<passed<<" missing="<<out["unreconstructed"].size()<<" seconds="<<seconds(start)<<"\n";
 return out["unreconstructed"].empty()?0:2;
}catch(const std::exception&e){std::cerr<<"reconstruct_partial: "<<e.what()<<"\n";return 1;}}
