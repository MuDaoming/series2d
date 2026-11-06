// Validate every coefficient, including implicit zeros, against fresh J-IBP.
#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"

int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("verify_connection CONNECTION POINT_MANIFEST OUTPUT");
 json answer,manifest;{std::ifstream f(argv[1]);f>>answer;}{std::ifstream f(argv[2]);f>>manifest;}
 if(answer.value("status",std::string("complete"))!="complete")throw std::runtime_error("connection is incomplete");
 prime=answer.at("prime");nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);
 auto masters=answer.at("masters").get<std::vector<std::string>>();int n=masters.size();
 std::map<std::string,int>index;for(int i=0;i<n;++i)index[masters[i]]=i;
 std::map<std::pair<int,int>,Fraction> functions;
 for(auto row=answer.at("de").begin();row!=answer.at("de").end();++row)
  for(auto col=row.value().at("rhs").begin();col!=row.value().at("rhs").end();++col){
   Fraction f;f.n=readp(col.value().at("numerator_tau"));f.d=readp(col.value().at("denominator_tau"));
   functions[{index.at(row.key()),index.at(col.key())}]=std::move(f);
  }
 U t0=constant(readr(answer.at("t0")));long long checks=0;json points=json::array();
 for(auto& path:manifest.at("points")){
  json point;{std::ifstream f(path.get<std::string>());f>>point;}
  if(point.at("prime").get<U>()!=prime||point.at("masters_at_point").get<std::vector<std::string>>()!=masters||!point.at("unclosed_derivatives").empty())throw std::runtime_error("point field/basis/closure mismatch");
  if(point.at("DE_values").size()!=masters.size())throw std::runtime_error("point DE coverage incomplete");
  U x=sub(constant(readr(point.at("t_point"))),t0);
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   U got=0;auto it=functions.find({i,j});
   if(it!=functions.end()){U d=nmod_poly_evaluate_nmod(it->second.d.v,x);if(!d)throw std::runtime_error("candidate has a pole at verification point");got=nmod_div(nmod_poly_evaluate_nmod(it->second.n.v,x),d,mod);}
   U expected=point.at("DE_values").at(masters[i]).value(masters[j],U(0));++checks;
   if(got!=expected)throw std::runtime_error("coefficient mismatch at "+point.at("t_point").dump()+" row="+masters[i]+" column="+masters[j]);
  }
  points.push_back(point.at("t_point"));
 }
 json result={{"all_pass",true},{"prime",prime},{"master_count",n},{"point_count",points.size()},
              {"scalar_checks_including_zeros",checks},{"points",points}};
 std::ofstream f(argv[3]);f<<result.dump(2)<<"\n";std::cout<<result.dump()<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<"verify_connection: "<<e.what()<<"\n";return 1;}}
