#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
using json=nlohmann::json;
int main(int argc,char**argv){try{
 if(argc!=4&&argc!=5)throw std::runtime_error("usage: compare_blade FI_REDUCTION BLADE_MODP REPORT [--allow-subset]");
 bool subset=argc==5&&std::string(argv[4])=="--allow-subset";if(argc==5&&!subset)throw std::runtime_error("unknown comparison option");
 json fi,blade;{std::ifstream f(argv[1]);f>>fi;}{std::ifstream f(argv[2]);f>>blade;}if(fi.at("prime")!=blade.at("prime"))throw std::runtime_error("prime mismatch");if(!fi.at("all_targets_reduced").get<bool>())throw std::runtime_error("unreduced FI targets");json report={{"prime",fi.at("prime")},{"matched",0},{"mismatches",json::array()}};int matched=0;for(auto&rec:blade.at("records")){std::string label="J{";for(size_t i=0;i<rec["indices"].size();++i)label+=(i?",":"")+rec["indices"][i].dump();label+="}";if(rec.at("unresolved").get<bool>())throw std::runtime_error("Blade unresolved original master: "+label);auto row=fi.at("reductions").at(label).at("fi_coefficients");json expected=rec.contains("coefficients")?rec.at("coefficients"):json{{"FI111",rec.at("sunset")},{"FI110",rec.at("tadpole_product")}};bool equal=true;for(auto it=row.begin();it!=row.end();++it)if(it.value()!=expected.value(it.key(),json(0)))equal=false;for(auto it=expected.begin();it!=expected.end();++it)if(it.value()!=row.value(it.key(),json(0)))equal=false;if(!equal)report["mismatches"].push_back({{"label",label},{"new",row},{"blade",rec}});else ++matched;}
 if(!subset&&size_t(matched)+report["mismatches"].size()!=fi.at("reductions").size())throw std::runtime_error("target coverage mismatch");
 report["full_target_coverage"]=size_t(matched)+report["mismatches"].size()==fi.at("reductions").size();report["fi_target_count"]=fi.at("reductions").size();
 report["matched"]=matched;report["all_match"]=report["mismatches"].empty();{std::ofstream f(argv[3]);f<<report.dump(2);}std::cout<<report.dump(2)<<"\n";return report["mismatches"].empty()?0:1;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
