#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <set>
using json=nlohmann::json;
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("usage: check_series REFERENCE_SERIES NATIVE_OUTPUT");
 json a,b;{std::ifstream f(argv[1]);f>>a;}{std::ifstream f(argv[2]);f>>b;}
 if(b.contains("prime")&&b["prime"]!=a["prime"])throw std::runtime_error("prime mismatch");
 int count=0;
 for(auto it=a.at("targets").begin();it!=a.at("targets").end();++it){const auto& actual=b.at("coordinates").at(it.key());size_t limit=std::min(size_t(a.at("order").get<int>()+1),actual.size());
  for(size_t n=0;n<limit;++n){json expected=json::object();for(size_t i=0;i<a.at("period_keys").size();++i){auto& key=a["period_keys"][i];std::string k=key[0].get<std::string>();for(size_t j=1;j<key.size();++j)k+=":"+key[j].dump();if(it.value()[i][n]!=0)expected[k]=it.value()[i][n];}if(expected!=actual[n])throw std::runtime_error("coordinate mismatch: "+it.key()+" order "+std::to_string(n));++count;}
 }
 if(b.at("coordinates").size()!=a.at("targets").size())throw std::runtime_error("target count mismatch");
 std::cout<<"matched "<<count<<" exact finite-field coefficient rows\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
