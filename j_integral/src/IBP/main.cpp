// Public, input-driven fixed-t J-IBP entry. No AB or archived code dependency.
#include "kernel/point.cpp"
#include <flint/ulong_extras.h>
#include <regex>
#include <unistd.h>

namespace fs=std::filesystem;
json load_json(const fs::path&path){std::ifstream f(path);if(!f)throw std::runtime_error("cannot read "+path.string());json j;f>>j;return j;}
void save_json(const fs::path&path,const json&j){std::ofstream f(path);if(!f)throw std::runtime_error("cannot write "+path.string());f<<j.dump(2)<<'\n';}
struct Scratch {
 fs::path path;
 Scratch(){path=fs::temp_directory_path()/("j-ibp-"+std::to_string(getpid())+"-"+std::to_string(Clock::now().time_since_epoch().count()));if(!fs::create_directory(path))throw std::runtime_error("cannot create temporary directory");}
 ~Scratch(){std::error_code ec;fs::remove_all(path,ec);}
};
std::string scalar_text(const json&value){if(!value.is_string()&&!value.is_number_integer())throw std::runtime_error("exact scalars must be integers or strings; floating point is unsupported");return value.is_string()?value.get<std::string>():value.dump();}
std::string specialize_text(const json&value,const json&parameters){
 std::string text=scalar_text(value),out;
 for(size_t i=0;i<text.size();){unsigned char c=text[i];if(std::isalpha(c)||c=='_'){
   size_t start=i++;while(i<text.size()&&(std::isalnum(static_cast<unsigned char>(text[i]))||text[i]=='_'))++i;
   std::string name=text.substr(start,i-start);if(!parameters.contains(name))throw std::runtime_error("unassigned scalar parameter: "+name);
   std::string v=scalar_text(parameters.at(name));if(!std::regex_match(v,std::regex(R"([0-9+\-*/^() \t]+)")))throw std::runtime_error("parameter assignments must be exact scalar expressions");out+='(';out+=v;out+=')';
 }else{if(!std::isdigit(c)&&!std::isspace(c)&&std::string("+-*/^()").find(c)==std::string::npos)throw std::runtime_error("unsupported scalar expression token");out+=text[i++];}}
 if(out.empty())throw std::runtime_error("empty scalar expression");return out;
}
json array_or_file(const json&value,const fs::path&base){json a=value.is_string()?load_json(base/value.get<std::string>()):value;if(!a.is_array())throw std::runtime_error("expected array or a relative JSON array filename");return a;}
Nu parse_label(const std::string&name,int n){std::smatch m;if(!std::regex_match(name,m,std::regex(R"((J|Bz0|Bz1)\{([0-9,]+)\})")))throw std::runtime_error("invalid integral label: "+name);Nu nu=json::parse("["+m[2].str()+"]").get<Nu>();if(int(nu.size())!=n)throw std::runtime_error("basis/derivative index length differs from propagators");return nu;}
int main(int argc,char**argv){try{
 std::string input,output;
 for(int i=1;i<argc;++i){std::string arg=argv[i];if(arg=="--help"){std::cout<<"j-ibp --input REQUEST.json --output RESULT.json|-\n";return 0;}if((arg!="--input"&&arg!="--output")||i+1>=argc)throw std::runtime_error("usage: j-ibp --input REQUEST.json --output RESULT.json|-");(arg=="--input"?input:output)=argv[++i];}
 if(input.empty()||output.empty())throw std::runtime_error("both --input and --output are required");
 auto started=Clock::now();fs::path base=fs::absolute(input).parent_path();json request=load_json(input);
 if(request.at("schema_version")!=1)throw std::runtime_error("unsupported schema_version");
 const auto&f=request.at("family");const auto&point=request.at("point");const auto&j=request.at("j");const auto&task=request.at("request");json options=request.value("options",json::object());
 if(f.at("loops").size()!=2)throw std::runtime_error("this engine requires two loop momenta");
 if(j.value("coordinates",std::string("tz"))!="tz")throw std::runtime_error("only X1=(1-z)t,X2=(1-z)(1-t),X3=z is supported");
 U p=point.at("prime").get<U>();if(p<3||!n_is_prime(p))throw std::runtime_error("point.prime must be an odd prime");
 json parameters=point.value("parameters",json::object());if(!parameters.is_object())throw std::runtime_error("point.parameters must be an object");
 std::string t=specialize_text(point.at("t"),parameters),dimension=specialize_text(point.at("dimension"),parameters);
 json cfg={{"kind","two_loop_family"},{"prime",p},{"dimension",dimension},{"gram",json::array()},{"propagators",f.at("propagators")}};
 int E=f.at("external_momenta").size(),n=cfg["propagators"].size();if(n<3||n>20)throw std::runtime_error("requires 3..20 propagators");
 if(int(f.at("gram").size())!=E)throw std::runtime_error("Gram matrix length differs from external_momenta");
 for(int a=0;a<E;++a){if(int(f["gram"][a].size())!=E)throw std::runtime_error("Gram matrix must be square");json row=json::array();for(int b=0;b<E;++b)row.push_back(specialize_text(f["gram"][a][b],parameters));cfg["gram"].push_back(row);}
 // The arithmetic kernel checks Gram symmetry after exact specialization.
 json branches=j.value("branches",json());if(!branches.is_null()&&(!branches.is_array()||int(branches.size())!=n))throw std::runtime_error("j.branches must have one entry per propagator");
 std::set<int>seen;
 for(int i=0;i<n;++i){auto&prop=cfg["propagators"][i];if(prop.at("loop").size()!=2||int(prop.at("shift").size())!=E)throw std::runtime_error("propagator loop/shift length differs");
  int branch;
  if(!branches.is_null())branch=branches[i].get<int>();
  else{auto loop=prop["loop"].get<std::vector<int>>();if((loop[0]==1||loop[0]==-1)&&loop[1]==0)branch=0;else if(loop[0]==0&&(loop[1]==1||loop[1]==-1))branch=1;else if((loop[0]==1||loop[0]==-1)&&(loop[1]==1||loop[1]==-1))branch=2;else throw std::runtime_error("noncanonical loop routing requires explicit j.branches");}
  if(branch<0||branch>2)throw std::runtime_error("branch labels must be 0,1,2");seen.insert(branch);prop["branch"]=branch;
  for(auto&v:prop["loop"])v=specialize_text(v,parameters);for(auto&v:prop["shift"])v=specialize_text(v,parameters);prop["mass2"]=specialize_text(prop.at("mass2"),parameters);
 }
 if(seen!=std::set<int>{0,1,2})throw std::runtime_error("ordinary J requires all three branches");
 std::string z0=j.at("boundary").at("z0"),z1=j.at("boundary").at("z1");
 if(z0!="regularized"&&z0!="fbi_analytic")throw std::runtime_error("z0 must be regularized or fbi_analytic");if(z1!="regulated_zero")throw std::runtime_error("only the declared regulated_zero boundary at z1 is supported");
 cfg["canonical_fbi_symmetries"]=true;cfg["include_fbi_boundaries"]=true;cfg["z0_fbi_analytic"]=z0=="fbi_analytic";cfg["endpoints"]={{"0","regular"},{"1","regulated_zero"}};
 if(j.contains("cut_sector")){if(z0=="fbi_analytic")throw std::runtime_error("fbi_analytic boundary cannot be imposed on a propagator cut");cfg["canonical_fbi_symmetries"]=false;cfg["cut_sector"]=j["cut_sector"];if(int(cfg["cut_sector"].size())!=n)throw std::runtime_error("cut_sector length differs");for(auto&v:cfg["cut_sector"])if(!v.is_number_integer()||(v!=0&&v!=1))throw std::runtime_error("cut_sector must contain only0/1");}
 auto valid_nu=[&](const Nu&nu){if(int(nu.size())!=n)throw std::runtime_error("target index length differs");std::set<int>bs;for(int i=0;i<n;++i){if(nu[i]<0)throw std::runtime_error("rank/numerator indices are unsupported");if(nu[i])bs.insert(cfg["propagators"][i]["branch"].get<int>());if(cfg.contains("cut_sector")&&cfg["cut_sector"][i]==1&&nu[i]==0)throw std::runtime_error("target is outside cut_sector");}if(bs!=std::set<int>{0,1,2})throw std::runtime_error("target must retain every branch");};
 if(task.contains("targets")==task.contains("generate_targets"))throw std::runtime_error("provide exactly one of request.targets and request.generate_targets");
 std::vector<Nu>targets;
 if(task.contains("targets")){for(auto&nu:array_or_file(task["targets"],base)){Nu v=nu.get<Nu>();valid_nu(v);targets.push_back(v);}}
 else{auto g=task["generate_targets"];int dot=g.at("max_dot");std::string scope=g.value("scope",std::string("all"));if(dot<0||!(scope=="all"||scope=="top"))throw std::runtime_error("invalid target generation");
  for(int mask=1;mask<(1<<n);++mask){if(scope=="top"&&mask!=(1<<n)-1)continue;Nu v(n);for(int i=0;i<n;++i)v[i]=(mask>>i)&1;std::set<int>bs;bool retained=true;for(int i=0;i<n;++i){if(v[i])bs.insert(cfg["propagators"][i]["branch"].get<int>());if(cfg.contains("cut_sector")&&cfg["cut_sector"][i]==1&&!v[i])retained=false;}if(!retained||bs.size()!=3)continue;
   std::function<void(int,int)>visit=[&](int i,int left){if(i==n){targets.push_back(v);return;}if(!((mask>>i)&1)){visit(i+1,left);return;}for(int k=0;k<=left;++k){v[i]=1+k;visit(i+1,left-k);}v[i]=1;};visit(0,dot);
  }
 }
 if(targets.empty())throw std::runtime_error("no requested targets");std::set<Nu>unique_targets(targets.begin(),targets.end());if(unique_targets.size()!=targets.size())throw std::runtime_error("duplicate targets");
 json basis;if(task.contains("basis")){basis=array_or_file(task["basis"],base);if(basis.empty())throw std::runtime_error("provided basis cannot be empty");std::set<std::string>names;for(auto&v:basis){std::string name=v;valid_nu(parse_label(name,n));if(!names.insert(name).second)throw std::runtime_error("duplicate basis label");}}
 json derivatives=task.value("derivatives",json("basis"));if(derivatives.is_string()&&derivatives!="basis")derivatives=array_or_file(derivatives,base);if(!derivatives.is_array()&&derivatives!="basis")throw std::runtime_error("derivatives must be basis or an array");if(derivatives.is_array())for(auto&name:derivatives)valid_nu(parse_label(name.get<std::string>(),n));
 // Native fixed-basis coordinates and derivatives need the corresponding J parents.
 std::set<Nu>parents=unique_targets;if(!basis.is_null())for(auto&v:basis)parents.insert(parse_label(v.get<std::string>(),n));json native_targets=json::array();for(auto&v:parents)native_targets.push_back(v);
 int completion=options.value("completion",16),germ=options.value("germ_order",16);if(completion<1||germ<1)throw std::runtime_error("completion and germ_order must be positive");
 Scratch scratch;save_json(scratch.path/"family.json",cfg);save_json(scratch.path/"targets.json",native_targets);
 std::vector<std::string>args={"j-ibp-kernel",(scratch.path/"family.json").string(),(scratch.path/"result.json").string(),"--point",t,"--targets-file",(scratch.path/"targets.json").string(),"--integration","atomic","--values","1","--compact","1","--completion",std::to_string(completion),"--germ-order",std::to_string(germ)};
 if(!basis.is_null()){save_json(scratch.path/"basis.json",basis);args.insert(args.end(),{"--basis-file",(scratch.path/"basis.json").string()});}
 if(derivatives.is_array()){save_json(scratch.path/"derivatives.json",derivatives);args.insert(args.end(),{"--derivative-names-file",(scratch.path/"derivatives.json").string()});}
 std::vector<char*>argv_native;for(auto&s:args)argv_native.push_back(s.data());
 auto*original=std::cout.rdbuf(std::cerr.rdbuf());int code=j_ibp_kernel_main(argv_native.size(),argv_native.data());std::cout.rdbuf(original);if(code)throw std::runtime_error("J-IBP kernel failed; see preceding diagnostic");
 json native=load_json(scratch.path/"result.json"),reductions=json::object();
 for(auto&nu:targets)reductions[label(nu)]=native.at("reduction_values").at(label(nu));
 if(task.value("include_boundary_reductions",true))for(auto it=native["reduction_values"].begin();it!=native["reduction_values"].end();++it)if(it.key().rfind("Bz",0)==0)reductions[it.key()]=it.value();
 bool closed=native["unclosed_derivatives"].empty();json result={{"schema_version",1},{"status",closed?"complete":"unclosed"},{"point",{{"prime",p},{"dimension",dimension},{"t",t},{"parameters",parameters}}},{"definition",{{"coordinates","tz"},{"branches",json::array()},{"boundary",j["boundary"]}}},{"basis",native["masters_at_point"]},{"basis_status",basis.is_null()?"candidate_at_this_point":"provided_basis_independent_and_spans_requested_targets"},{"reductions",reductions},{"derivatives",native["DE_values"]},{"unclosed_derivatives",native["unclosed_derivatives"]},{"diagnostics",{{"FBI_rank",native["FBI_rank"]},{"basis_size",native["master_count_at_point"]},{"requested_J_targets",targets.size()},{"auxiliary_J_parents",parents.size()-targets.size()},{"completion_relations",native["completion_relations"]},{"physical_boundary_relations",native["physical_boundary_relations"]},{"scalar_certificate_checks",native["scalar_certificate_checks"]},{"timings",native["timings"]},{"wall_seconds",seconds(started)}}}};
 for(auto&prop:cfg["propagators"])result["definition"]["branches"].push_back(prop["branch"]);if(j.contains("cut_sector"))result["definition"]["cut_sector"]=j["cut_sector"];
 if(output=="-")std::cout<<result.dump()<<'\n';else{fs::path path=fs::absolute(output);fs::create_directories(path.parent_path());auto temporary=path;temporary+=".tmp-"+std::to_string(getpid());save_json(temporary,result);fs::rename(temporary,path);}
 return closed?0:2;
}catch(const std::exception&e){std::cerr<<"j-ibp: "<<e.what()<<'\n';return 1;}}
