// Compute only the homogeneous J blocks, using scalar maximal-cut FBI systems.
// No global J connection or global J reduction is read by this worker.
#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"
#include "../../j_integral_global/native/triangular.hpp"

Nu indices(const std::string& name) {
    auto a=name.find('{'), b=name.find('}');
    if(a==std::string::npos || b==std::string::npos) throw std::runtime_error("invalid master label");
    return json::parse("["+name.substr(a+1,b-a-1)+"]").get<Nu>();
}
std::string jlabel(const Nu& nu) {
    std::string s="J{";
    for(size_t i=0;i<nu.size();++i) s+=(i?",":"")+std::to_string(nu[i]);
    return s+"}";
}

int main(int argc,char**argv) {
 try {
    if(argc!=5) throw std::runtime_error("cut_diagonal FAMILY BASIS POINT OUTPUT");
    json cfg,basis;
    {std::ifstream f(argv[1]);f>>cfg;}
    {std::ifstream f(argv[2]);f>>basis;}
    prime=cfg.at("prime");nmod_init(&mod,prime);
    nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);
    U point=constant(readr(std::string(argv[3])));
    std::map<Nu,std::vector<std::string>> groups;
    for(auto& value:basis) {
        std::string name=value;
        if(name[0]!='J')continue;
        Nu sector=indices(name);for(int& v:sector)v=v>0;
        groups[sector].push_back(name);
    }
    json out={{"prime",prime},{"t_point",argv[3]},{"blocks",json::array()},{"DE_values",json::object()}};
    auto start=Clock::now();
    for(auto&[sector,names]:groups) {
      json report={{"sector",sector},{"masters",names}};
      auto stamp=Clock::now();
      try {
        std::vector<int> active;
        json local=cfg;local["propagators"]=json::array();
        for(int i=0;i<int(sector.size());++i)if(sector[i]) {
            active.push_back(i);local["propagators"].push_back(cfg["propagators"][i]);
        }
        local["canonical_fbi_symmetries"]=false;local["z0_fbi_analytic"]=false;
        local.erase("fbi_corner_aliases");local.erase("specialize_t");
        local["cut_sector"]=Nu(active.size(),1);local["targets"]=json::array();
        std::map<std::string,std::string> trimmed;
        for(auto& name:names) {
            Nu nu=indices(name),v;for(int i:active)v.push_back(nu[i]);
            local["targets"].push_back(v);trimmed[name]=jlabel(v);
        }
        NativeFamily family(local);RMat At,Az;std::map<std::string,RRow> raw,drows;
        family.model(At,Az,raw);int r=family.rank;
        if(r!=1)throw std::runtime_error("maximal-cut FBI is not scalar");
        for(auto&[name,row]:raw) {
            RRow d(r);
            for(int i=0;i<r;++i) {
                d[i]=row[i].diff(0);
                for(int j=0;j<r;++j)d[i]+=row[j]*At[j][i];
            }
            drows[name]=std::move(d);
        }
        P q(1);auto include=[&](const std::pair<P,P>& x) {
            P g;nmod_poly_gcd(g.v,q.v,x.second.v);q=quotient(q,g)*x.second;
        };
        std::map<std::string,RationalRow> values,derivs;
        std::vector<RationalRow> all;
        for(auto&[name,row]:raw) {
            for(auto& x:row){auto s=specialize(x,point);include(s);values[name].push_back(s);}
            for(auto& x:drows.at(name)){auto s=specialize(x,point);include(s);derivs[name].push_back(s);}
            all.push_back(values[name]);all.push_back(derivs[name]);
        }
        auto az=specialize(Az[0][0],point);include(az);
        P g;nmod_poly_gcd(g.v,q.v,derivative(q).v);q=quotient(q,g);
        nmod_poly_make_monic(q.v,q.v);
        TriangularIntegrator integrator({{ZR(az.first,az.second)}},24);
        integrator.complete(24);
        // Homogeneous block: quotient every cut endpoint mode, including modes
        // represented by a derivative germ rather than a nonzero endpoint value.
        integrator.completion.insert(C{{"endpoint:0:0",1}});
        auto integral=[&](const RationalRow& row) {
            std::vector<ZR> v;for(auto& x:row)v.emplace_back(x.first,x.second);
            return integrator.integrate(std::move(v));
        };
        Span span;int boundary_rank=0;
        std::map<std::string,C> coords;
        for(int i=0;i<int(names.size());++i) {
            auto& name=names[i];coords[name]=integral(values.at(trimmed.at(name)));
            if(!span.insert(coords[name],boundary_rank+i))
                throw std::runtime_error("global J basis is dependent modulo cut boundaries");
        }
        for(auto& name:names) {
            auto row=span.solve(integral(derivs.at(trimmed.at(name))),boundary_rank+names.size());
            json coeff=json::object();
            for(int j=0;j<int(names.size());++j)if(row[boundary_rank+j])coeff[names[j]]=row[boundary_rank+j];
            out["DE_values"][name]=std::move(coeff);
        }
        report["status"]="ok";report["FBI_rank"]=r;report["boundary_rank"]=boundary_rank;
      } catch(const std::exception& e) {
        report["status"]="failed";report["error"]=e.what();
      }
      report["seconds"]=seconds(stamp);out["blocks"].push_back(report);
      std::cerr<<out["blocks"].size()<<"/"<<groups.size()<<" "<<report["status"]<<" "<<seconds(start)<<"s\n";
    }
    out["total_seconds"]=seconds(start);
    std::filesystem::path path(argv[4]);std::filesystem::create_directories(path.parent_path());
    std::ofstream f(path);f<<out.dump(2)<<"\n";
    return 0;
 }catch(const std::exception& e){std::cerr<<"cut_diagonal: "<<e.what()<<"\n";return 1;}
}
