"""Generate an instrumented COPY of the production kernel; no algorithm edits."""
from pathlib import Path
import shutil,subprocess
here=Path(__file__).resolve().parents[1]
production=here.parents[1]/'src/IBP/kernel'
out=here/'build/seed_audit_kernel';out.mkdir(parents=True,exist_ok=True)
for src in production.iterdir():
    if src.is_file():shutil.copyfile(src,out/src.name)
text=(out/'z_integrator.hpp').read_text()
log=r'''
std::string seed_audit_phase;
std::map<std::string,std::map<std::string,json>> seed_audit_entries;
void seed_audit_add(int row,const std::string&type,const P&f,int power,int degree){
 json key={{"component",row},{"kind",type},{"power",power},{"numerator_degree",degree}};if(type=="finite")key["factor"]=polyjson(f);
 auto name=key.dump();auto&entry=seed_audit_entries[seed_audit_phase][name];if(entry.is_null()){entry=key;entry["calls"]=0;}entry["calls"]=entry["calls"].get<long long>()+1;
}
void seed_audit_finite(int row,const P&f,int power,const P&a){for(int degree=0;degree<=a.deg();++degree)if(a.at(degree))seed_audit_add(row,"finite",f,power,degree);}
void seed_audit_polynomial(int row,int degree){seed_audit_add(row,"polynomial",P(1),0,degree);}
void seed_audit_primitive(int row,const ZR&g){if(g.d.deg()==0){for(int d=0;d<=g.n.deg();++d)if(g.n.at(d))seed_audit_polynomial(row,d);}else{Factors factors;auto&fs=factors.get(g.d);if(fs.size()!=1)throw std::runtime_error("unexpected multi-factor completion primitive");seed_audit_finite(row,fs[0].first,fs[0].second,g.n);}}
'''
def replace(old,new):
    global text
    assert text.count(old)==1,(old,text.count(old))
    text=text.replace(old,new)
replace('struct ScalarHermite {',log+'\nstruct ScalarHermite {\n int audit_row;')
replace('ScalarHermite(ZR A,Factors&f):A(std::move(A)),factors(&f){}','ScalarHermite(ZR A,Factors&f,int row):audit_row(row),A(std::move(A)),factors(&f){}')
replace('std::make_unique<ScalarHermite>(this->A[i][i],factors)','std::make_unique<ScalarHermite>(this->A[i][i],factors,i)')
replace('ZR g(a,ppow(f,k-1));v=', 'ZR g(a,ppow(f,k-1));seed_audit_finite(audit_row,f,k-1,a);v=')
replace('ZR g=zp(p);v=', 'ZR g=zp(p);seed_audit_polynomial(audit_row,d+1);v=')
replace('ZR g(nmod_div(c,lambda,mod));v=', 'ZR g(nmod_div(c,lambda,mod));seed_audit_polynomial(audit_row,0);v=')
replace('auto submit=[&](const ZR&g){std::vector<ZR>v(r);', 'auto submit=[&](const ZR&g){auto previous_phase=seed_audit_phase;seed_audit_phase="supplementary_outer";seed_audit_primitive(i,g);seed_audit_phase="supplementary_hermite";std::vector<ZR>v(r);')
replace('++completion_count;}};', '++completion_count;}seed_audit_phase=previous_phase;};')
(out/'z_integrator.hpp').write_text(text)
source=(here/'src/experiment.cpp').read_text().replace('#include "../../../src/IBP/kernel/point.cpp"','#include "point.cpp"')
extra=r'''
json seed_summary(const std::map<std::string,json>&entries){int finite=0;std::map<int,int>powers,degrees;long long calls=0;for(auto&[k,e]:entries){calls+=e.at("calls").get<long long>();if(e["kind"]=="finite"){++finite;++powers[e["power"].get<int>()];}else ++degrees[e["numerator_degree"].get<int>()];}json fs=json::object(),ps=json::object();for(auto&[k,v]:powers)fs[std::to_string(k)]=v;for(auto&[k,v]:degrees)ps[std::to_string(k)]=v;return {{"distinct_seeds",entries.size()},{"finite_seeds",finite},{"finite_primitive_power_histogram",fs},{"polynomial_primitive_degree_histogram",ps},{"calls",calls}};}
void audit_seeds(const std::string&input,const std::string&output){auto req=load_json(input);init(req["point"]["prime"]);auto cfg=native_config(req);auto start=Clock::now();NativeFamily family(cfg);RMat At,Az;std::map<std::string,RRow>raw;family.model(At,Az,raw);std::vector<std::vector<ZR>>A;for(auto&v:Az)A.push_back(asz(v));
 {TriangularIntegrator target(A,16);target.atomic_mode=true;seed_audit_phase="target_hermite";int count=0;for(auto&nu:req["request"]["targets"]){target.raw(asz(raw.at(label(nu.get<Nu>()))));if(++count%200==0)std::cerr<<"target "<<count<<" elapsed="<<seconds(start)<<"\n";}}
 {TriangularIntegrator supplement(A,16);supplement.atomic_mode=true;seed_audit_phase="none";supplement.complete(16);}
 json summaries=json::object(),entries=json::object();for(auto&[phase,rows]:seed_audit_entries){summaries[phase]=seed_summary(rows);entries[phase]=json::array();for(auto&[key,value]:rows)entries[phase].push_back(value);}
 auto&target=seed_audit_entries["target_hermite"];std::map<std::string,json>extra;for(auto phase:{"supplementary_outer","supplementary_hermite"})for(auto&[key,value]:seed_audit_entries[phase])if(!target.count(key))extra[key]=value;summaries["additional_to_target_hermite"]=seed_summary(extra);entries["additional_to_target_hermite"]=json::array();for(auto&[key,value]:extra)entries["additional_to_target_hermite"].push_back(value);
 json result={{"input",input},{"point",req["point"]},{"ordinary_targets",req["request"]["targets"].size()},{"FBI_masters",family.masters},{"completion_limit",16},{"summary",summaries},{"seeds",entries},{"seconds",seconds(start)},{"semantics","primitive basis e_i*z^degree/f^power or e_i*z^degree; scalar multiples deduplicated; target and supplementary runs have separate atomic caches"}};save_json(output,result);std::cerr<<summaries.dump(2)<<"\n";}
'''
source=source.replace('int main(int argc',extra+'\nint main(int argc',1)
source=source.replace('if(mode=="audit")audit_ranks','if(mode=="seeds")audit_seeds(argv[2],argv[3]);else if(mode=="audit")audit_ranks')
(out/'seed_audit.cpp').write_text(source)
subprocess.run(['g++','-O3','-std=c++17','-I'+str(out),'-I/usr/local/include',str(out/'seed_audit.cpp'),'-L/usr/local/lib','-Wl,-rpath,/usr/local/lib','-lflint','-lgmp','-o',str(here/'build/seed_audit')],check=True)
print(here/'build/seed_audit')
