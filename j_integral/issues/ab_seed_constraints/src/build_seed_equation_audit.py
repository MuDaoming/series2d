"""Build a read-only seed-equation variable auditor against production kernels.

The input seed manifest is produced by build_seed_audit.py. This tool does not
modify that auditor or the production source. Each primitive basis seed yields
its unreduced z-IBP row; P/Q variables are counted before Hermite simplification.
"""
from pathlib import Path
import shutil,subprocess
here=Path(__file__).resolve().parents[1]
production=here.parents[1]/'src/IBP/kernel'
out=here/'build/seed_equation_audit_kernel'
out.mkdir(parents=True,exist_ok=True)
for src in production.iterdir():
    if src.is_file():shutil.copyfile(src,out/src.name)
source=(here/'src/experiment.cpp').read_text().replace('#include "../../../src/IBP/kernel/point.cpp"','#include "point.cpp"')
extra=r'''
struct EquationVariableSet {
 std::set<std::string>finite,polynomial,boundary;
 long long equations=0,nonzero_equations=0,coefficient_nonzeros=0;
 void merge(const EquationVariableSet&other){++equations;if(!other.finite.empty()||!other.polynomial.empty()||!other.boundary.empty())++nonzero_equations;finite.insert(other.finite.begin(),other.finite.end());polynomial.insert(other.polynomial.begin(),other.polynomial.end());boundary.insert(other.boundary.begin(),other.boundary.end());coefficient_nonzeros+=other.finite.size()+other.polynomial.size()+other.boundary.size();}
 json summary()const{return {{"seed_equations",equations},{"nonzero_seed_equations",nonzero_equations},{"distinct_P_variables",finite.size()},{"distinct_Q_variables",polynomial.size()},{"distinct_P_plus_Q_variables",finite.size()+polynomial.size()},{"distinct_boundary_constants",boundary.size()},{"coefficient_nonzeros_before_elimination",coefficient_nonzeros}};}
 json variables()const{return {{"P",finite},{"Q",polynomial},{"boundary_constants",boundary}};}
};
std::string seed_key(json seed){seed.erase("calls");return seed.dump();}
void add_rational_variables(EquationVariableSet&out,int component,const ZR&value,Factors&factors){
 if(value.zero())return;P polynomial,proper;nmod_poly_divrem(polynomial.v,proper.v,value.n.v,value.d.v);
 for(int m=0;m<=polynomial.deg();++m)if(polynomial.at(m))out.polynomial.insert(json::array({"Q",component,m}).dump());
 for(auto&[f,k]:factors.get(value.d)){P fk=ppow(f,k),rest=quotient(value.d,fk);P digits=modp(proper*invert_mod(modp(rest,fk),fk),fk);for(int n=k;n>=1;--n){P hi,lo;nmod_poly_divrem(hi.v,lo.v,digits.v,f.v);for(int d=0;d<=lo.deg();++d)if(lo.at(d))out.finite.insert(json::array({"P",component,polyjson(f),n,d}).dump());digits=std::move(hi);}}
}
ZR seed_primitive(const json&seed){P numerator;numerator.set(seed.at("numerator_degree").get<int>(),1);if(seed.at("kind")=="polynomial")return zp(numerator);P f;for(size_t d=0;d<seed.at("factor").size();++d)f.set(d,seed["factor"][d].get<U>());return ZR(numerator,ppow(f,seed.at("power").get<int>()));}
void audit_seed_equations(const std::string&request_path,const std::string&seed_path,const std::string&output){
 auto req=load_json(request_path),manifest=load_json(seed_path);if(req.at("point")!=manifest.at("point"))throw std::runtime_error("request/seed point mismatch");init(req["point"]["prime"]);auto cfg=native_config(req);cfg["targets"]=json::array();auto start=Clock::now();NativeFamily family(cfg);if(manifest.at("FBI_masters")!=json(family.masters))throw std::runtime_error("seed/FBI basis mismatch");RMat At,Az;std::map<std::string,RRow>unused;family.model(At,Az,unused);std::vector<std::vector<ZR>>A;for(auto&v:Az)A.push_back(asz(v));TriangularIntegrator integrator(A,16);Factors factors;
 std::map<std::string,json>all_seeds;std::map<std::string,std::set<std::string>>phase_keys;for(auto it=manifest.at("seeds").begin();it!=manifest.at("seeds").end();++it){if(it.key()=="additional_to_target_hermite")continue;for(auto&seed:it.value()){auto key=seed_key(seed);all_seeds[key]=seed;phase_keys[it.key()].insert(key);}}
 for(auto phase:{"target_hermite","supplementary_outer","supplementary_hermite"})phase_keys["all_used_seeds"].insert(phase_keys[phase].begin(),phase_keys[phase].end());for(auto phase:{"supplementary_outer","supplementary_hermite"})phase_keys["supplementary_outer_and_internal"].insert(phase_keys[phase].begin(),phase_keys[phase].end());
 std::map<std::string,EquationVariableSet>memo;int done=0;for(auto&[key,seed]:all_seeds){int component=seed.at("component");ZR g=seed_primitive(seed);EquationVariableSet equation;for(int j=component;j<int(A.size());++j){ZR v=g*A[component][j];if(j==component)v=v+g.diff();add_rational_variables(equation,j,v,factors);}auto boundary=integrator.endpoint.finite_part(component,g);for(auto&[name,value]:boundary)if(value)equation.boundary.insert(name);memo.emplace(key,std::move(equation));if(++done%500==0)std::cerr<<"seed equations="<<done<<" elapsed="<<seconds(start)<<"\n";}
 family.cfg["targets"]=req["request"]["targets"];std::map<std::string,RRow>physical_rows;family.model(At,Az,physical_rows);EquationVariableSet target_variables;for(auto&[name,row]:physical_rows){EquationVariableSet one;for(int j=0;j<int(row.size());++j)add_rational_variables(one,j,asz(row[j]),factors);target_variables.merge(one);}
 json summaries=json::object(),variables=json::object();for(auto&[phase,keys]:phase_keys){EquationVariableSet aggregate;for(auto&key:keys)aggregate.merge(memo.at(key));summaries[phase]=aggregate.summary();variables[phase]=aggregate.variables();}
 auto target_summary=target_variables.summary();target_summary["target_rows"]=target_variables.equations;target_summary["nonzero_target_rows"]=target_variables.nonzero_equations;target_summary.erase("seed_equations");target_summary.erase("nonzero_seed_equations");summaries["target_density_representation"]=target_summary;variables["target_density_representation"]=target_variables.variables();EquationVariableSet combined;for(auto&key:phase_keys["all_used_seeds"])combined.merge(memo.at(key));combined.finite.insert(target_variables.finite.begin(),target_variables.finite.end());combined.polynomial.insert(target_variables.polynomial.begin(),target_variables.polynomial.end());auto combined_summary=combined.summary();combined_summary["ordinary_targets"]=physical_rows.size();combined_summary.erase("coefficient_nonzeros_before_elimination");summaries["all_seed_equations_plus_target_variables"]=combined_summary;variables["all_seed_equations_plus_target_variables"]=combined.variables();
 json result={{"request",request_path},{"seed_manifest",seed_path},{"point",req["point"]},{"FBI_rank",family.rank},{"summary",summaries},{"variables",variables},{"seconds",seconds(start)},{"semantics","unreduced nabla_z g coefficients are combined exactly, then partially fractioned over Fp; P=[component, monic irreducible factor, denominator power, numerator degree], Q=[component, monomial degree]; formal z0 finite-part constants counted separately, no physical endpoint constraints used; z1 declared regulated zero"}};save_json(output,result);std::cerr<<summaries.dump(2)<<"\n";
}
'''
assert source.count('int main(int argc')==1
source=source.replace('int main(int argc',extra+'\nint main(int argc',1)
source=source.replace('if(mode=="audit")audit_ranks','if(mode=="equation-vars"&&argc==5)audit_seed_equations(argv[2],argv[3],argv[4]);else if(mode=="audit")audit_ranks')
(out/'seed_equation_audit.cpp').write_text(source)
subprocess.run(['g++','-O3','-std=c++17','-I'+str(out),'-I/usr/local/include',str(out/'seed_equation_audit.cpp'),'-L/usr/local/lib','-Wl,-rpath,/usr/local/lib','-lflint','-lgmp','-o',str(here/'build/seed_equation_audit')],check=True)
print(here/'build/seed_equation_audit')
