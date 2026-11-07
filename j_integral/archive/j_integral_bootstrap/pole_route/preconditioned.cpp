#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"
#include <filesystem>
P jp(const json&v){P p;for(size_t i=0;i<v.size();++i)p.set(i,v[i].get<U>());return p;}
bool precfit(const Interpolator&I,const std::vector<U>&y,const std::vector<U>&xs,const P&D,Fraction&out,const P&N=P(1)){
 std::vector<U>v;for(size_t i=0;i<xs.size();++i){U nv=nmod_poly_evaluate_nmod(N.v,xs[i]);if(!nv)return false;v.push_back(nmod_div(mul(y[i],nmod_poly_evaluate_nmod(D.v,xs[i])),nv,mod));}
 P a;for(size_t i=0;i<I.basis.size();++i)a+=scale(I.basis[i],v[i]);P r0=I.modulus,r1=a,s0,s1(1);
 while(!r1.zero()){
  if(r1.deg()+s1.deg()<int(I.basis.size())){
   bool ok=true;for(size_t k=I.basis.size();k<xs.size();++k){U d=nmod_poly_evaluate_nmod(s1.v,xs[k]);if(!d||nmod_poly_evaluate_nmod(r1.v,xs[k])!=mul(v[k],d)){ok=false;break;}}
   if(ok){out.n=r1*N;out.d=s1*D;out.normal();for(size_t k=0;k<xs.size();++k){U d=nmod_poly_evaluate_nmod(out.d.v,xs[k]);if(!d||nmod_poly_evaluate_nmod(out.n.v,xs[k])!=mul(y[k],d))return false;}return true;}
  }
  P q,rem;nmod_poly_divrem(q.v,rem.v,r0.v,r1.v);P sn=s0-q*s1;r0=std::move(r1);r1=std::move(rem);s0=std::move(s1);s1=std::move(sn);
 }return false;
}
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("preconditioned MANIFEST OUTPUT");json cfg;{std::ifstream f(argv[1]);f>>cfg;}prime=cfg.at("prime").get<U>();nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);auto start=Clock::now();
 std::string section=cfg.at("section"),output=section=="DE_values"?"de":"reductions";U center=constant(readr(cfg.at("t0")));auto masters=cfg.at("masters");std::set<std::string>masterset;for(auto&m:masters)masterset.insert(m.get<std::string>());int holdout=cfg.value("holdout",3);bool ratios=cfg.value("ratios",false);int ratio_budget=cfg.value("ratio_budget",12);int max_den_degree=cfg.value("max_denominator_degree",512);
 std::map<U,json>points;for(auto&fn:cfg.at("points")){json p;{std::ifstream f(fn.get<std::string>());f>>p;}
  if(p.at("prime").get<U>()!=prime||p.at("masters_at_point")!=masters)throw std::runtime_error("point field/basis mismatch: "+fn.get<std::string>());
  U t=constant(readr(p.at("t_point")));if(!points.count(t))points[t]=json::object();
  for(auto r=p.at(section).begin();r!=p.at(section).end();++r){if(points[t].contains(r.key())&&points[t][r.key()]!=r.value())throw std::runtime_error("overlapping numerical row differs");points[t][r.key()]=r.value();}}
 std::map<std::string,json>known;for(auto&fn:cfg.at("known")){json p;{std::ifstream f(fn.get<std::string>());f>>p;}
  if(p.at("prime").get<U>()!=prime||constant(readr(p.at("t0")))!=center)throw std::runtime_error("known field/tau mismatch: "+fn.get<std::string>());
  if(p.contains("masters")&&p.at("masters")!=masters)throw std::runtime_error("known basis mismatch");
  if(!p.contains(output))continue;
  for(auto r=p.at(output).begin();r!=p.at(output).end();++r){json v=output=="de"?r.value().at("rhs"):r.value();for(auto c=v.begin();c!=v.end();++c){
   if(!masterset.count(c.key()))throw std::runtime_error("known column absent from declared basis: "+c.key());
   Fraction rf;rf.n=jp(c.value().at("numerator_tau"));rf.d=jp(c.value().at("denominator_tau"));if(rf.d.zero())throw std::runtime_error("known denominator zero");rf.normal();json nv=rf.value();
   if(known[r.key()].contains(c.key())&&known[r.key()][c.key()]!=nv)throw std::runtime_error("overlapping rational coefficient differs: "+r.key()+" / "+c.key());known[r.key()][c.key()]=nv;}}}
 struct CandidateWitness {std::string source,row,col;json value;};std::vector<CandidateWitness>external_witnesses;
 if(cfg.contains("witness_only"))for(auto&source:cfg.at("witness_only")){
  std::string path=source.at("path"),kind=source.at("section");if(kind!="de"&&kind!="reductions")throw std::runtime_error("witness section must be de or reductions");json w;{std::ifstream f(path);if(!f)throw std::runtime_error("cannot read witness source");f>>w;}
  if(w.at("prime").get<U>()!=prime||constant(readr(w.at("t0")))!=center)throw std::runtime_error("witness field/tau mismatch");if(w.contains("masters")&&w.at("masters")!=masters)throw std::runtime_error("witness basis mismatch");
  for(auto r=w.at(kind).begin();r!=w.at(kind).end();++r){json vs=kind=="de"?r.value().at("rhs"):r.value();for(auto c=vs.begin();c!=vs.end();++c){if(!masterset.count(c.key()))throw std::runtime_error("witness column absent from declared basis");Fraction f;f.n=jp(c.value().at("numerator_tau"));f.d=jp(c.value().at("denominator_tau"));if(f.d.zero())throw std::runtime_error("witness denominator zero");f.normal();if(!f.n.zero())external_witnesses.push_back({path,r.key(),c.key(),f.value()});}}
 }
 std::map<std::string,std::map<std::string,json>> column_values;for(auto&[r,vs]:known)for(auto c=vs.begin();c!=vs.end();++c)column_values[c.key()][r]=c.value();
 using Counts=std::map<std::string,int>;std::map<std::string,Counts>rd,cd;std::map<std::string,P>denoms;
 auto addknown=[&](const std::string&r,const std::string&c,const json&v){std::string d=v.at("denominator_tau").dump();rd[r][d]++;cd[c][d]++;if(!denoms.count(d))denoms.emplace(d,jp(v.at("denominator_tau")));};
 for(auto&[r,vs]:known)for(auto c=vs.begin();c!=vs.end();++c)addknown(r,c.key(),c.value());
 // Witness functions contribute candidates only. They never enter known or skip a target.
 for(size_t i=0;i<external_witnesses.size();++i){auto&w=external_witnesses[i];std::string d=w.value.at("denominator_tau").dump();cd[w.col][d]++;if(!denoms.count(d))denoms.emplace(d,jp(w.value.at("denominator_tau")));column_values[w.col]["witness_only:"+std::to_string(i)]=w.value;}

 auto top=[&](const Counts&cs,int n){std::vector<std::pair<int,std::string>>v;for(auto&[d,k]:cs)v.push_back({-k,d});std::sort(v.begin(),v.end());std::vector<std::string>o;for(auto&[k,d]:v){if(int(o.size())>=n)break;o.push_back(d);}return o;};
 std::set<std::string>rows;for(auto&[t,p]:points)for(auto r=p.begin();r!=p.end();++r)rows.insert(r.key());
 json out={{"prime",prime},{"t0",cfg.at("t0")},{"masters",masters},{"internal_holdout",holdout},{"witness_only_coefficients",external_witnesses.size()},{output,json::object()}};long long trials=0,passes=0;int rowno=0;
 auto save=[&](){out["new_count"]=passes;out["trials"]=trials;out["seconds"]=seconds(start);std::string tmp=std::string(argv[2])+".tmp";{std::ofstream f(tmp);f<<out.dump();}std::filesystem::rename(tmp,argv[2]);};
 for(auto&row:rows){std::vector<U>xs;std::vector<const json*>pr;std::set<std::string>cols;for(auto&[t,p]:points)if(p.contains(row)){xs.push_back(sub(t,center));pr.push_back(&p.at(row));for(auto c=p.at(row).begin();c!=p.at(row).end();++c)cols.insert(c.key());}if(xs.size()<=size_t(holdout))continue;Interpolator I(std::vector<U>(xs.begin(),xs.end()-holdout));
  auto rowds=top(rd[row],30);
  for(auto&col:cols){if(known[row].contains(col))continue;std::vector<U>v;for(auto*p:pr)v.push_back(p->value(col,U(0)));if(std::all_of(v.begin(),v.end(),[](U x){return !x;}))continue;
   auto cds=top(cd[col],20);std::vector<std::string>ds=rowds;ds.insert(ds.end(),cds.begin(),cds.end());for(size_t i=0;i<std::min(size_t(5),rowds.size());++i)for(size_t j=0;j<std::min(size_t(5),cds.size());++j){const P&a=denoms.at(rowds[i]);const P&b=denoms.at(cds[j]);P g;nmod_poly_gcd(g.v,a.v,b.v);P d=quotient(a*b,g);std::string key=polyjson(d).dump();if(!denoms.count(key))denoms.emplace(key,d);ds.push_back(key);}
   if(ratios){
    std::vector<json>witnesses;for(auto c=known[row].begin();c!=known[row].end();++c)witnesses.push_back(c.value());
    for(auto&[r,co]:column_values[col])witnesses.push_back(co);
    std::sort(witnesses.begin(),witnesses.end(),[](const json&a,const json&b){return a.at("numerator_tau").size()>b.at("numerator_tau").size();});
    std::set<std::string>used;int tried=0;bool solved=false;
    for(auto&w:witnesses){if(!used.insert(w.dump()).second)continue;if(tried++>=ratio_budget)break;
     P N=jp(w.at("numerator_tau")),D=jp(w.at("denominator_tau"));if(N.zero())continue;++trials;Fraction f;
     if(!precfit(I,v,xs,D,f,N))continue;json value=f.value();if(output=="de")out[output][row]["rhs"][col]=value;else out[output][row][col]=value;
     known[row][col]=value;column_values[col][row]=value;addknown(row,col,value);++passes;solved=true;break;
    }
    if(solved)continue;
   }
   std::set<std::string>seen;for(auto&d:ds){if(!seen.insert(d).second||denoms.at(d).deg()>max_den_degree)continue;++trials;Fraction f;if(!precfit(I,v,xs,denoms.at(d),f))continue;json value=f.value();if(output=="de")out[output][row]["rhs"][col]=value;else out[output][row][col]=value;known[row][col]=value;addknown(row,col,value);++passes;break;}
  }
  if(++rowno%20==0){save();std::cerr<<"rows="<<rowno<<" new="<<passes<<" trials="<<trials<<" seconds="<<seconds(start)<<"\n";}
 }save();std::cout<<"DONE new="<<passes<<" trials="<<trials<<" seconds="<<seconds(start)<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
