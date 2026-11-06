#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"
#include <cassert>
bool equal_at(const RF&a,const RF&b,U point){auto x=specialize(a,point),y=specialize(b,point);return (x.first*y.second-y.first*x.second).zero();}
int main(int argc,char**argv){assert(argc==2);json cfg;{std::ifstream f(argv[1]);f>>cfg;}prime=cfg["prime"];nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);U point=17;
 json fixed=cfg;fixed["specialize_t"]=point;fixed["targets"]=json::array();NativeFamily value(fixed);RMat At,Az;std::map<std::string,RRow>unused;value.model(At,Az,unused);
 cfg["fbi_corner_aliases"]=json::array();for(auto&[a,b]:value.corner_aliases)cfg["fbi_corner_aliases"].push_back({a,b});first_jet_center=point;first_jet_mode=true;NativeFamily jet(cfg);assert(jet.masters==value.masters);std::vector<Nu>nus;Nu corner(cfg["propagators"].size(),1);nus.push_back(corner);Nu raised=corner;raised[0]+=3;nus.push_back(raised);for(auto&[sector,info]:value.sectors){nus.push_back(sector);Nu target=sector;for(int i=int(target.size())-1;i>=0;--i)if(target[i]){++target[i];break;}nus.push_back(target);}
 for(auto&nu:nus){auto reduced=jet.reduce(nu,0);RF prefactor=jet.prefactor(nu);std::vector<RF>row;for(int j=0;j<jet.rank;++j)row.push_back(prefactor*reduced[j]/jet.prefactor(jet.masters[j]));std::vector<std::pair<P,P>>deriv;for(auto&x:row)deriv.push_back(specialize(x.diff(0),point));std::vector<std::pair<P,P>>vals;for(auto&x:row)vals.push_back(specialize(x,point));first_jet_mode=false;auto exact=value.target_derivative(nu,0);
  for(int i=0;i<jet.rank;++i){RF d=RF(MP(0)); // convert specialized univariate polynomials back into RF.
   auto back=[](const std::pair<P,P>&p){MP n,den;for(int j=0;j<=p.first.deg();++j){ulong e[2]={0,ulong(j)};nmod_mpoly_set_coeff_ui_ui(n.v,p.first.at(j),e,mc);}for(int j=0;j<=p.second.deg();++j){ulong e[2]={0,ulong(j)};nmod_mpoly_set_coeff_ui_ui(den.v,p.second.at(j),e,mc);}return RF(std::move(n),std::move(den));};d=back(deriv[i]);for(int j=0;j<jet.rank;++j)d+=back(vals[j])*At[j][i];assert(equal_at(d,exact[i],point));}
  first_jet_mode=true;
 }
 first_jet_mode=false;std::cout<<"Dual-number FBI coefficient derivatives versus raised-target identity, including dot3: PASS rank="<<jet.rank<<"\n";
}
