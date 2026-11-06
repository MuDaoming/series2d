#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"
#include <cassert>
bool same(const std::pair<P,P>&a,const std::pair<P,P>&b){return (a.first*b.second-b.first*a.second).zero();}
int main(int argc,char**argv){assert(argc==2);json source;{std::ifstream f(argv[1]);f>>source;}prime=source["prime"];nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);
 json cfg=source;cfg["propagators"]=json::array();for(int i:{1,2,3,6})cfg["propagators"].push_back(source["propagators"][i]);cfg.erase("cut_sector");cfg["targets"]={{1,1,1,1},{2,1,1,1},{1,2,1,1},{0,1,1,1},{1,0,1,1}};
 NativeFamily symbolic(cfg);RMat At,Az;std::map<std::string,RRow>rows;symbolic.model(At,Az,rows);
 assert(symbolic.rank>1);for(U point:{17UL,19UL}){json fixed=cfg;fixed["specialize_t"]=point;NativeFamily numerical(fixed);assert(symbolic.masters==numerical.masters);RMat Bt,Bz;std::map<std::string,RRow>r;numerical.model(Bt,Bz,r);
  for(int i=0;i<symbolic.rank;++i)for(int j=0;j<symbolic.rank;++j){assert(same(specialize(At[i][j],point),specialize(Bt[i][j],point)));assert(same(specialize(Az[i][j],point),specialize(Bz[i][j],point)));}
  for(auto&values:cfg["targets"]){Nu nu=values.get<Nu>();auto d=numerical.target_derivative(nu,0),sd=symbolic.target_derivative(nu,0);for(int j=0;j<symbolic.rank;++j)assert(same(specialize(d[j],point),specialize(sd[j],point)));}
 }
 std::cout<<"Point-before-FBI versus symbolic-then-point: rank="<<symbolic.rank<<", both connections and target t derivatives PASS\n";
}
