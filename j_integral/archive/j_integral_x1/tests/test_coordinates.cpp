#define main unused_series_main
#include "../../j_integral/native/series.cpp"
#undef main
#include <cassert>
int main(int argc,char**argv){
 assert(argc==2);json cfg;{std::ifstream f(argv[1]);f>>cfg;}
 prime=cfg.at("prime");nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);
 NativeFamily family(cfg);RF t=Expression("t").parse(),z=Expression("z").parse();
 assert((family.X[0]-t).zero());assert((family.X[1]-(RF(1)-t)*z).zero());
 assert((family.jacobian-(RF(1)-t)).zero());
 assert((family.U-(t*(RF(1)-t)+rpow(RF(1)-t,2)*z*(RF(1)-z))).zero());
 // The equal-mass vacuum corner has G=(1-t)U^(-d/2) up to a constant.
 RMat At,Az;std::map<std::string,RRow> targets;family.model(At,Az,targets);
 assert(At.size()==1);RF h=family.dimension/RF(2);
 assert((At[0][0]+RF(1)/(RF(1)-t)+h*family.U.diff(0)/family.U).zero());
 assert((Az[0][0]+h*family.U.diff(1)/family.U).zero());
 // Equal masses make the corner invariant under y -> 1-y: its endpoint values coincide.
 RF u0=t*(RF(1)-t);
 assert((family.U-u0-rpow(RF(1)-t,2)*z*(RF(1)-z)).zero());
 auto invalid=cfg;invalid["jacobian"]="1-z";bool rejected=false;
 try{NativeFamily bad(invalid);}catch(const std::exception&){rejected=true;}assert(rejected);
 auto legacy=cfg;legacy.erase("branch_parameters");legacy.erase("jacobian");NativeFamily old(legacy);
 assert((old.jacobian-(RF(1)-z)).zero());
 std::cout<<"X1 coordinates, Jacobian, exact vacuum connection, invalid map rejection, legacy default: PASS\n";
}
