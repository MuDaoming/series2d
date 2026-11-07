#define DIRECT_J_NO_MAIN
#include "../../j_integral_methods/native/direct_j.cpp"
#include <cassert>
int main(int argc,char**argv){assert(argc==2);json cfg;{std::ifstream f(argv[1]);f>>cfg;}prime=cfg["prime"];nmod_init(&mod,prime);nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);cfg["canonical_fbi_symmetries"]=true;cfg["specialize_t"]="17";NativeFamily family(cfg);int count=0;Nu source,target;for(auto&[a,b]:family.corner_aliases)if(a!=b){assert(family.constrained_corner_key(a)==family.constrained_corner_key(b));source=a;target=b;++count;}assert(count>0);
 int changed=-1;for(int i=0;i<family.n;++i)if(source[i]&&!target[i]){changed=i;break;}assert(changed>=0);auto mass=cfg["propagators"][changed]["mass2"];cfg["propagators"][changed]["mass2"]="("+(mass.is_string()?mass.get<std::string>():mass.dump())+")+1/11";NativeFamily different(cfg);assert(different.corner_aliases.at(source)!=different.corner_aliases.at(target));
 std::cout<<"Exact constrained FBI corner aliases="<<count<<"; changed-mass counterexample rejected PASS\n";
}
