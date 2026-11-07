// Finite-field relation discovery on exact fixed-period coefficient coordinates.
#include <flint/nmod_poly.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using json=nlohmann::json;using U=mp_limb_t;static U prime;static nmod_t mod;
U add(U a,U b){return nmod_add(a,b,mod);}U sub(U a,U b){return nmod_sub(a,b,mod);}U mul(U a,U b){return nmod_mul(a,b,mod);}
struct Poly{nmod_poly_t v;Poly(){nmod_poly_init(v,prime);}Poly(const Poly&b):Poly(){nmod_poly_set(v,b.v);}Poly(Poly&&b)noexcept:Poly(){nmod_poly_swap(v,b.v);}Poly&operator=(const Poly&b){nmod_poly_set(v,b.v);return *this;}~Poly(){nmod_poly_clear(v);}bool zero()const{return nmod_poly_is_zero(v);}};
struct Rat{Poly n,d;Rat(){nmod_poly_one(d.v);}explicit Rat(const std::vector<U>&v):Rat(){for(size_t i=0;i<v.size();++i)nmod_poly_set_coeff_ui(n.v,i,v[i]);}bool zero()const{return n.zero();}void normal(){if(n.zero()){nmod_poly_one(d.v);return;}Poly g;nmod_poly_gcd(g.v,n.v,d.v);nmod_poly_div(n.v,n.v,g.v);nmod_poly_div(d.v,d.v,g.v);U c=nmod_inv(*nmod_poly_lead(d.v),mod);nmod_poly_scalar_mul_nmod(n.v,n.v,c);nmod_poly_scalar_mul_nmod(d.v,d.v,c);}};
Rat operator*(const Rat&a,const Rat&b){Rat x;nmod_poly_mul(x.n.v,a.n.v,b.n.v);nmod_poly_mul(x.d.v,a.d.v,b.d.v);x.normal();return x;}
Rat operator/(const Rat&a,const Rat&b){if(b.zero())throw std::runtime_error("zero rational pivot");Rat x;nmod_poly_mul(x.n.v,a.n.v,b.d.v);nmod_poly_mul(x.d.v,a.d.v,b.n.v);x.normal();return x;}
Rat operator-(const Rat&a,const Rat&b){Rat x;Poly y;nmod_poly_mul(x.n.v,a.n.v,b.d.v);nmod_poly_mul(y.v,b.n.v,a.d.v);nmod_poly_sub(x.n.v,x.n.v,y.v);nmod_poly_mul(x.d.v,a.d.v,b.d.v);x.normal();return x;}
std::vector<U> coeffs(const Poly&p){std::vector<U>v(nmod_poly_length(p.v));for(size_t i=0;i<v.size();++i)v[i]=nmod_poly_get_coeff_ui(p.v,i);return v;}
json ratjson(const Rat&a,bool negate=false){auto n=coeffs(a.n);if(negate)for(auto&c:n)if(c)c=prime-c;return {{"numerator_tau",n},{"denominator_tau",coeffs(a.d)}};}
using Coords=std::map<std::string,U>;using Series=std::vector<Coords>;
struct Echelon{
 int cols;std::map<int,std::vector<U>>rows;
 explicit Echelon(int n):cols(n){}
 void insert(std::vector<U>a){for(auto&[pc,row]:rows)if(a[pc]){U c=a[pc];for(int j=pc;j<cols;++j)a[j]=sub(a[j],mul(c,row[j]));}int pc=0;while(pc<cols&&!a[pc])++pc;if(pc==cols)return;U inv=nmod_inv(a[pc],mod);for(int j=pc;j<cols;++j)a[j]=mul(a[j],inv);rows.emplace(pc,std::move(a));}
 std::vector<std::vector<U>> kernel()const{std::vector<std::vector<U>>out;for(int f=0;f<cols;++f)if(!rows.count(f)){std::vector<U>v(cols);v[f]=1;for(auto it=rows.rbegin();it!=rows.rend();++it){U c=0;for(int j=it->first+1;j<cols;++j)c=add(c,mul(it->second[j],v[j]));v[it->first]=c?prime-c:0;}out.push_back(std::move(v));}return out;}
};
U at(const Series&s,int n,const std::string&key){if(n<0||n>=int(s.size()))return 0;auto it=s[n].find(key);return it==s[n].end()?0:it->second;}
struct Found {std::vector<std::vector<U>>relations;int rank,candidates,rejected,degree,train,order;};
Found search(const std::vector<Series>&s,int degree,int train,int order){int vars=s.size()*(degree+1);std::set<std::string>keys;for(auto&v:s)for(auto&row:v)for(auto&[key,c]:row)keys.insert(key);Echelon e(vars);for(int n=0;n<=train;++n)for(auto&key:keys){std::vector<U>row(vars);for(int i=0;i<int(s.size());++i)for(int d=0;d<=degree;++d)row[i*(degree+1)+d]=at(s[i],n-d,key);e.insert(std::move(row));}
 auto basis=e.kernel();std::vector<std::vector<U>>accepted;int rejected=0;
 for(auto&v:basis){bool ok=true;for(int n=train+1;n<=order&&ok;++n)for(auto&key:keys){U x=0;for(int i=0;i<int(s.size());++i)for(int d=0;d<=degree;++d)x=add(x,mul(v[i*(degree+1)+d],at(s[i],n-d,key)));if(x){ok=false;break;}}if(ok)accepted.push_back(v);else++rejected;}
 // If individual training basis vectors fail, intersect the entire training kernel
 // with holdout equations. This avoids losing a valid linear combination.
 if(rejected){Echelon restriction(basis.size());for(int n=train+1;n<=order;++n)for(auto&key:keys){std::vector<U>row(basis.size());for(int a=0;a<int(basis.size());++a)for(int i=0;i<int(s.size());++i)for(int d=0;d<=degree;++d)row[a]=add(row[a],mul(basis[a][i*(degree+1)+d],at(s[i],n-d,key)));restriction.insert(std::move(row));}accepted.clear();for(auto&w:restriction.kernel()){std::vector<U>v(vars);for(int a=0;a<int(w.size());++a)if(w[a])for(int j=0;j<vars;++j)v[j]=add(v[j],mul(w[a],basis[a][j]));accepted.push_back(std::move(v));}}
 return {accepted,int(e.rows.size()),int(basis.size()),rejected,degree,train,order};
}
struct Reduction {std::vector<int>masters;json reductions;};
Reduction reduce_relations(const Found&f,const std::vector<std::string>&labels){int m=labels.size();std::vector<std::vector<Rat>>rows;for(auto&v:f.relations){std::vector<Rat>row;for(int i=0;i<m;++i)row.emplace_back(std::vector<U>(v.begin()+i*(f.degree+1),v.begin()+(i+1)*(f.degree+1)));rows.push_back(std::move(row));}std::vector<int>pivots;int row=0;for(int c=0;c<m&&row<int(rows.size());++c){int k=row;while(k<int(rows.size())&&rows[k][c].zero())++k;if(k==int(rows.size()))continue;std::swap(rows[k],rows[row]);Rat v=rows[row][c];for(int j=c;j<m;++j)rows[row][j]=rows[row][j]/v;for(int i=0;i<int(rows.size());++i)if(i!=row&&!rows[i][c].zero()){Rat v=rows[i][c];for(int j=c;j<m;++j)rows[i][j]=rows[i][j]-v*rows[row][j];}pivots.push_back(c);++row;}
 std::vector<int>free;for(int i=0;i<m;++i)if(std::find(pivots.begin(),pivots.end(),i)==pivots.end())free.push_back(i);json out=json::object();for(int i=0;i<int(pivots.size());++i){json r=json::object();for(int j:free)if(!rows[i][j].zero())r[labels[j]]=ratjson(rows[i][j],true);out[labels[pivots[i]]]=r;}return {free,out};}
int complexity(const std::string&l){int sum=0,n=0;for(char c:l){if(std::isdigit(c))n=n*10+c-'0';else {sum+=n;n=0;}}return sum+n;}
int main(int argc,char**argv){try{
 if(argc<3)throw std::runtime_error("usage: j_search SERIES OUTPUT [--degree D --train N --de-degree D --j-only]");
 json in;{std::ifstream f(argv[1]);f>>in;}prime=in.at("prime");nmod_init(&mod,prime);int order=in.at("order"),degree=4,train=order*2/3,de_degree=16;bool jonly=false;for(int i=3;i<argc;++i){std::string a=argv[i];if(a=="--j-only"){jonly=true;continue;}if(i+1>=argc)throw std::runtime_error("missing option value");int v=std::stoi(argv[++i]);if(a=="--degree")degree=v;else if(a=="--train")train=v;else if(a=="--de-degree")de_degree=v;else throw std::runtime_error("unknown option");}if(train>=order||train<0)throw std::runtime_error("heldout window required");
 std::vector<std::string>labels,zero;for(auto it=in.at("coordinates").begin();it!=in.at("coordinates").end();++it){if(jonly&&it.key()[0]!='J')continue;bool nonzero=false;for(auto&row:it.value())if(!row.empty())nonzero=true;(nonzero?labels:zero).push_back(it.key());}
 // Eliminate high-dot J first; retain simple boundary corners and simple J.
 std::sort(labels.begin(),labels.end(),[](const std::string&a,const std::string&b){bool ba=a[0]=='B',bb=b[0]=='B';if(ba!=bb)return !ba;int ca=complexity(a),cb=complexity(b);return ca!=cb?ca>cb:a>b;});
 std::vector<Series>series;for(auto&l:labels){Series s;for(auto&row:in["coordinates"][l])s.push_back(row.get<Coords>());series.push_back(std::move(s));}
 Found f=search(series,degree,train,order);auto reduced=reduce_relations(f,labels);json out={{"prime",prime},{"order",order},{"train_order",train},{"heldout_orders",{train+1,order}},{"relation_degree",degree},{"train_rank",f.rank},{"training_kernel_dimension",f.candidates},{"rejected_training_basis_vectors",f.rejected},{"relations_after_holdout",f.relations.size()},{"zero_objects",zero},{"masters",json::array()},{"reductions",reduced.reductions},{"de",json::object()},{"input",argv[1]},{"tested_objects",labels}};
 for(int i:reduced.masters)out["masters"].push_back(labels[i]);
 for(int mi:reduced.masters){std::vector<Series>d;Series deriv;for(int n=0;n<order;++n){Coords c;for(auto&[k,v]:series[mi][n+1])if(U x=mul(n+1,v))c[k]=x;deriv.push_back(std::move(c));}d.push_back(deriv);for(int j:reduced.masters){Series v=series[j];v.pop_back();d.push_back(std::move(v));}bool found=false;for(int b=0;b<=de_degree&&!found;++b){Found df=search(d,b,std::min(train,order-2),order-1);for(auto&v:df.relations){std::vector<U>dv(v.begin(),v.begin()+b+1);Rat pivot(dv);if(pivot.zero())continue;json rel=json::object();for(int j=0;j<int(reduced.masters.size());++j){Rat a(std::vector<U>(v.begin()+(j+1)*(b+1),v.begin()+(j+2)*(b+1)));if(!a.zero())rel[labels[reduced.masters[j]]]=ratjson(a/pivot,true);}out["de"][labels[mi]]={{"degree",b},{"rhs",rel},{"verification_order",order-1}};found=true;break;}}if(!found)out["de"][labels[mi]]={{"status","unresolved"}};}
 out["master_count"]=reduced.masters.size();out["caveat"]="Finite-order coordinate relations; DE closure is tested separately. Boundary assumptions are inherited from coefficient input.";{std::ofstream f(argv[2]);f<<out.dump(2);}std::cout<<json({{"objects",labels.size()+zero.size()},{"masters",out["masters"]},{"de",out["de"]}}).dump(2)<<"\n";return 0;
 }catch(const std::exception&e){std::cerr<<"j_search: "<<e.what()<<"\n";return 1;}}
