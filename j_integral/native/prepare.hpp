#pragma once
// Generic two-loop input and FBI identities, in F_p[t,z].
// Follows expand's U redefinition and the same four FBI sector cases.
#include <flint/nmod_mpoly.h>
#include <cctype>
#include <numeric>
#include <set>
static nmod_mpoly_ctx_t mc;
struct MP {
 nmod_mpoly_t v;
 MP(){nmod_mpoly_init(v,mc);} explicit MP(U x):MP(){nmod_mpoly_set_ui(v,x%prime,mc);}
 MP(const MP& b):MP(){nmod_mpoly_set(v,b.v,mc);} MP(MP&& b)noexcept:MP(){nmod_mpoly_swap(v,b.v,mc);}
 MP& operator=(const MP& b){nmod_mpoly_set(v,b.v,mc);return *this;} MP& operator=(MP&& b)noexcept{nmod_mpoly_swap(v,b.v,mc);return *this;}
 ~MP(){nmod_mpoly_clear(v,mc);} bool zero()const{return nmod_mpoly_is_zero(v,mc);} bool one()const{return nmod_mpoly_is_one(v,mc);}
};
MP operator+(const MP&a,const MP&b){MP c;nmod_mpoly_add(c.v,a.v,b.v,mc);return c;}
MP operator-(const MP&a,const MP&b){MP c;nmod_mpoly_sub(c.v,a.v,b.v,mc);return c;}
MP operator*(const MP&a,const MP&b){MP c;nmod_mpoly_mul(c.v,a.v,b.v,mc);return c;}
MP mscale(const MP&a,U c){MP b;nmod_mpoly_scalar_mul_ui(b.v,a.v,c,mc);return b;}
MP mquot(const MP&a,const MP&b){MP c;if(!nmod_mpoly_divides(c.v,a.v,b.v,mc))throw std::runtime_error("nonexact multivariate polynomial quotient");return c;}
MP mgcd(const MP&a,const MP&b){MP c;if(!nmod_mpoly_gcd(c.v,a.v,b.v,mc))throw std::runtime_error("multivariate gcd failed");return c;}
MP mlcm(const MP&a,const MP&b){return mquot(a,mgcd(a,b))*b;}
struct RF {
 MP n,d;
 RF():n(),d(U(1)){} RF(int c):n(c<0?(prime-U(-c)%prime)%prime:U(c)%prime),d(U(1)){}
 explicit RF(MP p):n(std::move(p)),d(U(1)){}
 RF(MP n,MP d):n(std::move(n)),d(std::move(d)){normal();}
 bool zero()const{return n.zero();}
 void normal(){if(d.zero())throw std::runtime_error("zero denominator in FBI preprocessing");if(n.zero()){d=MP(1);return;}if(!d.one()){MP g=mgcd(n,d);if(!g.one()){n=mquot(n,g);d=mquot(d,g);}U lc=nmod_mpoly_get_term_coeff_ui(d.v,0,mc);if(lc!=1){U c=nmod_inv(lc,mod);n=mscale(n,c);d=mscale(d,c);}}}
 RF diff(int j)const{MP a,b;nmod_mpoly_derivative(a.v,n.v,j,mc);nmod_mpoly_derivative(b.v,d.v,j,mc);return RF(a*d-n*b,d*d);}
};
RF operator+(const RF&a,const RF&b){if(a.zero())return b;if(b.zero())return a;if(nmod_mpoly_equal(a.d.v,b.d.v,mc))return RF(a.n+b.n,a.d);MP g=mgcd(a.d,b.d),ad=mquot(a.d,g),bd=mquot(b.d,g);return RF(a.n*bd+b.n*ad,ad*b.d);}
RF operator-(const RF&a){return RF(mscale(a.n,prime-1),a.d);}
RF operator-(const RF&a,const RF&b){return a+(-b);}
RF operator*(const RF&a,const RF&b){if(a.zero()||b.zero())return RF();return RF(a.n*b.n,a.d*b.d);}
RF operator/(const RF&a,const RF&b){return RF(a.n*b.d,a.d*b.n);}
RF& operator+=(RF&a,const RF&b){a=a+b;return a;}
RF rpow(RF a,int n){if(n<0)return rpow(RF(a.d,a.n),-n);RF b(1);while(n){if(n&1)b=b*a;n/=2;if(n)a=a*a;}return b;}
struct Expression {
 std::string text;size_t pos=0;
 explicit Expression(std::string s):text(std::move(s)){}
 void ws(){while(pos<text.size()&&std::isspace(text[pos]))++pos;}
 bool eat(char c){ws();if(pos<text.size()&&text[pos]==c){++pos;return true;}return false;}
 RF atom(){ws();if(eat('(')){RF a=sum();if(!eat(')'))throw std::runtime_error("missing ')' in input expression");return a;}if(eat('t')||eat('z')){int j=text[pos-1]=='t'?0:1;MP p;nmod_mpoly_gen(p.v,j,mc);return RF(p);}U c=0;size_t start=pos;while(pos<text.size()&&std::isdigit(text[pos])){c=add(mul(c,10),U(text[pos++]-'0'));}if(pos==start)throw std::runtime_error("unsupported token in rational input: "+text.substr(pos));return RF(MP(c));}
 RF unary(){if(eat('+'))return unary();if(eat('-'))return -unary();RF a=atom();ws();bool power=eat('^');if(!power&&pos+1<text.size()&&text.substr(pos,2)=="**"){pos+=2;power=true;}if(power){bool neg=eat('-');ws();int n=0;size_t start=pos;while(pos<text.size()&&std::isdigit(text[pos]))n=n*10+text[pos++]-'0';if(pos==start)throw std::runtime_error("integer exponent required");a=rpow(a,neg?-n:n);}return a;}
 RF product(){RF a=unary();while(true){ws();if(pos<text.size()&&text[pos]=='*'&&(pos+1==text.size()||text[pos+1]!='*')){++pos;a=a*unary();}else if(eat('/'))a=a/unary();else break;}return a;}
 RF sum(){RF a=product();while(true){if(eat('+'))a=a+product();else if(eat('-'))a=a-product();else break;}return a;}
 RF parse(){RF a=sum();ws();if(pos!=text.size())throw std::runtime_error("unparsed expression: "+text.substr(pos));return a;}
};
RF readr(const json& x){return Expression(x.is_string()?x.get<std::string>():x.dump()).parse();}
using RRow=std::vector<RF>;using RMat=std::vector<RRow>;using Nu=std::vector<int>;
struct Elimination {RMat matrix;std::vector<int> pivots;};
Elimination eliminate(RMat a,int cols){int row=0;std::vector<int>pivots;for(int c=0;c<cols&&row<int(a.size());++c){int k=row;while(k<int(a.size())&&a[k][c].zero())++k;if(k==int(a.size()))continue;std::swap(a[k],a[row]);RF v=a[row][c];for(auto& x:a[row])x=x/v;for(int i=0;i<int(a.size());++i)if(i!=row&&!a[i][c].zero()){RF v=a[i][c];for(int j=0;j<int(a[i].size());++j)a[i][j]=a[i][j]-v*a[row][j];}pivots.push_back(c);++row;}return {std::move(a),std::move(pivots)};}
RMat inverse(const RMat& a){int n=a.size();RMat b=a;for(int i=0;i<n;++i){b[i].resize(2*n);b[i][n+i]=RF(1);}auto e=eliminate(std::move(b),n);if(int(e.pivots.size())!=n)return {};RMat x(n,RRow(n));for(int i=0;i<n;++i)for(int j=0;j<n;++j)x[i][j]=e.matrix[i][n+j];return x;}
struct SectorInfo {int type;std::vector<int>active;RMat inv;RF C;RRow zs;};
struct NativeFamily {
 json cfg;int n,rank;std::vector<int>branch;Nu cut;RF dimension,U,jacobian;std::vector<RF>X;RMat R;std::map<Nu,SectorInfo> sectors;std::vector<Nu>masters;std::map<Nu,int>master_index;std::map<std::pair<Nu,int>,RRow> cache;
 explicit NativeFamily(json cfg):cfg(std::move(cfg)){
  n=this->cfg.at("propagators").size();if(n>20)throw std::runtime_error("sector enumeration limited to 20 propagators");dimension=readr(this->cfg.at("dimension"));
  RF t=Expression("t").parse(),z=Expression("z").parse();X={t*(RF(1)-z),(RF(1)-t)*(RF(1)-z),z};jacobian=RF(1)-z;
  if(this->cfg.contains("branch_parameters")||this->cfg.contains("jacobian")){
   if(!this->cfg.contains("branch_parameters")||!this->cfg.contains("jacobian")||this->cfg.at("branch_parameters").size()!=3)throw std::runtime_error("provide three branch_parameters and jacobian together");
   for(int b=0;b<3;++b)X[b]=readr(this->cfg.at("branch_parameters")[b]);
   jacobian=readr(this->cfg.at("jacobian"));
   RF determinant=X[0].diff(0)*X[1].diff(1)-X[0].diff(1)*X[1].diff(0);
   if(!(X[0]+X[1]+X[2]-RF(1)).zero())throw std::runtime_error("branch_parameters must sum to one");
   if(jacobian.zero()||(!(jacobian-determinant).zero()&&!(jacobian+determinant).zero()))throw std::runtime_error("jacobian must equal the coordinate determinant up to orientation");
  }
  auto props=this->cfg.at("propagators");int E=this->cfg.at("gram").size();RMat gram(E,RRow(E)),c(n,RRow(2)),shift(n,RRow(E));RRow mass(n),kin(n);for(int a=0;a<E;++a)for(int b=0;b<E;++b)gram[a][b]=readr(this->cfg["gram"][a][b]);
  std::set<int> bs;std::map<int,RRow>bc;
  for(int i=0;i<n;++i){int b=props[i].at("branch");branch.push_back(b);bs.insert(b);if(props[i].at("loop").size()!=2||int(props[i].at("shift").size())!=E)throw std::runtime_error("invalid loop/shift length");for(int j=0;j<2;++j)c[i][j]=readr(props[i]["loop"][j]);for(int j=0;j<E;++j)shift[i][j]=readr(props[i]["shift"][j]);mass[i]=readr(props[i].at("mass2"));
   RRow quad={c[i][0]*c[i][0],c[i][0]*c[i][1],c[i][1]*c[i][1]};if(bc.count(b)){for(int j=0;j<3;++j)if(!(bc[b][j]-quad[j]).zero())throw std::runtime_error("branch quadratic loop parts differ");}else bc[b]=quad;
   kin[i]=-mass[i];for(int a=0;a<E;++a)for(int b=0;b<E;++b)kin[i]+=shift[i][a]*gram[a][b]*shift[i][b];
  }
  if(bs!=std::set<int>({0,1,2}))throw std::runtime_error("expected branches 0,1,2");
  RF aa,ab,bb;for(int b=0;b<3;++b){aa+=X[b]*bc[b][0];ab+=X[b]*bc[b][1];bb+=X[b]*bc[b][2];}U=aa*bb-ab*ab;if(U.zero())throw std::runtime_error("degenerate loop quadratic form");
  R.assign(n,RRow(n));for(int i=0;i<n;++i)for(int j=0;j<n;++j){RF dot;for(int a=0;a<E;++a)for(int b=0;b<E;++b)dot+=shift[i][a]*gram[a][b]*shift[j][b];RF adj=bb*c[i][0]*c[j][0]-ab*(c[i][0]*c[j][1]+c[i][1]*c[j][0])+aa*c[i][1]*c[j][1];R[i][j]=RF(2)*X[branch[i]]*X[branch[j]]*adj*dot;if(branch[i]==branch[j])R[i][j]=R[i][j]-U*X[branch[i]]*(kin[i]+kin[j]);}
  if(this->cfg.contains("cut_sector"))cut=this->cfg["cut_sector"].get<Nu>();
  if(!cut.empty()&&int(cut.size())!=n)throw std::runtime_error("invalid cut length");
  for(int mask=1;mask<(1<<n);++mask){Nu nu(n);for(int i=0;i<n;++i)nu[i]=(mask>>i)&1;if(valid(nu)&&keep(nu))sectors.emplace(nu,analyze(nu));}
  for(auto& [nu,info]:sectors)if(info.type==0)masters.push_back(nu);
  std::sort(masters.begin(),masters.end(),[](const Nu&a,const Nu&b){int sa=std::accumulate(a.begin(),a.end(),0),sb=std::accumulate(b.begin(),b.end(),0);return sa!=sb?sa>sb:a>b;});rank=masters.size();if(!rank)throw std::runtime_error("no retained MFBI");for(int i=0;i<rank;++i)master_index[masters[i]]=i;
 }
 bool valid(const Nu&nu)const{if(int(nu.size())!=n)return false;std::set<int>b;for(int i=0;i<n;++i){if(nu[i]<0)return false;if(nu[i])b.insert(branch[i]);}return b.size()==3;}
 bool keep(const Nu&nu)const{for(int i=0;i<int(cut.size());++i)if(cut[i]&&!nu[i])return false;return true;}
 SectorInfo analyze(const Nu&nu){std::vector<int>active;for(int i=0;i<n;++i)if(nu[i])active.push_back(i);int nn=active.size(),size=nn+3;RMat S(size,RRow(size));for(int i=0;i<nn;++i){S[branch[active[i]]][3+i]=S[3+i][branch[active[i]]]=RF(1);for(int j=0;j<nn;++j)S[3+i][3+j]=R[active[i]][active[j]];}
  // Eliminate affine branch constraints, as in the prototype's KKT identity.
  std::vector<int>piv(3,-1),free;for(int i=0;i<nn;++i){int b=branch[active[i]];if(piv[b]<0)piv[b]=i;else free.push_back(i);}int m=free.size();RMat H(m,RRow(m));for(int a=0;a<m;++a)for(int b=0;b<m;++b){int i=free[a],j=free[b],pi=piv[branch[active[i]]],pj=piv[branch[active[j]]];H[a][b]=S[3+i][3+j]-S[3+pi][3+j]-S[3+i][3+pj]+S[3+pi][3+pj];}RMat hi=inverse(H),invS;
  if(!m||!hi.empty()){RMat v(nn,RRow(m)),K(nn,RRow(nn)),Bm(nn,RRow(3));for(int a=0;a<m;++a){v[free[a]][a]=RF(1);v[piv[branch[active[free[a]]]]][a]=RF(-1);}for(int i=0;i<nn;++i)for(int j=0;j<nn;++j)for(int a=0;a<m;++a)for(int b=0;b<m;++b)K[i][j]+=v[i][a]*hi[a][b]*v[j][b];for(int i=0;i<nn;++i)for(int b=0;b<3;++b){Bm[i][b]=RF(i==piv[b]);for(int j=0;j<nn;++j)Bm[i][b]=Bm[i][b]-K[i][j]*S[3+j][3+piv[b]];}
   invS.assign(size,RRow(size));for(int a=0;a<3;++a){for(int b=0;b<3;++b)for(int j=0;j<nn;++j)invS[a][b]=invS[a][b]-S[3+piv[a]][3+j]*Bm[j][b];for(int j=0;j<nn;++j)invS[a][3+j]=invS[3+j][a]=Bm[j][a];}for(int i=0;i<nn;++i)for(int j=0;j<nn;++j)invS[3+i][3+j]=K[i][j];RRow cz(size);for(int i=0;i<size;++i)for(int b=0;b<3;++b)cz[i]+=invS[i][b];RF C=cz[0]+cz[1]+cz[2];return {C.zero()?1:0,active,std::move(invS),C,RRow(cz.begin()+3,cz.end())};
  }
  auto elim=eliminate(S,size);std::set<int>pc(elim.pivots.begin(),elim.pivots.end());RRow selected;RF C;for(int f=0;f<size;++f)if(!pc.count(f)){RRow v(size);v[f]=RF(1);for(int i=0;i<int(elim.pivots.size());++i)v[elim.pivots[i]]=-elim.matrix[i][f];RF cc=v[0]+v[1]+v[2];if(selected.empty()||!cc.zero()){selected=v;C=cc;}if(!cc.zero())break;}if(selected.empty())throw std::runtime_error("singular sector nullspace missing");return {C.zero()?3:2,active,{},C,RRow(selected.begin()+3,selected.end())};
 }
 void accumulate(RRow&out,const RF&c,const RRow&v){if(!c.zero())for(int j=0;j<rank;++j)out[j]+=c*v[j];}
 RRow reduce(const Nu&nu,int sh){if(!valid(nu)||!keep(nu))return RRow(rank);auto key=std::make_pair(nu,sh);auto it=cache.find(key);if(it!=cache.end())return it->second;Nu sec(n);for(int i=0;i<n;++i)sec[i]=nu[i]>0;auto& info=sectors.at(sec);RRow out(rank);int total=std::accumulate(nu.begin(),nu.end(),0);RF delta=dimension+RF(sh);auto lower=[&](int i){Nu a=nu;--a[i];return a;};
  if(info.type==0){if(*std::max_element(nu.begin(),nu.end())<=1){if(!sh)out[master_index.at(nu)]=RF(1);else if(sh>0){RF f=RF(2)*delta-RF(total+3);accumulate(out,info.C/f,reduce(nu,sh-1));for(int j=0;j<int(info.active.size());++j)accumulate(out,-info.zs[j]/f,reduce(lower(info.active[j]),sh-1));}else{RF f=RF(2)*(delta+RF(1))-RF(total+3);accumulate(out,f/info.C,reduce(nu,sh+1));for(int j=0;j<int(info.active.size());++j)accumulate(out,info.zs[j]/info.C,reduce(lower(info.active[j]),sh));}}
   else{int i=info.active[0];for(int j:info.active)if(nu[j]>nu[i])i=j;Nu prev=lower(i);int a=std::find(info.active.begin(),info.active.end(),i)-info.active.begin();const auto&row=info.inv[3+a];RF c;for(int b=0;b<3;++b)c+=row[b];accumulate(out,-c/RF(nu[i]-1),reduce(prev,sh-1));for(int j=0;j<int(info.active.size());++j){Nu src=prev;--src[info.active[j]];accumulate(out,row[3+j]/RF(nu[i]-1),reduce(src,sh-1));}}
  }else if(info.type==1){for(int j=0;j<int(info.active.size());++j)accumulate(out,-info.zs[j]/(RF(2)*delta-RF(total+3)),reduce(lower(info.active[j]),sh-1));}
  else if(info.type==2){for(int j=0;j<int(info.active.size());++j)accumulate(out,info.zs[j]/info.C,reduce(lower(info.active[j]),sh));}
  else{int a=0;while(info.zs[a].zero())++a;for(int j=0;j<int(info.active.size());++j)if(j!=a){Nu src=nu;++src[info.active[a]];--src[info.active[j]];accumulate(out,-info.zs[j]/info.zs[a],reduce(src,sh));}}
  cache.emplace(key,out);return out;
 }
 RRow diff(const Nu&nu,int variable){RRow out(rank);for(int i=0;i<n;++i)if(nu[i])for(int j=0;j<n;++j)if(nu[j]){RF c=R[i][j].diff(variable);if(!c.zero()){Nu src=nu;++src[i];++src[j];accumulate(out,-c*RF(nu[i]*(nu[j]+(i==j)))/RF(2),reduce(src,1));}}return out;}
 RF prefactor(const Nu&nu){RF a=jacobian;for(int b=0;b<3;++b){int k=-1;for(int i=0;i<n;++i)if(branch[i]==b)k+=nu[i];a=a*rpow(X[b],k);}return a*rpow(U,std::accumulate(nu.begin(),nu.end(),0));}
 void model(RMat&At,RMat&Az,std::map<std::string,RRow>&targets){RRow p;for(auto&nu:masters)p.push_back(prefactor(nu));At.assign(rank,RRow(rank));Az=At;for(int v=0;v<2;++v)for(int i=0;i<rank;++i){auto row=diff(masters[i],v);for(int j=0;j<rank;++j)row[j]=p[i]*row[j]/p[j];row[i]+=p[i].diff(v)/p[i]-RF(3)*dimension*U.diff(v)/(RF(2)*U);(v?Az:At)[i]=std::move(row);}for(auto& values:cfg.at("targets")){Nu nu=values.get<Nu>();if(!valid(nu))throw std::runtime_error("invalid target indices");RRow row=reduce(nu,0);RF fac=prefactor(nu);std::string label="J{";for(int i=0;i<n;++i)label+=(i?",":"")+std::to_string(nu[i]);label+="}";for(int j=0;j<rank;++j)row[j]=fac*row[j]/p[j];targets[label]=std::move(row);}}
};
P univariate(const MP& a){P p;for(slong i=0;i<nmod_mpoly_length(a.v,mc);++i){ulong e[2];nmod_mpoly_get_term_exp_ui(e,a.v,i,mc);if(e[0])throw std::runtime_error("expected specialized polynomial");p.set(e[1],nmod_mpoly_get_term_coeff_ui(a.v,i,mc));}return p;}
std::vector<P> shiftjet(const MP& a,U t0){std::vector<P>out(1);for(slong i=0;i<nmod_mpoly_length(a.v,mc);++i){ulong e[2];nmod_mpoly_get_term_exp_ui(e,a.v,i,mc);if(e[0]>=prime)throw std::runtime_error("polynomial degree exceeds characteristic");out.resize(std::max(out.size(),size_t(e[0]+1)));U c=nmod_mpoly_get_term_coeff_ui(a.v,i,mc),bin=1;for(ulong j=0;j<=e[0];++j){U v=mul(mul(c,bin),nmod_pow_ui(t0,e[0]-j,mod));out[j].set(e[1],add(out[j].at(e[1]),v));if(j<e[0])bin=nmod_div(mul(bin,e[0]-j),j+1,mod);}}return out;}
json polyjson(const P&a){json x=json::array();for(int j=0;j<=a.deg();++j)x.push_back(a.at(j));return x;}
json jetjson(const std::vector<P>& a){json x=json::array();for(auto& p:a)x.push_back(polyjson(p));return x;}
U constant(const RF&a){if(nmod_mpoly_length(a.n.v,mc)>1||nmod_mpoly_length(a.d.v,mc)>1)throw std::runtime_error("expected scalar");ulong e[2];if(!a.n.zero()){nmod_mpoly_get_term_exp_ui(e,a.n.v,0,mc);if(e[0]||e[1])throw std::runtime_error("expected scalar");}nmod_mpoly_get_term_exp_ui(e,a.d.v,0,mc);if(e[0]||e[1])throw std::runtime_error("expected scalar");U n=a.n.zero()?0:nmod_mpoly_get_term_coeff_ui(a.n.v,0,mc),d=nmod_mpoly_get_term_coeff_ui(a.d.v,0,mc);return nmod_div(n,d,mod);}
json prepare(const json&cfg){auto before=Clock::now();nmod_mpoly_ctx_init(mc,2,ORD_LEX,prime);RMat At,Az;std::map<std::string,RRow>targets;json metadata;
 if(cfg.at("kind")=="two_loop_family"){json full=cfg;if(cfg.contains("max_dot")){int n=cfg.at("propagators").size(),dot=cfg.at("max_dot");if(dot<0)throw std::runtime_error("negative max_dot");full["targets"]=json::array();Nu nu(n,1);std::function<void(int,int)>visit=[&](int i,int left){if(i==n){full["targets"].push_back(nu);return;}for(int k=0;k<=left;++k){nu[i]=1+k;visit(i+1,left-k);}nu[i]=1;};visit(0,dot);}NativeFamily family(full);family.model(At,Az,targets);metadata={{"fbi_masters",family.masters},{"cut_sector",family.cut},{"dimension",cfg.at("dimension")},{"target_indices",full.at("targets")}};}
 else if(cfg.at("kind")=="fbi_connection"){int r=cfg.at("At").size();At.assign(r,RRow(r));Az=At;for(int i=0;i<r;++i)for(int j=0;j<r;++j){At[i][j]=readr(cfg["At"][i][j]);Az[i][j]=readr(cfg["Az"][i][j]);}for(auto it=cfg.at("targets").begin();it!=cfg.at("targets").end();++it){RRow v;for(auto&x:it.value())v.push_back(readr(x));targets[it.key()]=std::move(v);}}
 else throw std::runtime_error("unknown input kind");
 int r=At.size();if(!r||Az.size()!=At.size())throw std::runtime_error("empty or mismatched connection");
 // Exact flatness validation in F_p(t,z), independent of J.
 for(int i=0;i<r;++i)for(int j=0;j<r;++j){RF f=At[i][j].diff(1)-Az[i][j].diff(0);for(int k=0;k<r;++k)f+=At[i][k]*Az[k][j]-Az[i][k]*At[k][j];if(!f.zero())throw std::runtime_error("FBI connection flatness failed");}
 U t0=constant(readr(cfg.at("t0")));MP qt(1);for(auto&row:At)for(auto&x:row)qt=mlcm(qt,x.d);json out={{"prime",prime},{"rank",r},{"order",cfg.at("order")},{"verify_all",cfg.value("verify_all",false)},{"include_z_boundaries",cfg.value("include_z_boundaries",false)},{"qt",jetjson(shiftjet(qt,t0))},{"targets",json::array()},{"nt",json::array()},{"endpoints",{cfg.at("endpoints").at("0"),cfg.at("endpoints").at("1")}}};
 for(auto& row:At){json rr=json::array();for(auto&x:row)rr.push_back(jetjson(shiftjet(x.n*mquot(qt,x.d),t0)));out["nt"].push_back(rr);}
 P q(1);auto addden=[&](const P&d){if(d.zero())throw std::runtime_error("nonordinary t0");P g;nmod_poly_gcd(g.v,q.v,d.v);q=quotient(q,g)*d;};
 addden(shiftjet(qt,t0)[0]);
 std::vector<std::vector<std::pair<P,P>>> az0(r);for(int i=0;i<r;++i)for(auto&x:Az[i]){P n=shiftjet(x.n,t0)[0],d=shiftjet(x.d,t0)[0];P g;nmod_poly_gcd(g.v,n.v,d.v);n=quotient(n,g);d=quotient(d,g);addden(d);az0[i].push_back({std::move(n),std::move(d)});}
 for(auto&[label,row]:targets){if(int(row.size())!=r)throw std::runtime_error("target rank mismatch");MP den(1);for(auto&x:row)den=mlcm(den,x.d);auto dj=shiftjet(den,t0);addden(dj[0]);json jj=json::array();for(auto&x:row)jj.push_back(jetjson(shiftjet(x.n*mquot(den,x.d),t0)));out["targets"].push_back({{"label",label},{"den",jetjson(dj)},{"numer",jj}});}
 P g;nmod_poly_gcd(g.v,q.v,derivative(q).v);q=quotient(q,g);nmod_poly_make_monic(q.v,q.v);if(q.deg()<1)throw std::runtime_error("polynomial-only connection unsupported");out["q"]=polyjson(q);out["B"]=json::array();
 for(auto&row:az0){json rr=json::array();for(auto&[n,d]:row){P b=quotient(q*n,d);if(b.deg()>=q.deg())throw std::runtime_error("non-Fuchsian/infinity polynomial connection unsupported");rr.push_back(polyjson(b));}out["B"].push_back(rr);}
 metadata["branch_parameters"]=cfg.value("branch_parameters",json::array({"t*(1-z)","(1-t)*(1-z)","z"}));metadata["jacobian"]=cfg.value("jacobian",json("1-z"));
 out["metadata"]=metadata;out["t0"]=cfg.at("t0");out["preparation_seconds"]=seconds(before);std::cerr<<"FBI preparation "<<seconds(before)<<"s rank="<<r<<" q_degree="<<q.deg()<<"\n";return out;
}
