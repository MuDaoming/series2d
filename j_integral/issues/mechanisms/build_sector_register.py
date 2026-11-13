"""Track what is actually understood, without treating raw sector labels as ranks."""
import json,collections,gzip
from pathlib import Path
base=Path(__file__).resolve().parent;root=base.parents[2]
source=json.load(open(root/'j_integral/issues/sector_symmetries/audit_with_cross_sector.json'))
fi=json.load(open(root/'j_integral/archive/j_integral_bootstrap/reference/fi_master_counts/dp.json'));FC=collections.Counter(''.join(str(int(v>0)) for v in nu[:8]) for nu in fi['master_indices'])
fam=json.load(open(root/'j_integral/tests/dp/input/family_request.json'));br=fam['j']['branches'];props=fam['family']['propagators'];rows=[]
threshold={q['mask'] for q in json.load(open(base/'threshold_bubble/dp_matches.json'))}
def counts(path):
 f=json.load(open(path)) if path.exists() else json.loads(gzip.decompress(path.with_suffix(path.suffix+'.gz').read_bytes()));return collections.Counter(''.join(str(int(int(v)>0)) for v in n.split('{')[1][:-1].split(',')) for n in f['basis'] if n.startswith('J'))
eold=counts(base/'branch_factorization/e6_full/j_original_dot2.json');enew=counts(base/'branch_factorization/e6_full/j_dot2.json');dp=counts(base/'dp_dot2_audit/j_result.json')
for x in source['e6']:
 assert eold[x['mask']]==x['J'];assert enew[x['mask']]==x['FI']
 rows.append({'family':'e6','mask':x['mask'],'FI_raw':x['FI'],'J_original_dot2':eold[x['mask']],'J_separated_dot2':enew[x['mask']],'status':'resolved','mechanism':'Two retained branches have zero intrinsic polynomials; explicit Beta separation makes J/FI reduction equivalent. Old projection retains t inside the kernel.','evidence':'branch_factorization/README.md'})
for x in source['dp']:
 mask=x['mask'];assert dp[mask]==x['J'];groups=[[i for i,q in enumerate(mask) if q=='1' and br[i]==b] for b in range(3)];single=[g[0] for g in groups if len(g)==1];bubble=[]
 for i in single:
  for j in single:
   if i<j and props[i]['mass2']==props[j]['mass2']:bubble.append([i+1,j+1])
 maps=[q['mask'] for q in x['equivalent_sectors'] if q['mask']!=mask];otherFI=[m for m in maps if FC[m]>0]
 status='unresolved';mechanism='No specific complete identity explaining the excess has been derived.';evidence='';remaining='Derive a concrete FI identity and show why it does not preserve the fixed-t slice.'
 if mask=='00100101':status='resolved';mechanism='m,0,m zero-invariant vacuum: homogeneity plus momentum IBP reduces I111 to tadpole product.';evidence='vacuum_m0m/README.md';remaining=''
 elif mask=='00110001':status='resolved';mechanism='Equal-mass sunrise at s=m2: full permutation plus disappearance of an explicit (s/m2-1)-weighted moment.';evidence='sunrise_deformation/README.md';remaining=''
 elif mask in ['01110001','01110010']:status='resolved';mechanism='Equal-mass massive-bubble exchange moves t; representative mass splitting increases FI3->4 while J stays4.';evidence='dp_massive_bubble/README.md';remaining=''
 elif mask in threshold:status='resolved';mechanism='A massive-massless bubble threshold coincides with a parent propagator; its explicit subloop DE and momentum IBP remove the topsector.';evidence='threshold_bubble/README.md';remaining=''
 elif x['FI']==0 and otherFI:status='representatives_pending';mechanism='FI representatives occur in an equivalent sector, so raw FI=0 is not topsector vanishing.';remaining='Merge equivalent sector representatives before diagnosing any physical J-FI excess.'
 elif bubble and len(groups[2])==1:status='partial';mechanism='Equal-mass bubble exchange involving the integrated branch exists; its contribution to the whole excess has not been completed.';remaining='Write the explicit transported FI relation and account for any remaining independent excess.'
 row={'family':'dp','mask':mask,'FI_raw':x['FI'],'J_original_dot2':dp[mask],'branch_sizes':list(map(len,groups)),'masses_squared':[props[i]['mass2'] for i,q in enumerate(mask) if q=='1'],'status':status,'mechanism':mechanism,'remaining':remaining,'equal_mass_singleton_pairs':bubble,'equivalent_sectors':maps,'equivalent_sectors_with_FI_representatives':otherFI,'evidence':evidence};rows.append(row)
summary={family:dict(collections.Counter(r['status'] for r in rows if r['family']==family)) for family in ['e6','dp']}
(base/'sector_register.json').write_text(json.dumps({'scope':'Raw sectors with J>FI in the saved lists; cross-sector representatives must not be mistaken for physical rank differences. New full dot<=2 J tests preserve these original counts.','summary':summary,'sectors':rows},indent=2)+'\n')
labels={'resolved':'已有具体解析机制','partial':'部分：bubble交换，差额未全解释','representatives_pending':'先合并跨sector代表','unresolved':'未解释'}
text='''# 仍需研究的 DP sector 清单

每个 mask 从左到右对应原 DP 的 D1 至 D8，完整定义见 [DP 家族输入](../../tests/dp/input/family_request.json)。

本表明确保留尚未解释的情况，不用“有 t-IBP”替代具体机制。原 FI 列是 Blade 清单的逐 sector 代表数，不是无条件的物理最小数；尤其 FI=0 且存在等价 sector 的 FI 代表时，不能称该顶层为零。

DP 已重新生成全部保留三个 branch、dot<=2 的 2952 个目标，仍为 496 维（481 普通 J、15 边界），DE 闭合；这些 sector 的原始 J 数均未变化。这个检查排除了本次少给 dot<=2 目标造成的数目变化，不是任意指标完备性定理。

完整状态、质量、等价 sector 和剩余任务见 [sector_register.json](sector_register.json)。e6 的 14 个正差 sector 已有统一的解析分离机制，见 [完整推导](branch_factorization/README.md)。复杂 DP 的多数 sector 仍未完成同等级的解释。

| DP sector | FI原代表数 | J | 当前状态 |
| --- | ---: | ---: | --- |
'''
for r in rows:
 if r['family']=='dp':text+=f"| {r['mask']} | {r['FI_raw']} | {r['J_original_dot2']} | {labels[r['status']]} |\n"
text+='\n“部分”仅表示存在明确的 bubble 对称关系，尚未证明它解释了全部差额；不计为已完整解决。跨 sector 代表状态也不计为一种新的积分机制。\n'
(base/'unresolved_sectors.md').write_text(text);print(json.dumps(summary,indent=2))
