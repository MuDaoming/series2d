import json,gzip,collections,hashlib
from pathlib import Path
from analytic_bootstrap.upstream.model import read_matrix
b=Path(__file__).resolve().parent;root=b.parents[3]
load=lambda p:json.load(open(p))
def size(tag):
 f=load(b/f'j_{tag}_common_prime.json');a=f['diagnostics']['linear_algebra'];systems=a['FBI_elimination_systems'];return dict(basis=len(f['basis']),prime=f['point']['prime'],FBI_rank=f['diagnostics']['FBI_rank'],FBI_systems=len(systems),FBI_equation_rows_sum=sum(s['equations'] for s in systems),FBI_max_equations=max(s['equations'] for s in systems),FBI_shapes=[dict(rows=k[0],columns=k[1],count=v) for k,v in sorted(collections.Counter((s['equations'],s['coefficient_columns']) for s in systems).items())],**{k:v for k,v in a.items() if k!='FBI_elimination_systems'})
def ab(variant,n,anchor_calls):
 folder=b/variant;anchor=json.loads(gzip.decompress((folder/'anchor.json.gz').read_bytes()));result=load(folder/'result.json');a=read_matrix(anchor['reductions']['dot2'],1000003,(1111,n));r=read_matrix(result['matrix'],1000003,(1111,n));assert a.entries==r.entries
 events=[json.loads(line) for line in (folder/'numeric/events.jsonl').read_text().splitlines()];calls=[e for e in events if e['operation']=='IBP-complete'];assert len({e['t'] for e in calls})==len(calls);assert not any(e['operation']=='IBP-failed' for e in events)
 phases=collections.Counter(e['mode'] for e in calls);check=result['checks'];assert check['status']=='unique-with-holdouts'
 return dict(actual_IBP_calls=len(calls),DE_stage=phases['de'],reduction_stage_new=phases['reduction'],reduction_training_points=len(check['training_points']),reduction_holdout_points=len(check['holdout_points']),learning_anchor_IBP_calls=anchor_calls,total_including_this_reduction_anchor=len(calls)+anchor_calls,excludes='Historical DE-anchor acquisition, if any. The separated DE anchor is analytic.',rational_entries_verified=len(a.entries),ansatz_dimension=check['dimension'],rejected_points=check['rejected_points'])
fi=load(b/'fi/run_rank1/statistics.json');fi_check=load(b/'verification_fi_j.json');assert fi_check['all_FI_vs_J_reductions_equal'];d=load(b/'j_direct/statistics.json');ds=load(b/'j_direct_separated/statistics.json')
for directory,stat in [('j_direct',d),('j_direct_separated',ds)]:
 events=[json.loads(x) for x in (b/directory/'events.jsonl').read_text().splitlines()];assert len(events)==stat['actual_IBP_calls']==len({e['t'] for e in events})
report={'scope':{'all_dot_le_2_targets':1695,'common_three_branch_targets':1111,'additional_scaleless_zero_targets':584,'dimension':'26/7','kinematics':'unchanged e6 input; all FI parameters numerical','size_comparison_prime':fi['prime'],'sampling_comparison_prime':1000003,'sampling_variable':'tau=t-2/5','count_rule':'Only actual new public J-IBP executions; cache hits and repeated block requests are excluded. AB learning-anchor costs are reported separately.'},'FI':fi,'J_original_size':size('original'),'J_separated_size':size('separated'),'J_original_direct':d,'J_separated_direct_generic':ds,'J_separated_direct_using_proved_factorization':{'IBP_calls':1,'note':'The complete rational t dependence follows from the proved factorization, so no t reconstruction is needed.'},'J_original_AB':ab('j_ab',79,148),'J_separated_AB':ab('j_ab_separated',61,1),'FI_J_full_numeric_verification':fi_check,'limitations':['The unmodified AB second-layer existing-DE route hit KeyError:block0 in residue witness setup before any IBP. We used its default reduce flow, which establishes/checks the DE; DE and reduction point counts are split above.','An initial AB interface variant unnecessarily computed derivatives on reduction requests. It is retained locally as a diagnostic and excluded from the final run counts. The final interface requests derivatives only in DE mode.','FI metadata final AssociateTo calls had a logging arity error; saved matrix dimensions were validated and the metadata recovered without repeating the numerical solve.']}
refs=b/'reference';refs.mkdir(exist_ok=True)
for source,name in [(b/'fi/run_rank1/reductions.json','fi_reductions.json.gz'),(b/'fi/run_rank1/worker.log','fi_worker.log.gz'),(b/'j_original_common_prime.json','j_original_common_prime.json.gz'),(b/'j_separated_common_prime.json','j_separated_common_prime.json.gz'),(b/'j_direct/matrix.json','j_original_direct_matrix.json.gz'),(b/'j_direct/events.jsonl','j_original_direct_events.jsonl.gz'),(b/'j_direct_separated/matrix.json','j_separated_direct_matrix.json.gz'),(b/'j_direct_separated/events.jsonl','j_separated_direct_events.jsonl.gz'),(b/'j_ab/result.json','j_original_AB_result.json.gz'),(b/'j_ab/numeric/events.jsonl','j_original_AB_events.jsonl.gz'),(b/'j_ab_separated/result.json','j_separated_AB_result.json.gz'),(b/'j_ab_separated/numeric/events.jsonl','j_separated_AB_events.jsonl.gz')]:
 (refs/name).write_bytes(gzip.compress(source.read_bytes()))
for directory,name in [('j_direct','j_original_holdouts.json.gz'),('j_direct_separated','j_separated_holdouts.json.gz')]:
 cfg=load(b/directory/'fit_input.json');train=cfg['train'];snap=[]
 for x,file in zip(cfg['xs'][train:],cfg['files'][train:]):
  v=load(file);snap.append({'x':x,'basis':v['basis'],'reductions':{key:row for key,row in v['reductions'].items() if key.startswith('J')}})
 (refs/name).write_bytes(gzip.compress(json.dumps(snap,separators=(',',':')).encode()))
ranges={}
for variant in ['j_direct','j_direct_separated','j_ab','j_ab_separated']:
 path=b/variant/('numeric/events.jsonl' if 'ab' in variant else 'events.jsonl')
 records=[json.loads(line) for line in path.read_text().splitlines()]
 records=[r for r in records if r['operation'] in ('numeric_IBP','IBP-complete')]
 ranges[variant]={key:[min(r['diagnostics']['linear_algebra'][key] for r in records),max(r['diagnostics']['linear_algebra'][key] for r in records)] for key in ['z_nonzero_relations_submitted','z_independent_relations','z_relation_columns','reduction_coordinate_columns']}
report['all_sampling_point_size_ranges']=ranges
report['provenance_sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'j_integral/src/IBP/build/j-ibp',b/'backend_worker.py',b/'reconstruct.cpp',b/'fi/numeric.wl',root/'j_integral/archive/j_integral_bootstrap/extern/analytic_bootstrap-master/analytic_bootstrap/upstream/pipeline.py']}
(b/'comparison.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k in ['J_original_AB','J_separated_AB','FI_J_full_numeric_verification']},indent=2))
