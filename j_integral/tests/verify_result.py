"""Verify the original AB result and summarize its own checks/timings plus IBP receipts."""
import argparse,json,time,gzip
from pathlib import Path
from analytic_bootstrap.upstream import SolvedDE
from analytic_bootstrap.upstream.model import read_matrix
p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--run',required=True);a=p.parse_args();case=Path(a.case).resolve();run=Path(a.run).resolve()
load=lambda path:json.loads(path.read_text())
result=SolvedDE.from_dict(load(run/'solution.json'));anchor=json.loads(gzip.decompress((case/'input/anchor.json.gz').read_bytes()));expected=read_matrix(anchor['de'],anchor['prime'],(len(anchor['masters']),)*2);assert result.context['masters']==anchor['masters']
zero=None;keys=set(result.matrix.entries)|set(expected.entries)
for key in keys:
 left=result.matrix.entries.get(key);right=expected.entries.get(key)
 if left is None:assert right.is_zero(),key
 elif right is None:assert left.is_zero(),key
 else:assert (left-right).is_zero(),key
receipts=[load(path) for path in (run/'numeric/receipts').glob('point*.json')];assert receipts and all(r['returncode']==0 and r['validated'] for r in receipts)
training=set().union(*(set(c['training_points']) for c in result.checks));sector_checks=set().union(*(set(c['holdout_points']) for c in result.checks));final=set(result.global_check['holdout_points']);needed=training|sector_checks|final
assert len(receipts)==len(needed),'IBP output count must equal distinct actually required points'
events=[json.loads(s) for s in (run/'numeric/events.jsonl').read_text().splitlines()];assert not any(e.get('task','').startswith('__prefix_de_holdouts__') for e in events)
report=dict(status='complete',all_rational_entries_equal=True,masters=len(anchor['masters']),blocks=len(result.checks),matrix_entries_including_zeros=len(anchor['masters'])**2,explicit_union_entries=len(keys),training_distinct_points=len(training),training_plus_sector_checks_distinct_points=len(training|sector_checks),final_holdout_points=len(final),total_distinct_points=len(needed),actual_IBP_invocations=len(receipts),failed_IBP_invocations=sum(e['operation']=='IBP-failed' for e in events),IBP_seconds_sum=sum(r['wall_seconds'] for r in receipts),IBP_seconds_min=min(r['wall_seconds'] for r in receipts),IBP_seconds_max=max(r['wall_seconds'] for r in receipts),learn=load(run/'learn_timings.json'),solve=load(run/'solve_timings.json'),wall_seconds=time.time()-load(run/'run.json')['start'],sampling=load(run/'backend.json')['sampling'])
(run/'verification.json').write_text(json.dumps(report,indent=2)+'\n');metadata=load(run/'run.json');metadata.update(end=time.time(),scientific_complete=True);(run/'run.json').write_text(json.dumps(metadata,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
