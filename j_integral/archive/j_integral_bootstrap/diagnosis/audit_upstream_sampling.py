"""Audit unique AB points from saved checks; optionally reconcile worker files.

Reads local receipts only. It does not request IBP, access Tianhe, or validate R.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re


def audit(connection_path,file_list=None,event_path=None):
    raw=Path(connection_path).read_bytes();answer=json.loads(raw);prime=answer['prime']
    center=Fraction(answer['t0']);center=center.numerator*pow(center.denominator,-1,prime)%prime
    checks=answer['checks'];global_check=answer['global_check']
    if not checks or global_check.get('status')!='fresh-full-de-holdouts':
        raise ValueError('completed sector checks and full-DE holdouts required')
    train={int(x)%prime for c in checks for x in c.get('training_points',[])}
    held={int(x)%prime for c in checks for x in c.get('holdout_points',[])}
    final={int(x)%prime for x in global_check['holdout_points']}
    blocks=train|held;normal=blocks|final
    if not final or final&blocks:raise ValueError('final holdouts are not independent of every sector')
    absolute=lambda xs:sorted((x+center)%prime for x in xs)
    result=dict(connection=str(connection_path),connection_sha256=hashlib.sha256(raw).hexdigest(),
        prime=prime,t0=answer['t0'],sector_count=len(checks),
        training_points_tau=sorted(train),training_unique_count=len(train),
        block_holdout_points_tau=sorted(held),block_holdout_unique_count=len(held),
        block_union_points_tau=sorted(blocks),block_union_unique_count=len(blocks),
        block_holdouts_not_used_for_training_tau=sorted(held-train),
        observed_final_holdouts_tau=sorted(final),final_holdout_unique_count=len(final),
        observed_main_points_t=absolute(normal),observed_main_unique_count=len(normal),
        old_driver_expected_extra_prefix_checks=len(checks)-1,
        old_driver_expected_extra_prefix_points=(len(checks)-1)*len(final),
        warning='Observed final holdout coordinates belong to the old run. Removing prefix calls changes the seen set; a corrected run must obtain its own final full-DE checks.')
    if file_list:
        full=[];sub=[]
        for line in Path(file_list).read_text().splitlines():
            line=line.strip()
            if not line:continue
            match=re.search(r'/point(\d+)\.json$',line)
            if not match:raise ValueError('unrecognized numeric output path: '+line)
            record=(int(match.group(1))%prime,line)
            (sub if '/numeric/sectors/' in line else full).append(record)
        main=set(absolute(normal));all_t={t for t,_ in full+sub};extra=all_t-main
        counts=dict(total_successful_worker_outputs=len(full)+len(sub),
            full_outputs=len(full),restricted_outputs=len(sub),
            unique_t_count=len(all_t),main_unique_t_count=len(main),
            full_outputs_at_main_points=sum(t in main for t,_ in full),
            restricted_outputs_at_main_points=sum(t in main for t,_ in sub),
            extra_prefix_full_outputs=sum(t not in main for t,_ in full),
            extra_prefix_restricted_outputs=sum(t not in main for t,_ in sub),
            extra_unique_t_count=len(extra),
            all_main_points_have_full_outputs=main<={t for t,_ in full})
        counts['decomposition_sum']=(counts['full_outputs_at_main_points']+
            counts['restricted_outputs_at_main_points']+counts['extra_prefix_full_outputs']+
            counts['extra_prefix_restricted_outputs'])
        result['worker_output_reconciliation']=counts
        if event_path:
            events=[json.loads(line) for line in Path(event_path).read_text().splitlines() if line.strip()]
            by_prefix={};first=[]
            for event in events:
                task=event.get('task','')
                if not task.startswith('__prefix_de_holdouts__'):continue
                if event.get('route')=='restricted' and 't' in event:t=int(event['t'])%prime
                elif event.get('route')=='full' and 'x' in event:t=(int(event['x'])+center)%prime
                else:continue
                # Event streams include failed attempts: only successful output
                # filenames certify that the expensive point actually exists.
                if t not in all_t:continue
                by_prefix.setdefault(task,set()).add(t)
                if task=='__prefix_de_holdouts__1' and t not in first:first.append(t)
            successful_prefix=set().union(*by_prefix.values()) if by_prefix else set()
            result['prefix_event_reconciliation']=dict(tracked_prefix_task_count=len(by_prefix),
                tracked_successful_prefix_unique_t_count=len(successful_prefix),
                tracked_prefix_points_are_all_extra=successful_prefix<=extra,
                extra_output_point_count_matches_expected_prefix_points=len(extra)==(len(checks)-1)*len(final),
                event_coverage='Partial: restricted/fallback routing events do not enumerate every direct full evaluation.',
                first_prefix_successful_points_t=first,
                first_prefix_points_disjoint_from_block_union=not set(first)&set(absolute(blocks)),
                counterfactual_note='These first-prefix coordinates are candidates for the corrected final draw, not a replacement for executing and checking that full-DE draw.')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('connection');parser.add_argument('--numeric-file-list');parser.add_argument('--events')
    parser.add_argument('--output',required=True);args=parser.parse_args()
    result=audit(args.connection,args.numeric_file_list,args.events)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['sector_count','training_unique_count','block_union_unique_count',
        'final_holdout_unique_count','observed_main_unique_count','old_driver_expected_extra_prefix_points']},indent=2))
    if 'worker_output_reconciliation'in result:print(json.dumps(result['worker_output_reconciliation'],indent=2))
    if 'prefix_event_reconciliation'in result:print(json.dumps(result['prefix_event_reconciliation'],indent=2))


if __name__=='__main__':main()
