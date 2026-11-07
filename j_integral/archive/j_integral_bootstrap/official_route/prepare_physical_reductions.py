"""Prepare R-only physical restrictions without enlarging support for DE aliases."""
import argparse,json,subprocess
from pathlib import Path
from restricted_numeric import indices,label
p=argparse.ArgumentParser();p.add_argument('--family',required=True);p.add_argument('--checkpoint',required=True);p.add_argument('--reference-point',required=True);p.add_argument('--binary',required=True);p.add_argument('--keep',required=True);p.add_argument('--exclude-file');p.add_argument('work');a=p.parse_args()
read=lambda path:json.loads(Path(path).read_text())
f=read(a.family);S=read(a.checkpoint);V=read(a.reference_point);active=list(map(int,a.keep.split(',')));n=len(f['propagators']);done=set(S['known_targets']);exclude=set(read(a.exclude_file)) if a.exclude_file else set()
targets=[name for name in V['reduction_values'] if name.startswith('J') and name not in done|exclude and all(not v or i in active for i,v in enumerate(indices(name)[1]))]
if not targets:raise ValueError('no pending targets')
O=Path(a.work).resolve();O.mkdir(parents=True,exist_ok=True);(O/'points').mkdir(exist_ok=True)
compact=lambda name:label(indices(name)[0],[indices(name)[1][i] for i in active])
f['propagators']=[f['propagators'][i] for i in active];parents=targets+[name for name in V['masters_at_point'] if name.startswith('J') and all(not v or i in active for i,v in enumerate(indices(name)[1]))]
f['targets']=sorted({tuple(indices(compact(name))[1]) for name in parents})
f['targets']=[list(nu) for nu in f['targets']]
for filename,obj in [('family.json',f),('targets.json',f['targets']),('no_derivatives.json',[]),('original_targets.json',targets)]: (O/filename).write_text(json.dumps(obj))
out=O/'points'/('point'+str(V['t_point'])+'.json')
cmd=[str(Path(a.binary).resolve()),str(O/'family.json'),str(out),'--point',str(V['t_point']),'--targets-file',str(O/'targets.json'),'--derivative-names-file',str(O/'no_derivatives.json'),'--integration','atomic','--values','1','--compact','1']
with (O/'probe.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True)
local=read(out)
def expand(name):
 kind,nu=indices(name);z=[0]*n
 for i,v in zip(active,nu):z[i]=v
 return label(kind,z)
original=[expand(name) for name in local['masters_at_point']]
# Keep the transformation explicit; a nontrivial projection must be separately verified.
if not set(original)<=set(V['masters_at_point']):raise ValueError(('natural basis needs projection',sorted(set(original)-set(V['masters_at_point']))))
mp=dict(zip(local['masters_at_point'],original))
for name in targets:
 if {mp[c]:v for c,v in local['reduction_values'][compact(name)].items()}!=V['reduction_values'][name]:raise ValueError('reference reduction differs: '+name)
mapping=dict(active_positions=active,original_size=n,original_masters=original,compressed_masters=local['masters_at_point'])
(O/'mapping.json').write_text(json.dumps(mapping));(O/'basis.json').write_text(json.dumps(local['masters_at_point']));print(json.dumps(dict(targets=len(targets),rank=len(original),active_positions=active)))
