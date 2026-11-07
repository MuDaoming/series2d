"""Topology-independent physical-index compression for J numerical workers."""
import json
from pathlib import Path
from restricted_numeric import indices,label


def make_restriction(family, masters, targets, *, active_positions=None):
    """Return compressed family/basis/parent targets and reversible labels.

    `masters` is an already chosen closed subspace, not inferred from kinematics.
    targets are original-index vectors. No C coefficients enter this function.
    """
    n=len(family['propagators']);targets=[tuple(t) for t in targets]
    for t in targets:
        if len(t)!=n or any(v<0 for v in t):raise ValueError('unsupported target index')
    active={i for t in targets for i,v in enumerate(t) if v}
    for m in masters:
        _,nu=indices(m)
        if len(nu)!=n:raise ValueError('basis length differs')
        active.update(i for i,v in enumerate(nu) if v)
    if active_positions is not None:
        if not active<=set(active_positions):raise ValueError('active positions omit indices')
        active=set(active_positions)
    active=sorted(active)
    if {family['propagators'][i]['branch'] for i in active}!={0,1,2}:raise ValueError('three branches must be preserved')
    def compress(name):
        kind,nu=indices(name)
        if any(v and i not in active for i,v in enumerate(nu)):raise ValueError('label outside restriction')
        return label(kind,[nu[i] for i in active])
    cfg=dict(family);cfg['propagators']=[family['propagators'][i] for i in active]
    parents={tuple(t[i] for i in active) for t in targets}
    parents.update(tuple(indices(compress(m))[1]) for m in masters)
    cfg['targets']=[list(t) for t in sorted(parents)]
    return dict(family=cfg,basis=[compress(m) for m in masters],targets=cfg['targets'],mapping=dict(original_size=n,active_positions=active,original_masters=list(masters),compressed_masters=[compress(m) for m in masters]))


def expand_name(name,mapping):
    kind,nu=indices(name);active=mapping['active_positions'];out=[0]*mapping['original_size']
    if len(nu)!=len(active):raise ValueError('compressed label length differs')
    for i,v in zip(active,nu):out[i]=v
    return label(kind,out)


def write_restriction(directory,packet):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    for key in ['family','basis','targets','mapping']:
        (directory/(key+'.json')).write_text(json.dumps(packet[key],indent=2)+'\n')
