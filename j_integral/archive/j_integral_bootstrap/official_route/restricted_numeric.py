"""Support-restricted J numeric adapter; the upstream package is unchanged."""
import json,re,fcntl,subprocess
from pathlib import Path
from time import perf_counter
from upstream_j import NumericBackend
from analytic_bootstrap.storage import digest


def indices(name):
    m=re.fullmatch(r'(J|Bz[01])\{([0-9,]+)\}',name)
    if m is None:raise ValueError('unsupported master label: '+str(name))
    return m[1],tuple(map(int,m[2].split(',')))


def label(kind,nu):return kind+'{'+','.join(map(str,nu))+'}'


class RestrictedUnavailable(Exception):pass


def closure_tasks(bundle):
    """Extract names from learned structural lower_indices only; no C values."""
    record=bundle._record if hasattr(bundle,'_record') else bundle
    sectors=record['definition']['sectors'];masters=[m for sector in sectors for m in sector['masters']]
    result={}
    for sector in sectors:
        lower=record['de_templates'].get(sector['id'],{}).get('lower_indices',[])
        result[sector['id']]=[masters[i] for i in lower]+list(sector['masters'])
    return result


class SectorRestrictedNumericBackend(NumericBackend):
    """Embed a closed three-branch subfamily back into all requested columns.

    Independence follows from taking a subset of the fixed independent full
    basis and exact physical label embedding. Completeness is tested by the
    native worker on every subfamily point, including all retained derivatives.
    Only recognized incompleteness/degeneracy triggers a recorded full fallback;
    malformed data and unrelated executable failures are never swallowed.
    """
    def __init__(self,*args,closure_by_task=None,**kwargs):
        super().__init__(*args,**kwargs)
        self.closure_by_task=closure_by_task
        self.restricted={};self.events=[]

    def event(self,record):
        self.events.append(record)
        with (self.work/'restricted_events.jsonl').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX);f.write(json.dumps(record)+'\n')

    def plan(self,req):
        if self.cfg.get('cut') or self.cfg.get('cuts'):raise RestrictedUnavailable('propagator cuts are unsupported')
        active=set();wanted=[];n=len(self.cfg['propagators'])
        for target in req['targets']:
            if target.get('kind')!='derivative' or target.get('variable')!='tau':raise ValueError('unsupported target kind')
            m=target['master'];_,nu=indices(m)
            if len(nu)!=n:raise ValueError('master index length differs from family')
            active.update(i for i,v in enumerate(nu) if v);wanted.append(m)
        closure=None
        if self.closure_by_task is not None and req.get('task') in self.closure_by_task:
            closure=set(self.closure_by_task[req['task']])
            if not closure<=set(self.names) or not set(wanted)<=closure:raise ValueError('invalid structural task closure')
            for name in closure:
                _,nu=indices(name);active.update(i for i,v in enumerate(nu) if v)
        keep=sorted(active)
        if len(keep)==n:raise RestrictedUnavailable('target support is the full family')
        if {self.cfg['propagators'][i]['branch'] for i in keep}!={0,1,2}:raise RestrictedUnavailable('restriction does not preserve all three branches')
        names=[];small=[]
        for name in self.names:
            kind,nu=indices(name)
            if len(nu)!=n:raise ValueError('basis index length differs from family')
            if (closure is None or name in closure) and all(not v or i in active for i,v in enumerate(nu)):
                names.append(name);small.append(label(kind,[nu[i] for i in keep]))
        if not set(wanted)<=set(names):raise RestrictedUnavailable('requested target is not in the retained basis')
        if len(set(small))!=len(small):raise ValueError('restriction mapping is not injective')
        key=digest(keep) if closure is None else digest([keep,names]);folder=self.work/'sectors'/key;folder.mkdir(parents=True,exist_ok=True)
        if key not in self.restricted:
            cfg=dict(self.cfg);cfg['propagators']=[self.cfg['propagators'][i] for i in keep]
            parents=sorted({tuple(indices(m)[1]) for m in small})
            cfg['targets']=[list(nu) for nu in parents]
            # Index-sensitive cut/rank layouts are not inferred or silently changed.
            for field in ['cut','cuts','top_sector','numerators','irreducible_numerators']:
                if self.cfg.get(field):raise RestrictedUnavailable('index-dependent configuration '+field)
            for filename,obj in [('family.json',cfg),('basis.json',small),('targets.json',parents),('mapping.json',dict(active_positions=keep,original_masters=names,compressed_masters=small))]:
                path=folder/filename
                if path.exists() and json.loads(path.read_text())!=json.loads(json.dumps(obj)):raise ValueError('restriction identity changed')
                if not path.exists():path.write_text(json.dumps(obj,indent=2)+'\n')
            backend=NumericBackend(folder/'family.json',folder/'basis.json',folder/'targets.json',self.binary,folder/'numeric',gauge=self.gauge)
            self.restricted[key]=(backend,dict(zip(names,small)),dict(zip(small,names)),keep)
        return self.restricted[key]

    def prepare(self,req):
        if req['mode']!='de' or req['prime']!=self.prime or set(req['masters'])!=set(self.names):raise ValueError('unsupported numeric request')
        full=super().prepare(req)
        try:small,forward,backward,keep=self.plan(req)
        except RestrictedUnavailable as e:
            self.event(dict(route='full',stage='plan',reason=str(e),task=req.get('task')));return full
        owner=self
        class Session:
            def evaluate(self,x):
                start=perf_counter();before=small.calls
                try:
                    try:d=small.point(x)
                    except subprocess.CalledProcessError:
                        t=(x+small.center)%small.prime;log=(small.work/f'point{t}.log').read_text()
                        if any(s in log for s in ['dependent resume basis','chosen J basis does not span target/derivative']):raise RestrictedUnavailable(log[-1500:])
                        raise
                    except ValueError as e:
                        # NumericBackend has one aggregate validation exception;
                        # examine native closure only, never hide other mismatches.
                        t=(x+small.center)%small.prime;path=small.work/f'point{t}.json'
                        if path.exists():
                            invalid=json.loads(path.read_text())
                            if invalid.get('prime')==small.prime and int(invalid.get('t_point',-1))==t and invalid.get('masters_at_point')==small.names and invalid.get('unclosed_derivatives'):
                                raise RestrictedUnavailable('retained derivatives do not close') from e
                        raise
                    t=(x+owner.center)%owner.prime
                    if owner.gauge and t in (0,1):
                        from analytic_bootstrap.upstream import BadSamplePoint
                        raise BadSamplePoint('normalization endpoint')
                    shift=owner.gauge*(pow(t,-1,owner.prime)+pow((t-1)%owner.prime,-1,owner.prime))%owner.prime if owner.gauge else 0
                    rows=[]
                    for target in req['targets']:
                        name=target['master'];row=d['DE_values'][forward[name]]
                        embedded={backward[m]:v for m,v in row.items()}
                        rows.append([(embedded.get(m,0)+(shift if m==name else 0))%owner.prime for m in req['masters']])
                    owner.event(dict(route='restricted',task=req.get('task'),t=t,active_positions=keep,rank=len(small.names),FBI_rank=d['FBI_rank'],seconds=perf_counter()-start))
                    return dict(request_sha256=digest(req),x=x,values=rows)
                except RestrictedUnavailable as e:
                    owner.event(dict(route='full',stage='point',task=req.get('task'),x=x,reason=str(e)))
                    return full.evaluate(x)
                finally:owner.calls+=small.calls-before
            def close(self):full.close()
        return Session()
