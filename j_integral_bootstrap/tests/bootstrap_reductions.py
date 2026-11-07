"""Sectorwise J reduction using a closed known DE, pole witnesses and bootstrap."""
import argparse
from collections import defaultdict
import fcntl
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from time import perf_counter

from reference_space import BootstrapSpec,RatFunc,RatMat,build_pl_ansatz
from analytic_bootstrap.sampling import InconsistentSampleError
from rational_learning import make_learner
from fixed_fit import fixed_entry_fitter
from native_fit import joint_fit
from analytic_bootstrap.upstream.assembly import assemble_matrix
from sector_de import components


def read(path):return json.loads(Path(path).read_text())
def mask(name):return tuple(v>0 for v in json.loads('['+name[name.index('{')+1:name.index('}')]+']'))
def record(r):return dict(numerator_tau=r.num,denominator_tau=r.den)


def main(args):
    work=Path(args.work);work.mkdir(parents=True,exist_ok=True);start=perf_counter()
    connection=read(args.connection);p=connection['prime'];masters=connection['masters'];index={m:i for i,m in enumerate(masters)}
    zero=RatFunc.make([0],[1],p);center=Fraction(connection['t0']);center=center.numerator*pow(center.denominator,-1,p)%p
    L=make_learner(p,center);one,ep=L.one,L.ep
    known_masters=set(connection.get('known_rows',masters));available={};graph={i:set() for i in range(len(masters))}
    for directory in args.points:
        for path in Path(directory).glob('point*.json'):
            if '.tmp.' in path.name:continue
            d=read(path)
            if d['prime']!=p or d['masters_at_point']!=masters or d['unclosed_derivatives']:raise ValueError('point identity/closure differs')
            t=int(d['t_point']);rows=available.setdefault(t,{})
            for name,row in d['reduction_values'].items():
                if name in rows and rows[name]!=row:raise ValueError('reduction providers disagree')
                rows[name]=row
            for m,row in d.get('DE_values',{}).items():graph[index[m]].update(index[n] for n in row if n!=m)
    targets=sorted({name for rows in available.values() for name in rows});support={name:set() for name in targets}
    for rows in available.values():
        for name,row in rows.items():support[name].update(index[m] for m,v in row.items() if v)
    comp={i:k for k,c in enumerate(components(graph)) for i in c}
    def rf(v):return RatFunc.make(v['numerator_tau'],v['denominator_tau'],p)
    D={(index[m],index[n]):rf(v) for m,row in connection['de'].items() if m in known_masters
       for n,v in row['rhs'].items()}
    learned={}
    for path in args.partial:
        data=read(path)
        if data['prime']!=p or Fraction(data['t0'])!=Fraction(connection['t0']):raise ValueError('partial field/coordinate differs')
        for name,row in data['reductions'].items():
            for m,v in row.items():
                key=(name,index[m]);value=rf(v)
                if key in learned and not (learned[key]-value).is_zero():raise ValueError('partial reconstructions disagree')
                learned[key]=value
    done=set(masters);reductions={};reports=[];used_rows={};checkpoint=work/'reduction.json'
    if checkpoint.exists():
        saved=read(checkpoint)
        if (saved['prime']!=p or saved['masters']!=masters or 't0' not in saved
                or Fraction(saved['t0'])!=Fraction(connection['t0'])):
            raise ValueError('checkpoint field/basis/coordinate differs')
        done=set(saved['known_targets']);reports=saved['groups'];used_rows=saved.get('connection_rows_used',{})
        for m,digest in used_rows.items():
            now=hashlib.sha256(json.dumps(connection['de'][m]['rhs'],sort_keys=True).encode()).hexdigest()
            if now!=digest:raise ValueError('previously used connection row changed')
        for name,row in saved['reductions'].items():reductions[name]={index[m]:rf(v) for m,v in row.items()}
    def save(status,**extra):
        result=dict(status=status,prime=p,t0=connection['t0'],masters=masters,master_count=len(masters),
            known_targets=sorted(done),target_count=len(targets),known_target_count=len(done),
            tested_objects=targets,zero_objects=sorted(name for name,row in reductions.items() if not row),
            reductions={name:{masters[j]:record(r) for j,r in row.items() if not r.is_zero()}
                        for name,row in reductions.items()},groups=reports,connection_source=args.connection,
            connection_rows_used=used_rows,total_seconds=perf_counter()-start,
            DE_complete=connection.get('status')=='complete',**extra)
        if status=='complete' and result['DE_complete']:result['de']=connection['de']
        temp=checkpoint.with_suffix('.tmp.json');temp.write_text(json.dumps(result,separators=(',',':'))+'\n');temp.replace(checkpoint)
    for name in targets:
        if name in done:continue
        if all((name,j) in learned for j in support[name]):
            reductions[name]={j:learned[name,j] for j in support[name]};done.add(name)
    save('running')
    groups=defaultdict(list)
    for name in targets:
        if name not in done:groups[name[0],mask(name)].append(name)
    skipped=[];failed=[]
    for (kind,sector),names in sorted(groups.items(),key=lambda item:(sum(item[0][1]),item[0])):
        deps=set().union(*(support[name] for name in names))
        frontier=list(deps)
        while frontier:
            i=frontier.pop()
            for j in graph[i]:
                if j not in deps:deps.add(j);frontier.append(j)
        sub=sorted(deps)
        if any(masters[i] not in known_masters for i in sub):skipped.extend(names);continue
        pos={v:i for i,v in enumerate(sub)}
        B={(pos[i],pos[j]):r for (i,j),r in D.items() if i in pos and j in pos}
        source=one
        source_rows={index[m] for m in masters if m[0]==kind and mask(m)==sector}
        for (i,j),r in D.items():
            if i in source_rows and comp[i]==comp[j]:source=L.lcm(source,L.poly(r.den))
        colq={j:one for j in sub}
        fullcol={j:one for j in sub};fullsource=source
        for (i,j),r in D.items():
            if i in pos and j in pos and comp[i]==comp[j]:colq[j]=L.lcm(colq[j],L.poly(r.den))
            if i in pos and j in pos:fullcol[j]=L.lcm(fullcol[j],L.poly(r.den))
            if i in source_rows and j in pos:fullsource=L.lcm(fullsource,L.poly(r.den))
        anchors={name:one for name in names};group_anchor=one
        for name in names:
            for j in support[name]:
                if (name,j) in learned:
                    anchors[name]=L.lcm(anchors[name],L.poly(learned[name,j].den))
                    group_anchor=L.lcm(group_anchor,L.poly(learned[name,j].den))
        print(json.dumps(dict(stage='sector',sector=sector,kind=kind,targets=len(names),sub_rank=len(sub),
                              completed=len(done),total=len(targets))),flush=True)
        expensive_attempts={name:0 for name in names}
        # Start from simple coefficient columns to obtain useful pole witnesses.
        for passno in range(2):
            for name in names:
                times=sorted(t for t,rows in available.items() if name in rows)
                for j in sorted(support[name],key=lambda k:(int(colq[k].degree()),k)):
                    if (name,j) in learned:continue
                    values=[available[t][name].get(masters[j],0) for t in times]
                    a=anchors[name];candidates=[a,source,L.lcm(a,source),L.lcm(a,colq[j]),
                        L.lcm(L.lcm(a,source),colq[j]),a*colq[j]*ep,source*colq[j]*ep,group_anchor,
                        L.monic(fullsource.gcd(fullcol[j])),L.lcm(a,fullcol[j]),fullsource]
                    fitted=None;seen=set()
                    for den in candidates:
                        key=tuple(int(c) for c in L.monic(den))
                        if key in seen:continue
                        seen.add(key);fitted=L.fit(den,times,values)
                        if fitted is not None:break
                    if fitted is None and passno and expensive_attempts[name]<args.witness_budget:
                        expensive_attempts[name]+=1
                        for den in [L.lcm(L.lcm(a,source),colq[j]),a*colq[j]*ep,
                                    L.monic(fullsource.gcd(fullcol[j]))]:
                            fitted=L.prune(den,times,values,budget=32)
                            if fitted is None:fitted=L.boost(den,times,values)
                            if fitted is not None:break
                    if fitted is not None:
                        learned[name,j]=fitted;anchors[name]=L.lcm(a,L.poly(fitted.den))
                        group_anchor=L.lcm(group_anchor,L.poly(fitted.den))
        pending=[]
        for name in names:
            if all((name,j) in learned for j in support[name]):
                reductions[name]={j:learned[name,j] for j in support[name]};done.add(name)
            else:pending.append(name)
        if not pending:
            reports.append(dict(sector=sector,kind=kind,targets=len(names),route='pole-witness-learning'))
            save('running');continue
        times=sorted(t for t,rows in available.items() if all(name in rows for name in pending))
        train=times[:-args.internal_holdouts];holdout=times[-args.internal_holdouts:];attempts=[];solved=False
        for power,growth,gauge in [(0,0,1),(0,2,2),(0,4,3),(1,2,3),(1,4,3),(2,4,3)]:
            ds=[];caps=[]
            for name in pending:
                denrow=[];caprow=[]
                for j in sub:
                    if j not in support[name]:denrow.append([1]);caprow.append(-1)
                    elif (name,j) in learned:
                        r=learned[name,j];denrow.append(r.den);caprow.append(len(r.num)-1)
                    else:
                        base=L.lcm(L.lcm(anchors[name],source),L.lcm(colq[j],ep))
                        den=base if power==0 else (base*colq[j])**power
                        denrow.append([int(c) for c in L.monic(den)]);caprow.append(int(den.degree())+growth)
                ds.append(denrow);caps.append(caprow)
            h=ep**gauge;shift=RatFunc.make([int(c) for c in h.derivative()],[int(c) for c in h],p)
            entries=dict(B)
            for i in range(len(sub)):entries[i,i]=entries.get((i,i),zero)+shift
            spec=BootstrapSpec.from_inputs(family={'source':'full-J-reduction'},variable='tau',assignments={},
                target=pending,masters=[masters[j] for j in sub],de=RatMat(entries,(len(sub),len(sub)),p),
                denominators=ds,numerator_caps=caps,prime=p,mode='reduction')
            print(json.dumps(dict(stage='construct-reduction',targets=len(pending),power=power,growth=growth,
                                  ambient=sum(c+1 for row in caps for c in row))),flush=True)
            ansatz=build_pl_ansatz(spec,backend='PL-Hasse',laurent_mode='structural-refined-all')
            fixed={(a,b):learned[name,j] for a,name in enumerate(pending) for b,j in enumerate(sub) if (name,j) in learned}
            def sample(t):return dict(prime=p,x=(t-center)%p,
                values=[[available[t][name].get(masters[j],0) for j in sub] for name in pending])
            try:
                fit=joint_fit(ansatz,[sample(t) for t in train],fixed,[(t-center)%p for t in times])
                if not fit['compatible']:raise InconsistentSampleError('joint equations inconsistent')
                if fit['parameters'] is None:
                    attempts.append(dict(power=power,growth=growth,gauge=gauge,status='need_points',
                        dimension=ansatz.dimension,rank=fit['rank']));continue
                for t in holdout:
                    if ansatz.evaluate((t-center)%p,fit['parameters'])!=sample(t)['values']:
                        raise ValueError('holdout disagrees with fitted ansatz')
                answer=assemble_matrix(ansatz,fit['parameters'])
                for a,name in enumerate(pending):
                    reductions[name]={j:answer.entries[a,b] for b,j in enumerate(sub) if (a,b) in answer.entries and not answer.entries[a,b].is_zero()};done.add(name)
                reports.append(dict(sector=sector,kind=kind,targets=len(names),bootstrap_targets=len(pending),
                    route='analyticity-bootstrap',ambient=sum(c+1 for row in caps for c in row),
                    dimension=ansatz.dimension,joint_fit={k:v for k,v in fit.items() if k!='parameters'},training_points=len(train),
                    holdout_points=len(holdout),power=power,growth=growth,gauge=gauge,attempts=attempts))
                for j in sub:
                    m=masters[j];used_rows[m]=hashlib.sha256(json.dumps(connection['de'][m]['rhs'],sort_keys=True).encode()).hexdigest()
                solved=True;save('running');break
            except ValueError as e:
                if not isinstance(e,InconsistentSampleError) and str(e)!='holdout disagrees with fitted ansatz':raise
                attempts.append(dict(power=power,growth=growth,gauge=gauge,status='inconsistent',error=str(e)))
        if not solved:
            save('need_points_or_envelope',pending_targets=pending,attempts=attempts)
            print(json.dumps(dict(status='need_points_or_envelope',known_targets=len(done),pending=pending)),flush=True)
            failed.append(dict(targets=pending,attempts=attempts));continue
    missing_targets=[name for name in targets if name not in done]
    if missing_targets:
        save('need_points_or_envelope' if failed else 'waiting_for_connection',
             pending_targets=missing_targets,failed_groups=failed,waiting_for_connection=skipped)
        return 2
    save('complete',all_reconstructed=True)
    print(json.dumps(dict(status='complete',target_count=len(targets),nonmaster_targets=len(reductions))),flush=True)
    return 0


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('work');p.add_argument('--connection',required=True)
    p.add_argument('--points',action='append',required=True);p.add_argument('--partial',action='append',default=[])
    p.add_argument('--witness-budget',type=int,default=3)
    p.add_argument('--internal-holdouts',type=int,default=3)
    args=p.parse_args();Path(args.work).mkdir(parents=True,exist_ok=True)
    with (Path(args.work)/'run.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);raise SystemExit(main(args))
