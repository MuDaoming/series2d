"""Bootstrap a full J connection from independently learned diagonals and IBP points.

Unknown coupling envelopes are adaptive candidates, never imported from a full
reference connection. Every constructed block is checked at held-out points.
"""
import argparse
import fcntl
from fractions import Fraction
import json
from pathlib import Path
from time import perf_counter

from reference_space import BootstrapSpec,RatFunc,RatMat,SampleFitter,build_pl_ansatz,build_de_ansatz
from sector_de import components
from analytic_bootstrap.sampling import InconsistentSampleError
from rational_learning import make_learner
from fixed_fit import fixed_entry_fitter
from native_fit import joint_fit
from analytic_bootstrap.upstream.assembly import assemble_matrix


def read(path):return json.loads(Path(path).read_text())


def main(args):
    work=Path(args.work);work.mkdir(parents=True,exist_ok=True)
    cfg=read(args.family);p=cfg['prime'];masters=read(args.basis);n=len(masters)
    index={m:i for i,m in enumerate(masters)};zero=RatFunc.make([0],[1],p)
    center=Fraction(cfg['t0']);center=center.numerator*pow(center.denominator,-1,p)%p
    all_points={};full_points={};graph={i:set() for i in range(n)}
    for directory in args.points:
        for path in sorted(Path(directory).glob('point*.json')):
            if '.tmp.' in path.name:continue
            d=read(path)
            if d['prime']!=p or d['masters_at_point']!=masters or d['unclosed_derivatives']:
                raise ValueError('point field, basis or closure differs')
            t=int(d['t_point']);target=all_points.setdefault(t,{})
            for m,row in d['DE_values'].items():
                if m in target and target[m]!=row:raise ValueError('point providers disagree')
                target[m]=row
                graph[index[m]].update(index[j] for j,v in row.items() if v and j!=m)
            if set(d['DE_values'])==set(masters):full_points[t]=d
    if len(full_points)<9:raise ValueError('need at least nine complete pilot points')
    groups=components(graph)
    component={i:k for k,c in enumerate(groups) for i in c}
    def rf(v):return RatFunc.make(v['numerator_tau'],v['denominator_tau'],p)
    diagonal={};learned={}
    for path in args.diagonals:
        d=read(path)
        if d['prime']!=p or 't0' not in d or Fraction(d['t0'])!=Fraction(cfg['t0']):
            raise ValueError('diagonal field/coordinate differs')
        for m,row in d['de'].items():
            for j,v in row['rhs'].items():
                key=(index[m],index[j]);value=rf(v)
                if key in diagonal and not (diagonal[key]-value).is_zero():raise ValueError('diagonal sources disagree')
                diagonal[key]=value
    for path in args.partial:
        d=read(path)
        if d['prime']!=p or Fraction(d['t0'])!=Fraction(cfg['t0']):raise ValueError('partial reconstruction field/coordinate differs')
        for m,row in d['de'].items():
            for j,v in row['rhs'].items():
                key=(index[m],index[j]);value=rf(v)
                if key in learned and not (learned[key]-value).is_zero():raise ValueError('partial reconstructions disagree')
                learned[key]=value
    tools=make_learner(p,center)
    poly,one,ep,monic,lcm=tools.poly,tools.one,tools.ep,tools.monic,tools.lcm
    known_denominator_fit,prune_denominator,boost_denominator=tools.fit,tools.prune,tools.boost
    built={};known=[];reports=[];begin=perf_counter()
    checkpoint=work/'connection.json'
    def save(status,**extra):
        de={m:{'rhs':{}} for m in masters}
        for (i,j),r in built.items():
            if not r.is_zero():de[masters[i]]['rhs'][masters[j]]={'numerator_tau':r.num,'denominator_tau':r.den}
        result=dict(status=status,prime=p,t0=cfg['t0'],masters=masters,master_count=n,
                    known_rows=[masters[i] for i in known],block_count=len(groups),blocks=reports,
                    de=de,total_seconds=perf_counter()-begin,
                    envelope_source='candidate from known diagonal and lower-column poles',**extra)
        temp=checkpoint.with_suffix('.tmp.json');temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(checkpoint)
    saved=read(checkpoint) if checkpoint.exists() and not args.restart else None
    if saved:
        if (saved['prime']!=p or saved['masters']!=masters or 't0' not in saved
                or Fraction(saved['t0'])!=Fraction(cfg['t0'])):
            raise ValueError('checkpoint field/basis/coordinate differs')
        known=[index[m] for m in saved['known_rows']];reports=saved['blocks']
        for m,row in saved['de'].items():
            for j,v in row['rhs'].items():built[index[m],index[j]]=rf(v)
    unresolved=[]
    for block in groups:
        if all(i in known for i in block):continue
        if any(i in known for i in block):raise ValueError('checkpoint cuts through a block')
        missing=[]
        for i in block:
            for j in block:
                if j in graph[i] or i==j and any(all_points[t].get(masters[i],{}).get(masters[j],0) for t in all_points):
                    if (i,j) not in diagonal:missing.append([masters[i],masters[j]])
                    else:built[i,j]=diagonal[i,j]
        if missing:
            save('need_diagonal',missing_diagonal=missing)
            print(json.dumps(dict(status='need_diagonal',block=[masters[i] for i in block],missing=len(missing))),flush=True)
            unresolved.append(dict(block=[masters[i] for i in block],reason='need_diagonal',missing_diagonal=missing));continue
        deps=set()
        def follow(i):
            for j in graph[i]:
                if j not in block and j not in deps:deps.add(j);follow(j)
        for i in block:follow(i)
        sub=[i for i in known if i in deps]
        if any(i not in known for i in deps):
            unresolved.append(dict(block=[masters[i] for i in block],reason='waiting_for_dependencies'));continue
        active={(i,j) for i in block for j in sub if j in graph[i]}
        print(json.dumps(dict(stage='block',block=[masters[i] for i in block],sub_rank=len(sub),
                              active_entries=len(active),unknown_entries=sum(k not in learned for k in active))),flush=True)
        if not active:
            known.extend(block);reports.append(dict(top=[masters[i] for i in block],route='diagonal-only'))
            save('running');continue
        if all(k in learned for k in active):
            for k in active:built[k]=learned[k]
            known.extend(block);reports.append(dict(top=[masters[i] for i in block],route='pilot-rational-reconstruction'))
            save('running');continue
        qs=one
        for i in block:
            for j in block:
                if (i,j) in built:qs=lcm(qs,poly(built[i,j].den))
        # A unipotent change of master basis can put a pole solely in C.
        # Already reconstructed entries of this same block are extra witnesses.
        for k in active:
            if k in learned:qs=lcm(qs,poly(learned[k].den))
        cols={j:one for j in sub}
        diagonal_cols={j:one for j in sub}
        for i in sub:
            for j in sub:
                if (i,j) in built and not built[i,j].is_zero():
                    cols[j]=lcm(cols[j],poly(built[i,j].den))
                    if component[i]==component[j]:diagonal_cols[j]=lcm(diagonal_cols[j],poly(built[i,j].den))
        subpos={v:i for i,v in enumerate(sub)};toppos={v:i for i,v in enumerate(block)}
        baseB={(subpos[i],subpos[j]):r for (i,j),r in built.items() if i in subpos and j in subpos and not r.is_zero()}
        baseA={(toppos[i],toppos[j]):built[i,j] for i in block for j in block if (i,j) in built and not built[i,j].is_zero()}
        times=sorted(t for t,rows in all_points.items() if all(masters[i] in rows for i in block))
        if len(times)<9:raise ValueError('too few row-complete points')
        # Learn exact entry envelopes by fitting a numerator over known pole
        # candidates. This uses all pilot data; independent global validation
        # therefore remains mandatory after the complete connection is built.
        rowqs={i:one for i in block};bare_rows={i:one for i in block}
        for i in block:
            for j in block:
                if (i,j) in built:
                    rowqs[i]=lcm(rowqs[i],poly(built[i,j].den));bare_rows[i]=rowqs[i]
            for k in active:
                if k[0]==i and k in learned:rowqs[i]=lcm(rowqs[i],poly(learned[k].den))
        learned_now=0;expensive_attempts={i:0 for i in block}
        learning_order=sorted(active,key=lambda k:(int(diagonal_cols[k[1]].degree()),k[0],k[1]))
        for i,j in learning_order+learning_order:
            if (i,j) in learned:continue
            witness_allowed=expensive_attempts[i]<args.witness_budget
            values=[all_points[t][masters[i]].get(masters[j],0) for t in times]
            candidates=[lcm(lcm(bare_rows[i],diagonal_cols[j]),ep),
                        lcm(lcm(rowqs[i],diagonal_cols[j]),ep),rowqs[i]*diagonal_cols[j]*ep,
                        lcm(lcm(rowqs[i],cols[j]),ep),rowqs[i]*cols[j]*ep,
                        qs*cols[j]*ep,lcm(lcm(rowqs[i]**2,cols[j]**2),ep**2)]
            seen=set()
            for den in candidates:
                den=monic(den);key=tuple(int(c) for c in den)
                if key in seen:continue
                seen.add(key)
                fitted=known_denominator_fit(den,times,values)
                if fitted is not None:
                    learned[i,j]=fitted;learned_now+=1
                    rowqs[i]=lcm(rowqs[i],poly(fitted.den));qs=lcm(qs,poly(fitted.den))
                    break
            if (i,j) not in learned and witness_allowed:
                expensive_attempts[i]+=1
                for den in candidates[:3]:
                    fitted=prune_denominator(den,times,values)
                    if fitted is not None:
                        learned[i,j]=fitted;learned_now+=1
                        rowqs[i]=lcm(rowqs[i],poly(fitted.den));qs=lcm(qs,poly(fitted.den))
                        break
            if (i,j) not in learned and witness_allowed:
                for den in candidates[:2]:
                    fitted=boost_denominator(den,times,values)
                    if fitted is not None:
                        learned[i,j]=fitted;learned_now+=1
                        rowqs[i]=lcm(rowqs[i],poly(fitted.den));qs=lcm(qs,poly(fitted.den))
                        break
        if learned_now:print(json.dumps(dict(stage='known-pole-learning',new_entries=learned_now,
            remaining=sum(k not in learned for k in active))),flush=True)
        if all(k in learned for k in active):
            for k in active:built[k]=learned[k]
            known.extend(block);reports.append(dict(top=[masters[i] for i in block],route='known-pole-numerator-learning',
                pilot_points=len(times),new_entries=learned_now,verification_margin_minimum=3))
            save('running');continue
        train=times[:-args.internal_holdouts];holdout=times[-args.internal_holdouts:]
        attempts=[];solved=False
        settings=[(a,b,c,'endpoints','row-diagonal') for a,b,c in
                  [(0,0,1),(0,2,2),(0,3,3),(1,3,3),(2,3,3)]]
        settings += [(a,3,3,'endpoints',scope) for scope in ['row-column','block-column'] for a in [0,1,2]]
        if args.relax_all_poles:settings += [(0,3,3,'all-candidate-poles','block-column')]
        # Reuse a successful envelope class as a first trial, rather than
        # repeating all previously rejected small classes at every block.
        previous=[r for r in reports if r.get('route')=='analyticity-bootstrap' and len(r['top'])==len(block)]
        if previous:
            r=previous[-1];preferred=(r['power'],r['growth'],r['gauge'],'endpoints',r.get('pole_scope','block-column'))
            settings=[preferred]+[s for s in settings if s!=preferred]
        for power,growth,gauge,policy,pole_scope in settings:
            denominators=[];caps=[]
            for i in block:
                ds=[];cs=[]
                for j in sub:
                    if (i,j) not in active:ds.append([1]);cs.append(-1)
                    elif (i,j) in learned:
                        r=learned[i,j];ds.append(r.den);cs.append(len(r.num)-1)
                    else:
                        qa=bare_rows[i] if pole_scope=='row-diagonal' else rowqs[i] if pole_scope=='row-column' else qs
                        qb=diagonal_cols[j] if pole_scope=='row-diagonal' else cols[j]
                        den=lcm(lcm(qa,qb),ep) if power==0 else (qa*qb*ep)**power
                        ds.append([int(c) for c in monic(den)]);cs.append(int(den.degree())+growth)
                denominators.append(ds);caps.append(cs)
            h=ep
            if policy=='all-candidate-poles':
                h=one
                for row in denominators:
                    for d in row:h=lcm(h,poly(d))
            h=h**gauge;shift=RatFunc.make([int(v) for v in h.derivative()],[int(v) for v in h],p)
            be=dict(baseB);ae=dict(baseA)
            for i in range(len(sub)):be[i,i]=be.get((i,i),zero)+shift
            for i in range(len(block)):ae[i,i]=ae.get((i,i),zero)+shift
            spec=BootstrapSpec.from_inputs(family={'source':'full-J-native-points'},variable='tau',assignments={},
                masters=[masters[i] for i in sub+block],target=None,
                de={'A':RatMat(ae,(len(block),len(block)),p),'B':RatMat(be,(len(sub),len(sub)),p)},
                denominators=denominators,numerator_caps=caps,prime=p,mode='de')
            started=perf_counter();print(json.dumps(dict(stage='construct',power=power,growth=growth,gauge=gauge,policy=policy,pole_scope=pole_scope,
                ambient=sum(c+1 for row in caps for c in row))),flush=True)
            if args.backend=='jets':ansatz=build_de_ansatz(spec)
            elif args.backend=='PL-structural':ansatz=build_pl_ansatz(spec,backend='PL-Hasse',laurent_mode='structural-refined-all',extra_order=getattr(args,'extra_order',0))
            else:ansatz=build_pl_ansatz(spec,backend='PL-Hasse',extra_order=getattr(args,'extra_order',0))
            fixed={(a,b):learned[i,j] for a,i in enumerate(block) for b,j in enumerate(sub) if (i,j) in learned and (i,j) in active}
            def sample(t):return dict(prime=p,x=(t-center)%p,
                values=[[all_points[t][masters[i]].get(masters[j],0) for j in sub] for i in block])
            try:
                fit=joint_fit(ansatz,[sample(t) for t in train],fixed,[(t-center)%p for t in times])
                if not fit['compatible']:raise InconsistentSampleError('joint equations inconsistent')
                if fit['parameters'] is None:
                    attempts.append(dict(power=power,growth=growth,gauge=gauge,policy=policy,dimension=ansatz.dimension,
                        rank=fit['rank'],status='need_points',construction_seconds=perf_counter()-started))
                    if policy=='all-candidate-poles':break
                    continue
                for t in holdout:
                    if ansatz.evaluate((t-center)%p,fit['parameters'])!=sample(t)['values']:
                        raise ValueError('holdout disagrees with fitted ansatz')
                answer=assemble_matrix(ansatz,fit['parameters'])
                for a,i in enumerate(block):
                    for b,j in enumerate(sub):
                        r=answer.entries.get((a,b),zero)
                        if not r.is_zero():built[i,j]=r
                report=dict(top=[masters[i] for i in block],route='analyticity-bootstrap',
                    power=power,growth=growth,gauge=gauge,extra_order=getattr(args,'extra_order',0),ambient=sum(c+1 for row in caps for c in row),
                    normalization_policy=policy,
                    pole_scope=pole_scope,
                    fixed_entries=len(fixed),fixed_equations=fit['fixed_equation_count'],
                    dimension=ansatz.dimension,training_points=len(train),holdout_points=len(holdout),
                    joint_fit={k:v for k,v in fit.items() if k!='parameters'},
                    seconds=perf_counter()-started,attempts=attempts,status=ansatz.status)
                reports.append(report);known.extend(block);save('running');solved=True
                print(json.dumps(report),flush=True);break
            except ValueError as e:
                if not isinstance(e,InconsistentSampleError) and str(e)!='holdout disagrees with fitted ansatz':raise
                attempts.append(dict(power=power,growth=growth,gauge=gauge,policy=policy,status='inconsistent',
                    dimension=ansatz.dimension,seconds=perf_counter()-started,error=str(e)))
        if not solved:
            save('need_points_or_envelope',pending_block=[masters[i] for i in block],attempts=attempts,
                 available_points=len(times))
            print(json.dumps(dict(status='need_points_or_envelope',attempts=attempts)),flush=True)
            unresolved.append(dict(block=[masters[i] for i in block],reason='need_points_or_envelope',attempts=attempts));continue
    if unresolved:
        first=next((v for v in unresolved if v['reason']!='waiting_for_dependencies'),unresolved[0])
        save('need_points_or_envelope',pending_block=first['block'],pending_blocks=unresolved,attempts=first.get('attempts',[]))
        return 2
    save('complete',all_reconstructed=True)
    print(json.dumps(dict(status='complete',masters=n,blocks=len(groups))),flush=True)
    return 0


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('work');p.add_argument('--family',required=True)
    p.add_argument('--basis',required=True);p.add_argument('--points',action='append',required=True)
    p.add_argument('--diagonals',action='append',required=True);p.add_argument('--partial',action='append',default=[])
    p.add_argument('--backend',choices=['PL-Hasse','PL-structural','jets'],default='PL-Hasse')
    p.add_argument('--extra-order',type=int,default=0)
    p.add_argument('--witness-budget',type=int,default=4)
    p.add_argument('--internal-holdouts',type=int,default=3)
    p.add_argument('--relax-all-poles',action='store_true')
    p.add_argument('--restart',action='store_true');args=p.parse_args()
    Path(args.work).mkdir(parents=True,exist_ok=True)
    with (Path(args.work)/'run.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        raise SystemExit(main(args))
