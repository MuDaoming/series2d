"""Check exact U,F-preserving propagator permutations into other sectors."""
from audit import ROOT,polynomials
import json,itertools
from pathlib import Path
import sympy as s
report=json.load(open(Path(__file__).with_name('audit.json')))
for case,rows in report.items():
    f=json.load(open(ROOT/f'j_integral/tests/{case}/input/family_request.json'));br=f['j']['branches'];n=len(br);x,U,F=polynomials(f);terms=s.Poly(F,*x).terms();cache={}
    def data(active):
        key=tuple(active)
        if key not in cache:
            rest=set(range(n))-set(active)
            ts=[(tuple(ex[i] for i in active),c) for ex,c in terms if all(ex[j]==0 for j in rest)]
            cache[key]=(dict(ts),[[i for i in active if br[i]==b] for b in range(3)])
        return cache[key]
    for row in rows:
        source=[i for i,q in enumerate(row['mask']) if q=='1'];original,groups=data(source);matches=[]
        for target in itertools.combinations(range(n),len(source)):
            targetterms,tg=data(target);pos={i:k for k,i in enumerate(target)};found=None
            for bp in itertools.permutations(range(3)):
                if any(len(groups[b])!=len(tg[bp[b]]) for b in range(3)):continue
                for orders in itertools.product(*(itertools.permutations(tg[bp[b]]) for b in range(3))):
                    mapping={i:j for b in range(3) for i,j in zip(groups[b],orders[b])};per=[pos[mapping[i]] for i in source];trans={}
                    for ex,c in original.items():
                        out=[0]*len(source)
                        for k,v in enumerate(ex):out[per[k]]=v
                        trans[tuple(out)]=c
                    if trans==targetterms:
                        found={'mask':''.join('1' if i in target else '0' for i in range(n)),'branch_map':list(bp),'propagator_map':{str(i+1):j+1 for i,j in mapping.items()}};break
                if found:break
            if found:matches.append(found)
        row['equivalent_sectors']=matches
        print(case,row['mask'],'equivalent',[q['mask'] for q in matches],flush=True)
Path(__file__).with_name('audit_with_cross_sector.json').write_text(json.dumps(report,indent=2)+'\n')
