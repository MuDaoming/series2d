"""Build and sample one requested physical sector using generic compression."""
import sys,json,subprocess,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from upstream_j import ROOT
from restriction_inputs import make_restriction,write_restriction,expand_name
from restricted_numeric import indices
O=ROOT/'results/dp_reduction_assist';BINARY=ROOT/'build/global_point_subset'

def main():
    p=argparse.ArgumentParser();p.add_argument('--sector',required=True);p.add_argument('--work',required=True);a=p.parse_args();O=Path(a.work);O.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((ROOT.parent/'j_integral_global/examples/dp_physical.json').read_text());full=json.loads((ROOT/'results/dp_full/points/point17.json').read_text());D=json.loads((ROOT/'results/dp_final/connection_frozen.json').read_text());S={i for i,c in enumerate(a.sector) if c=='1'}
    targets=[n for n in full['reduction_values'] if n.startswith('J') and {i for i,v in enumerate(indices(n)[1]) if v}==S]
    closure=set().union(*(set(full['reduction_values'][n]) for n in targets));front=list(closure)
    while front:
        for j in D['de'][front.pop()]['rhs']:
            if j not in closure:closure.add(j);front.append(j)
    names=[n for n in D['masters'] if n in closure];packet=make_restriction(cfg,names,[indices(n)[1] for n in targets]);write_restriction(O,packet);(O/'original_targets.json').write_text(json.dumps(targets));(O/'no_derivatives.json').write_text('[]')
    (O/'points').mkdir(exist_ok=True)
    cmd=[str(BINARY),str(O/'family.json'),str(O/'natural17.json'),'--derivative-names-file',str(O/'no_derivatives.json'),'--point','17','--targets-file',str(O/'targets.json'),'--integration','atomic','--values','1','--compact','1']
    subprocess.run(cmd,check=True)
    d=json.loads((O/'natural17.json').read_text());natural=d['masters_at_point'];expanded=[expand_name(n,packet['mapping']) for n in natural]
    if not set(expanded)<=set(D['masters']):raise ValueError('natural basis requires a true basis transformation')
    (O/'closure_derivatives.json').write_text(json.dumps(packet['basis']))
    packet['basis']=natural;packet['mapping'].update(original_masters=expanded,compressed_masters=natural)
    # Keep original closure parents; extra boundary coordinates need no DE.
    # Extra natural boundary objects are R coordinates, not derivative targets;
    # adding their J parents would enlarge the requested R space unnecessarily.
    write_restriction(O,packet)
    checks=0
    for t in [17,18]:
        out=O/'points'/f'point{t}.json';subprocess.run([str(BINARY),str(O/'family.json'),str(out),'--derivative-names-file',str(O/'closure_derivatives.json'),'--point',str(t),'--basis-file',str(O/'basis.json'),'--targets-file',str(O/'targets.json'),'--integration','atomic','--values','1','--compact','1'],check=True)
        d=json.loads(out.read_text());assert not d['unclosed_derivatives'];old=json.loads((ROOT/f'results/dp_full/points/point{t}.json').read_text())
        for local,row in d['reduction_values'].items():
            name=expand_name(local,packet['mapping'])
            if name not in targets:continue
            expanded_row={expand_name(c,packet['mapping']):v for c,v in row.items()}
            assert expanded_row==old['reduction_values'][name],name;checks+=len(D['masters'])
    report=dict(sector=a.sector,targets=targets,closure_rank=len(names),natural_rank=len(natural),FBI_rank=d['FBI_rank'],active_positions=packet['mapping']['active_positions'],full_point_checks=checks)
    (O/'preparation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
if __name__=='__main__':main()
