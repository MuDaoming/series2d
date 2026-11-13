"""Collect independent Blade maximal-cut and J point outputs without rerunning jobs."""
import json
from pathlib import Path
base=Path(__file__).resolve().parent;rows=[]
for path in sorted(base.rglob('result.json')):
 if 'cache' in path.parts:continue
 if not path.with_name('request.json').exists() or not path.with_name('j_result.json').exists():continue
 fi=json.load(open(path));j=json.load(open(path.with_name('j_result.json')))
 ordinary=[n for n in j['basis'] if n.startswith('J')];top=[n for n in ordinary if all(int(v)>0 for v in n.split('{')[1][:-1].split(','))]
 assert j['status']=='complete' and not j['unclosed_derivatives']
 row={'case':str(path.parent.relative_to(base)),'FI_top':fi['count'],'J_top':len(top),'J_ordinary_all_sectors':len(ordinary),'boundary_count':sum(n.startswith('B') for n in j['basis']),'Blade_seconds':fi['seconds'],'prime':j['point']['prime'],'dimension':j['point']['dimension'],'t':j['point']['t'],'J_DE_closed':True};rows.append(row);print(row,flush=True)
(base/'counts.json').write_text(json.dumps(rows,indent=2)+'\n')
