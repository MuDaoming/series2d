"""Find exact (m,0) singleton bubbles with a matching threshold propagator."""
import json,itertools
from fractions import Fraction
from pathlib import Path
base=Path(__file__).resolve().parent;root=base.parents[3]
f=json.load(open(root/'j_integral/tests/dp/input/family_request.json'));br=f['j']['branches'];P=f['family']['propagators'];rows=json.load(open(root/'j_integral/issues/sector_symmetries/audit_with_cross_sector.json'))['dp'];found=[]
for row in rows:
 groups=[[i for i,v in enumerate(row['mask']) if v=='1' and br[i]==b] for b in range(3)];matches=[]
 for b1,b2 in itertools.combinations(range(3),2):
  if len(groups[b1])!=1 or len(groups[b2])!=1:continue
  i,j=groups[b1][0],groups[b2][0];mi,mj=Fraction(P[i]['mass2']),Fraction(P[j]['mass2'])
  if not ((mi==0 and mj>0) or (mj==0 and mi>0)):continue
  remain=next(b for b in range(3) if b not in [b1,b2]);m=max(mi,mj)
  for k in groups[remain]:
   if Fraction(P[k]['mass2'])!=m:continue
   tri={b1:i,b2:j,remain:k};vec=[Fraction(P[tri[2]]['shift'][a])-Fraction(P[tri[0]]['shift'][a])-Fraction(P[tri[1]]['shift'][a]) for a in range(4)]
   if not any(vec):matches.append({'bubble':[i+1,j+1],'threshold_line':k+1,'m2':str(m)})
 if matches:
  assert row['FI']==0
  found.append({'mask':row['mask'],'FI':row['FI'],'J':row['J'],'matches':matches})
(base/'dp_matches.json').write_text(json.dumps(found,indent=2)+'\n');print('Matched raw sectors:',len(found))
