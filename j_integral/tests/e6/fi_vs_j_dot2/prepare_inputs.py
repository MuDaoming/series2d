import json,itertools
from pathlib import Path
base=Path(__file__).resolve().parent;root=base.parents[3]
f=json.load(open(root/'j_integral/tests/e6/input/family_request.json'));blade=json.load(open(root/'j_integral/archive/j_integral_bootstrap/reference/fi_master_counts/e6_input.json'));mi=json.load(open(root/'j_integral/archive/j_integral_bootstrap/reference/fi_master_counts/e6.json'));alltargets=[];common=[]
for mask in itertools.product([0,1],repeat=7):
 if not any(mask):continue
 active=[i for i,x in enumerate(mask) if x]
 for dot in range(3):
  for ds in itertools.combinations_with_replacement(active,dot):
   nu=list(mask)
   for i in ds:nu[i]+=1
   alltargets.append(nu)
   if all(any(mask[i] for i,b in enumerate(f['j']['branches']) if b==k) for k in range(3)):common.append(nu)
assert len(common)==1111
for name,obj in [('fi_family',blade),('fi_masters',mi['master_indices']),('all_targets',alltargets),('common_targets',common)]:
 (base/'input'/f'{name}.json').write_text(json.dumps(obj,indent=2)+'\n')
print('all nonempty sectors dot<=2',len(alltargets),'three-branch common',len(common))
