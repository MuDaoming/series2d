"""Complete the scalar-product basis and run Blade maximal-cut MI discovery."""
import argparse,json,itertools,subprocess,shutil,hashlib
from fractions import Fraction as Q
from pathlib import Path

def rank(rows):
 a=[[Q(v) for v in r] for r in rows];k=0
 for j in range(len(a[0])):
  pivot=next((i for i in range(k,len(a)) if a[i][j]),None)
  if pivot is None:continue
  a[k],a[pivot]=a[pivot],a[k];v=a[k][j];a[k]=[x/v for x in a[k]]
  for i in range(k+1,len(a)):
   v=a[i][j]
   if v:a[i]=[x-v*y for x,y in zip(a[i],a[k])]
  k+=1
  if k==len(a):break
 return k

def linear(coefficients,names):
 terms=[]
 for c,n in zip(coefficients,names):
  c=Q(c)
  if not c:continue
  terms.append('('+str(c)+')*'+n)
 return '+'.join(terms) or '0'

def main(a):
 cfg=json.loads(Path(a.family).read_text());e=len(cfg['gram']);names=['l1','l2']+['p'+str(i+1) for i in range(e)]
 def row(loop,shift):
  x,y=map(Q,loop);return [x*x,y*y,2*x*y]+[2*x*Q(q) for q in shift]+[2*y*Q(q) for q in shift]
 props=[];rows=[]
 for p in cfg['propagators']:
  props.append('('+linear(p['loop']+p['shift'],names)+')^2-('+str(p['mass2'])+')');rows.append(row(p['loop'],p['shift']))
 if rank(rows)!=len(rows):raise ValueError('physical inverse propagators are linearly dependent in scalar products')
 for b,i in itertools.product(range(2),range(e)):
  loop=[int(b==0),int(b==1)];shift=[int(j==i) for j in range(e)];candidate=row(loop,shift)
  if rank(rows+[candidate])>len(rows):props.append('('+linear(loop+shift,names)+')^2');rows.append(candidate)
 if len(rows)!=3+2*e:raise ValueError('could not complete irreducible scalar-product basis')
 out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=True)
 if (out/'masters.json').exists():print((out/'masters.json').read_text());return
 packet=dict(name=a.name,dimension=cfg['dimension'],propagators=props,loops=['l1','l2'],external_momenta=names[2:],scalar_products=[['p'+str(i+1)+'*p'+str(j+1),str(cfg['gram'][i][j])] for i in range(e) for j in range(i,e)],top=[1]*len(cfg['propagators'])+[0]*(len(rows)-len(cfg['propagators'])),physical_count=len(cfg['propagators']),blade_root=str(Path(a.blade).resolve()),source_sha256=hashlib.sha256(Path(a.family).read_bytes()).hexdigest())
 (out/'input.json').write_text(json.dumps(packet,indent=2)+'\n');script=out/'run.wl';shutil.copy2(Path(__file__).with_suffix('.wl'),script)
 with (out/'run.log').open('w') as log:proc=subprocess.run([a.kernel,'-noprompt','-script',str(script)],cwd=out,stdout=log,stderr=subprocess.STDOUT)
 if proc.returncode or not (out/'masters.json').exists():raise RuntimeError('MI-only discovery failed: '+str(out/'run.log'))
 print((out/'masters.json').read_text())

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',required=True);p.add_argument('--name',required=True);p.add_argument('--output',required=True);p.add_argument('--blade',default='/home/mudaoming/Documents/Packages/blade/blade-main');p.add_argument('--kernel',default='/usr/local/bin/WolframKernel');main(p.parse_args())
