"""Build portable project-owned native executables against a chosen FLINT."""
import argparse,json,hashlib,os,subprocess,fcntl
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent
p=argparse.ArgumentParser();p.add_argument('--prefix',default=os.environ.get('FLINT_PREFIX','/usr/local'));p.add_argument('--json-include');a=p.parse_args();prefix=Path(a.prefix).resolve();out=ROOT/'build';out.mkdir(exist_ok=True)
source_deps=[PROJECT/'j_integral/native/series.cpp',PROJECT/'j_integral/native/prepare.hpp',PROJECT/'j_integral_methods/native/direct_j.cpp',PROJECT/'j_integral_global/native/triangular.hpp']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
with (out/'backend-build.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 for name,source in [('global_point_subset',PROJECT/'j_integral_global/native/global_point.cpp'),('reconstruct_partial',ROOT/'native/reconstruct_partial.cpp'),('verify_connection',ROOT/'native/verify_connection.cpp')]:
  target=out/name;receipt=out/(name+'.build.json');tmp=out/(name+'.tmp-'+str(os.getpid()))
  command=[os.environ.get('CXX','g++'),'-O3','-std=c++17','-I'+str(prefix/'include')]
  if a.json_include:command+=['-I'+str(Path(a.json_include).resolve())]
  command += [str(source),'-o',str(tmp),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lflint','-lgmp']
  compiler_major=int(subprocess.check_output([command[0],'-dumpversion'],text=True).strip().split('.')[0])
  if compiler_major<9:command+=['-lstdc++fs']
  identity=dict(source_hashes={str(v):sha(v) for v in [source,*source_deps]},command=command)
  if receipt.exists() and target.exists():
   saved=json.loads(receipt.read_text())
   old=dict(saved);old.pop('binary_sha256',None)
   # Temporary PID does not affect the build identity.
   old['command'][old['command'].index('-o')+1]=str(tmp)
   if old==identity and saved['binary_sha256']==sha(target):continue
  subprocess.run(command,check=True);os.replace(tmp,target);identity['binary_sha256']=sha(target)
  meta=receipt.with_suffix('.tmp.json');meta.write_text(json.dumps(identity,indent=2)+'\n');meta.replace(receipt)
  print(target,flush=True)
