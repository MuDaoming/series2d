"""Run unmodified upstream builders in an isolated tree; atomically publish.

No upstream source file or executing shared kernel is overwritten in place.
"""
from pathlib import Path
import os,sys,shutil,tempfile,subprocess,json,hashlib,re
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT/'extern/analytic_bootstrap-master';OUT=ROOT/'build/production'
OUT.mkdir(parents=True,exist_ok=True)
PREFIX=Path(os.environ.get('FLINT_PREFIX','/usr/local')).resolve()
FLINT_VERSION=re.search(r'#define\s+FLINT_VERSION\s+"([^"]+)"',(PREFIX/'include/flint/flint.h').read_text()).group(1)
with tempfile.TemporaryDirectory(prefix='ab-runtime-build-',dir=ROOT/'build') as tmp:
    stage=Path(tmp);(stage/'tools').mkdir()
    for name in ['analytic_bootstrap','native']:(stage/name).symlink_to(SOURCE/name,target_is_directory=True)
    for name in ['build_native_pole_kernel.py','build_laurent_kernel.py','build_coordinate_kernel.py','build_pure_runtime.py']:
        shutil.copyfile(SOURCE/'tools'/name,stage/'tools'/name)
    pc=stage/'pkgconfig';pc.mkdir()
    # Upstream requires pkg-config metadata; this local file describes the
    # same installed FLINT headers/libraries used by the existing J workers.
    (pc/'flint.pc').write_text('prefix='+str(PREFIX)+'\nlibdir=${prefix}/lib\nincludedir=${prefix}/include\nName: FLINT\nDescription: Installed FLINT\nVersion: '+FLINT_VERSION+'\nLibs: -L${libdir} -lflint -lgmp -lmpfr\nCflags: -I${includedir}\n')
    env=dict(os.environ,FLINT_PREFIX=str(PREFIX),PKG_CONFIG_PATH=str(pc))
    commands=[]
    for flag in [None,'--hasse','--derivative','--support','--support-legacy','--conditions']:
        commands.append([sys.executable,str(stage/'tools/build_native_pole_kernel.py'),*([flag] if flag else [])])
    for flag in [None,'--higher','--fused','--higher-fused']:
        commands.append([sys.executable,str(stage/'tools/build_laurent_kernel.py'),*([flag] if flag else [])])
    commands.extend([[sys.executable,str(stage/'tools'/name)] for name in ['build_coordinate_kernel.py','build_pure_runtime.py']])
    for command in commands:subprocess.run(command,env=env,check=True)
    output=stage/'build/hjj_native';records={}
    for path in sorted(output.iterdir()):
        if path.is_file():
            records[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
            os.replace(path,OUT/path.name)
    manifest=OUT/'runtime.json.tmp';manifest.write_text(json.dumps(dict(files=records,flint_version=FLINT_VERSION,builders={name:hashlib.sha256((stage/'tools'/name).read_bytes()).hexdigest() for name in ['build_native_pole_kernel.py','build_laurent_kernel.py','build_coordinate_kernel.py','build_pure_runtime.py']}),indent=2));os.replace(manifest,OUT/'runtime.json')
print(OUT,flush=True)
