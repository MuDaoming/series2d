"""Build the collaborator e6 source/input bundle without system libraries or caches."""
import argparse,hashlib,json,shutil,tarfile,tempfile,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
p=argparse.ArgumentParser();p.add_argument('--output',default=str(HERE.parent/'e6_ab_20261112.tar.gz'));p.add_argument('--directory',help='also retain the unpacked delivery tree');a=p.parse_args()
output=Path(a.output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
package_name='e6_ab_20261112'
ignore=shutil.ignore_patterns('__pycache__','*.pyc','build','.venv','.deps','runs')
with tempfile.TemporaryDirectory(prefix='e6-package-') as tmp:
 root=Path(tmp)/package_name;root.mkdir();src=PROJECT/'j_integral'
 for filename in ['README.md','setup.sh','build.sh','run_e6.sh','requirements.txt']:shutil.copy2(HERE/filename,root/filename)
 shutil.copytree(src/'src/IBP',root/'j_integral/src/IBP',ignore=ignore)
 tests=root/'j_integral/tests';tests.mkdir(parents=True)
 for filename in ['backend_worker.py','prepare_run.py','verify_result.py','export_physical.py']:shutil.copy2(src/'tests'/filename,tests/filename)
 shutil.copy2(HERE/'runtime.sh',tests/'runtime.sh')
 shutil.copytree(src/'tests/e6/input',tests/'e6/input')
 ref=tests/'e6/reference';ref.mkdir(parents=True)
 for filename in ['physical_de.json.gz','verification.json','rational_reconstruction_degrees.json','learn_timings.json','solve_timings.json']:shutil.copy2(src/'tests/e6/reference'/filename,ref/filename)
 hashes={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in ref.iterdir()};(ref/'SHA256.json').write_text(json.dumps(hashes,indent=2)+'\n')
 text=(src/'tests/e6/run.sh').read_text().replace('python ', '"$E6_PYTHON" ');(tests/'e6/run.sh').write_text(text)
 case_readme=(src/'tests/e6/README.md').read_text();case_readme=case_readme.replace('运行 `bash j_integral/tests/e6/run.sh RUN_NAME`。','运行包根目录的 `bash run_e6.sh RUN_NAME`。');case_readme+='\n交付包中所有 `archive/...` 字样仅记述历史数据来源，不是运行时路径；所需anchor、对角块、主基和物理参考DE均已置于本例input/reference。\n';(tests/'e6/README.md').write_text(case_readme)
 (tests/'README.md').write_text('# e6 原包 AB 测试\n\n本交付包的输入、依赖、运行和验证方式见根目录README.md；家族具体传播子与运动学见e6/README.md。数学求解完全调用原包learn/solve-de，接口仅作数据转换和数值点缓存。\n')
 ab=src/'archive/j_integral_bootstrap/extern/analytic_bootstrap-master';dest=root/'vendor/analytic_bootstrap';dest.mkdir(parents=True)
 upstream={}
 for directory in ['analytic_bootstrap','native']:shutil.copytree(ab/directory,dest/directory,ignore=ignore)
 (dest/'tools').mkdir()
 for filename in ['build_native_pole_kernel.py','build_laurent_kernel.py','build_coordinate_kernel.py','build_pure_runtime.py']:shutil.copy2(ab/'tools'/filename,dest/'tools'/filename)
 for filename in ['pyproject.toml','README.md','README-SHARE.md']:shutil.copy2(ab/filename,dest/filename)
 (dest/'docs').mkdir()
 for filename in ['USER_GUIDE.md','PL_PRODUCTION_UPSTREAM_IMPLEMENTATION_2026-10-05.md','CONTRACT.md','INSTALL.md']:shutil.copy2(ab/'docs'/filename,dest/'docs'/filename)
 for path in sorted(dest.rglob('*')):
  if path.is_file():
   relative=path.relative_to(dest);assert path.read_bytes()==(ab/relative).read_bytes();upstream[str(relative)]=hashlib.sha256(path.read_bytes()).hexdigest()
 (root/'UPSTREAM_SOURCE_SHA256.json').write_text(json.dumps(upstream,indent=2)+'\n')
 (root/'SOURCE.json').write_text(json.dumps(dict(project_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),case='e6',AB_version='0.3.0',upstream_source_files=len(upstream),system_dependencies_bundled=False,compiled_binaries_bundled=False),indent=2)+'\n')
 for path in [root/'setup.sh',root/'build.sh',root/'run_e6.sh',tests/'e6/run.sh']:path.chmod(0o755)
 (root/'MANIFEST.sha256').write_text(''.join(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(root))+'\n' for path in sorted(root.rglob('*')) if path.is_file()))
 if a.directory:
  directory=Path(a.directory).resolve()
  if directory.exists():raise FileExistsError(directory)
  shutil.copytree(root,directory)
 with tarfile.open(output,'w:gz') as tar:tar.add(root,arcname=package_name)
 digest=hashlib.sha256(output.read_bytes()).hexdigest();output.with_suffix(output.suffix+'.sha256').write_text(digest+'  '+output.name+'\n')
 print(json.dumps(dict(archive=str(output),bytes=output.stat().st_size,sha256=digest,upstream_files=len(upstream)),indent=2))
