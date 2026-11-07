"""Exercise production publication logic without recompiling shared libraries."""
from pathlib import Path
import json,tempfile
from unittest.mock import patch
REAL=Path(__file__).resolve().parents[1]
source=(REAL/'tests/build_kernels.py').read_text()
with tempfile.TemporaryDirectory(prefix='ab-atomic-regression-') as tmp:
    root=Path(tmp);(root/'tests').mkdir();(root/'extern').mkdir()
    (root/'extern/analytic_bootstrap-master').symlink_to(REAL/'extern/analytic_bootstrap-master',target_is_directory=True)
    out=root/'build/kernels';out.mkdir(parents=True)
    calls=[]
    def compile_stub(command,**kwargs):
        calls.append(command)
        if command[0]=='cc':
            temporary=Path(command[command.index('-o')+1]);assert '.tmp-' in temporary.name
            target=temporary.with_name(temporary.name.split('.tmp-')[0]);target.write_bytes(b'previous-complete-binary')
            temporary.write_bytes(b'new-complete-binary')
            assert target.read_bytes()==b'previous-complete-binary'
        else:
            staged=Path(command[command.index('--output')+1]);assert staged.name.startswith('upstream-build-')
            for name in ['upstream_evaluate.so','upstream_fit.so']: (staged/name).write_bytes(b'new-complete-library')
            (staged/'build.json').write_text('[]\n')
    with patch('subprocess.run',side_effect=compile_stub):
        exec(compile(source,str(root/'tests/build_kernels.py'),'exec'),{'__file__':str(root/'tests/build_kernels.py')})
    assert len(calls)==7
    assert not list(out.glob('*.tmp-*'))
    assert all((out/name).read_bytes()==b'new-complete-binary' for name in ['hasse_kernel','derivative_kernel','support_kernel','conditions_kernel','laurent_kernel','laurent_higher_kernel'])
    assert (root/'build/upstream/upstream_fit.so').read_bytes()==b'new-complete-library'
print('atomic-publication regression passed: 6 executables + 2 shared libraries')
