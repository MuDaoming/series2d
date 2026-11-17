"""Verify preserved full e6 reduction matrices without rerunning AB or IBP."""
import gzip,json
from pathlib import Path
results=Path(__file__).resolve().parents[1]/'results'
for case in ['original63','original_small','separated63']:
    summary=json.loads((results/(case+'.json')).read_text())
    module=json.loads((results/(case+'_module.json')).read_text())
    with gzip.open(results/(case+'_reductions.json.gz'),'rt') as stream:
        reference=json.load(stream)
    assert len(reference)==1111
    assert not module['unrepresented_physical_rows']
    for candidate in ['image_ab','image_trivial']:
        with gzip.open(results/(case+'_'+candidate+'.json.gz'),'rt') as stream:
            data=json.load(stream)
        assert data['reductions']==reference,(case,candidate)
        assert data['basis_independent'] and not data['unresolved']
        assert not data['unrepresented_relations']
        for key in ['global_nonzero_rows','global_rank','global_columns','global_nnz']:
            assert data[key]==summary['candidates'][candidate][key]
    print(case+': 1111 exact reduction rows verified for both candidates')
checks=json.loads((results/'local_condition_checks.json').read_text())
assert checks['physical_depth_margin_16']['AB_basis_identical_to_margin_4']
print('Saved checks passed; no AB/FLINT runtime required for this replay')
