#!/usr/bin/env python3
"""Compare the polynomial backend against saved reference J jets, then time it."""
from pathlib import Path
import argparse
import json
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from series_engine import Model
from native_engine import generate_native


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('reference',type=Path,help='directory with fbi_model.json and series.json')
    ap.add_argument('output',type=Path)
    ap.add_argument('--orders',type=int,nargs='+',default=[20,100,500,2000])
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    before=time.perf_counter();model=Model(json.loads((args.reference/'fbi_model.json').read_text()))
    prep=time.perf_counter()-before
    saved=json.loads((args.reference/'series.json').read_text());keys=[tuple(x) for x in saved['period_keys']]
    reference={label:[{key:columns[i][n] for i,key in enumerate(keys) if columns[i][n]}
                       for n in range(saved['order']+1)] for label,columns in saved['targets'].items()}
    records=[]
    for order in args.orders:
        directory=args.output/f'order_{order}';directory.mkdir(exist_ok=True)
        result=generate_native(model,order,verify_all=order<=24,output_dir=directory)
        overlap=min(order,saved['order'])+1
        assert set(result['coordinates'])==set(reference)
        for label in reference:
            assert result['coordinates'][label][:overlap]==reference[label][:overlap],label
        record={'order':order,'reference_coefficients_matched_per_target':overlap,
                'timings':result['timings'],'operations':result['operations']}
        records.append(record);print(json.dumps(record),flush=True)
        (args.output/'benchmark.json').write_text(json.dumps({'model_preparation_seconds':prep,'runs':records},indent=2))


if __name__=='__main__':main()
