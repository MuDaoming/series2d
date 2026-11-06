#!/usr/bin/env python3
"""Dynamic point dispatcher. Each C++ worker owns exactly one t point."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / 'j_integral_global'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze_input(work, source, role, expected_sha):
    target=work/'input'/f'{role}.json'
    target.parent.mkdir(exist_ok=True)
    if not target.exists():
        shutil.copy2(source,target)
    assert sha(target)==expected_sha, f'frozen input differs: {target}'
    return target


def dispatch(tasks, workers, function):
    answers = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(function, task): task for task in tasks}
        for future in as_completed(futures):
            answers.append(future.result())
    return answers


def main():
    p = argparse.ArgumentParser()
    p.add_argument('family', type=Path)
    p.add_argument('basis', type=Path)
    p.add_argument('targets', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--max-points', type=int, default=129)
    args = p.parse_args()
    assert args.workers > 0 and args.max_points >= 17
    family, basis, targets = (x.resolve() for x in [args.family, args.basis, args.targets])
    source_binary = HERE/'build/global_point'
    work = args.output.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = work/'worker/global_point'
    binary.parent.mkdir(exist_ok=True)
    if not binary.exists():
        shutil.copy2(source_binary, binary)
    fingerprint = {str(x): sha(x) for x in [family, basis, targets]}
    fingerprint[str(source_binary)] = sha(binary)
    identity = work/'input_identity.json'
    if identity.exists():
        assert json.loads(identity.read_text()) == fingerprint, 'immutable run inputs changed; use a new output directory'
    else:
        identity.write_text(json.dumps(fingerprint, indent=2)+'\n')
    family, basis, targets = [freeze_input(work, source, role, fingerprint[str(source)]) for source, role in zip([family,basis,targets], ['family','basis','targets'])]

    def worker(point):
        output = work/f'point{point}.json'
        receipt = work/f'point{point}.receipt.json'
        with (work/f'point{point}.lock').open('w') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if receipt.exists():
                data = json.loads(receipt.read_text())
                assert data['returncode'] == 0 and data['closed'] and output.exists(), data
                assert data['output_sha256'] == sha(output)
                return output
            start = time.time()
            temporary = output.with_suffix('.tmp.json')
            command = [str(binary), str(family), str(temporary), '--point', str(point),
                       '--targets-file', str(targets), '--basis-file', str(basis),
                       '--values', '1', '--compact', '1']
            with (work/f'point{point}.log').open('w') as log:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
                returncode = process.wait()
            data = {'point': point, 'worker_pid': process.pid, 'dispatcher_pid': os.getpid(),
                    'hostname': os.uname().nodename, 'job_id': os.environ.get('SLURM_JOB_ID'), 'start': start,
                    'end': time.time(), 'returncode': returncode, 'closed': False}
            if returncode == 0:
                answer = json.loads(temporary.read_text())
                data['closed'] = not answer['unclosed_derivatives']
                data['master_count'] = answer['master_count_at_point']
                data['target_count'] = answer['target_count']
                assert set(answer['DE_values']) == set(answer['masters_at_point'])
                if data['closed']:
                    temporary.replace(output)
                    data['output_sha256'] = sha(output)
            receipt.write_text(json.dumps(data, indent=2)+'\n')
            if returncode or not data['closed']:
                raise RuntimeError(f'point {point} failed or did not close; preserved {receipt}')
            print(f'point {point} closed in {data["end"]-start:.2f}s', flush=True)
            return output

    counts = [n for n in [17, 33, 65, 129, 257, 513] if n <= args.max_points]
    if args.max_points not in counts:
        counts.append(args.max_points)
    previous = 0
    for count in counts:
        dispatch(list(range(17+previous, 17+count)), args.workers, worker)
        paths = [work/f'point{point}.json' for point in range(17, 17+count)]
        manifest = work/'reconstruction_manifest.json'
        manifest.write_text(json.dumps({'family': str(family), 'points': [str(x) for x in paths]}, indent=2)+'\n')
        result = subprocess.run([str(HERE/'build/reconstruct'), str(manifest), str(work/'reduction.json')])
        if result.returncode == 0:
            return
        if result.returncode != 2:
            raise RuntimeError('reconstruction program failed')
        previous = count
    print('More sample points required; validated points retained.', flush=True)
    return 2


if __name__ == '__main__':
    raise SystemExit(main())
