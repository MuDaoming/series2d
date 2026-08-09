#!/usr/bin/env python3
"""Fixed non-interactive helper for mdm@tianhe.

This helper deliberately uses the user's working SSH config alias `tianhe` and
BatchMode. It must fail fast instead of falling back to password prompts.
"""

import argparse
import subprocess
import sys
from pathlib import Path


HOST = "tianhe"
TIMEOUT = 60

SSH = [
    "ssh",
    "-T",
    "-o", "BatchMode=yes",
    "-o", "ConnectTimeout=10",
    HOST,
]

SCP_BASE = [
    "scp",
    "-B",
    "-o", "BatchMode=yes",
    "-o", "ConnectTimeout=10",
]


def run_checked(argv, timeout=TIMEOUT):
    proc = subprocess.run(argv, text=True, timeout=timeout)
    return proc.returncode


def remote(command):
    return run_checked(SSH + [command])


def put(local, remote_path):
    return run_checked(SCP_BASE + [local, f"{HOST}:{remote_path}"])


def get(remote_path, local):
    return run_checked(SCP_BASE + [f"{HOST}:{remote_path}", local])


def job(job_id):
    return remote(f"squeue -j {job_id}")


def submit(slurm_path):
    slurm = Path(slurm_path)
    return remote(f"cd {slurm.parent} && sbatch {slurm.name}")


def main():
    parser = argparse.ArgumentParser(description="BatchMode SSH helper for tianhe")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run")
    p_run.add_argument("command")

    p_put = sub.add_parser("put")
    p_put.add_argument("local")
    p_put.add_argument("remote")

    p_get = sub.add_parser("get")
    p_get.add_argument("remote")
    p_get.add_argument("local")

    p_job = sub.add_parser("job")
    p_job.add_argument("job_id")

    p_submit = sub.add_parser("submit")
    p_submit.add_argument("slurm_path")

    args = parser.parse_args()
    if args.cmd == "run":
        rc = remote(args.command)
    elif args.cmd == "put":
        rc = put(args.local, args.remote)
    elif args.cmd == "get":
        rc = get(args.remote, args.local)
    elif args.cmd == "job":
        rc = job(args.job_id)
    elif args.cmd == "submit":
        rc = submit(args.slurm_path)
    else:
        raise AssertionError(args.cmd)
    sys.exit(rc)


if __name__ == "__main__":
    main()
