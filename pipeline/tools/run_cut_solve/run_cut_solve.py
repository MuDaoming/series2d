#!/usr/bin/env python3
import argparse
import concurrent.futures
import subprocess
import time
from pathlib import Path


def run_one(root, solver, cut_dir, delta, force):
    out = cut_dir / "integral_solution"
    log_path = cut_dir / "integral_solver.log"
    if out.exists() and out.stat().st_size > 0 and not force:
        return cut_dir.name, 0, 0.0, "skip"
    t0 = time.time()
    cmd = [
        str(solver),
        str(cut_dir / "G"),
        str(cut_dir / "poly_relation"),
        str(delta),
        str(out),
    ]
    with log_path.open("w") as log:
        rc = subprocess.run(cmd, cwd=root, stdout=log, stderr=subprocess.STDOUT).returncode
    elapsed = time.time() - t0
    (cut_dir / "solve_elapsed_sec").write_text(f"{elapsed:.6f}\n")
    (cut_dir / "solve_return_code").write_text(f"{rc}\n")
    return cut_dir.name, rc, elapsed, "run"


def main():
    ap = argparse.ArgumentParser(description="Run integral_solver for every cut search output.")
    ap.add_argument("--root", default=".")
    ap.add_argument("--family", required=True)
    ap.add_argument("--search-dir-name", required=True)
    ap.add_argument("--delta", default="571")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    run_root = root / "pipeline" / "runs" / args.family
    search_root = run_root / args.search_dir_name
    solver = root / "search" / "tools" / "integral_solver" / "integral_solver"

    cut_dirs = sorted([d for d in search_root.iterdir() if d.is_dir()])
    print(
        f"start family={args.family} cuts={len(cut_dirs)} "
        f"delta={args.delta} jobs={args.jobs} dir={search_root}",
        flush=True,
    )

    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futures = [ex.submit(run_one, root, solver, d, args.delta, args.force) for d in cut_dirs]
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            name, rc, elapsed, mode = fut.result()
            done += 1
            print(f"{mode} {name} rc={rc} elapsed={elapsed:.2f} done={done}/{len(cut_dirs)}", flush=True)
            if rc != 0:
                failed.append(name)

    print("finished failed=" + ",".join(failed), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
