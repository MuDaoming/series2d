#!/usr/bin/env python3
import argparse
import concurrent.futures
import subprocess
import time
from pathlib import Path


def sector_name(vec_text):
    return "".join(x.strip() for x in vec_text.strip().strip("{}").split(","))


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def run_one(root, searcher, cut_dir, force):
    out = cut_dir / "poly_relation"
    log_path = cut_dir / "poly_relation.log"
    if out.exists() and out.stat().st_size > 0 and not force:
        return cut_dir.name, 0, 0.0, "skip"
    t0 = time.time()
    cmd = [
        str(searcher),
        str(cut_dir / "config"),
        str(cut_dir / "G"),
        str(cut_dir / "series_list"),
        str(out),
    ]
    with log_path.open("w") as log:
        rc = subprocess.run(cmd, cwd=root, stdout=log, stderr=subprocess.STDOUT).returncode
    elapsed = time.time() - t0
    (cut_dir / "search_elapsed_sec").write_text(f"{elapsed:.6f}\n")
    (cut_dir / "search_return_code").write_text(f"{rc}\n")
    return cut_dir.name, rc, elapsed, "run"


def main():
    ap = argparse.ArgumentParser(description="Run poly_relation_searcher for every cut.")
    ap.add_argument("--root", default=".")
    ap.add_argument("--family", required=True)
    ap.add_argument("--series-dir-name", required=True)
    ap.add_argument("--out-name", required=True)
    ap.add_argument("--deg", type=int, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--ncheck", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    run_root = root / "pipeline" / "runs" / args.family
    series_root = run_root / args.series_dir_name
    out_root = run_root / args.out_name
    searcher = root / "search" / "tools" / "poly_relation_searcher" / "poly_relation_searcher"

    p_value = None
    n_value = None
    for line in (run_root / "input" / "base_config").read_text().splitlines():
        if line.startswith("p ="):
            p_value = line.split("=", 1)[1].strip()
        elif line.startswith("N ="):
            n_value = line.split("=", 1)[1].strip()
    if p_value is None or n_value is None:
        raise RuntimeError("missing p or N in base_config")

    cut_sectors = [
        line.strip()
        for line in (run_root / "metadata" / "cut_sectors").read_text().splitlines()
        if line.strip()
    ]

    cut_dirs = []
    for sec in cut_sectors:
        name = sector_name(sec)
        sdir = series_root / name
        odir = out_root / name
        odir.mkdir(parents=True, exist_ok=True)
        write(odir / "G", (run_root / "search" / name / "G").read_text())
        write(
            odir / "config",
            f"N = {n_value}\n"
            f"deg = {args.deg}\n"
            f"m = {args.m}\n"
            f"ncheck = {args.ncheck}\n"
            f"p = {p_value}\n",
        )
        write(odir / "series_list", f"{sdir / 'series'} {sdir / 'target'}\n")
        cut_dirs.append(odir)

    print(
        f"start family={args.family} cuts={len(cut_dirs)} deg={args.deg} "
        f"m={args.m} jobs={args.jobs} out={out_root}",
        flush=True,
    )
    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futures = [ex.submit(run_one, root, searcher, d, args.force) for d in cut_dirs]
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
