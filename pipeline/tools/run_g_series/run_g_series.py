#!/usr/bin/env python3
import argparse
import concurrent.futures
import subprocess
import time
from pathlib import Path


def sector_name(vec_text):
    return "".join(x.strip() for x in vec_text.strip().strip("{}").split(","))


def write_lines(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def copy_config(base_config, out_config, bc_index, bc_count, deg, sectormap_rel):
    lines = []
    for line in base_config.read_text().splitlines():
        if line.startswith("bc ="):
            bc = [0] * bc_count
            bc[bc_index] = 1
            lines.append("bc = {" + ", ".join(map(str, bc)) + "}")
        elif line.startswith("deg ="):
            lines.append(f"deg = {deg}")
        elif line.startswith("reduceMode ="):
            lines.append("reduceMode = cut")
        elif line.startswith("sectorMap ="):
            lines.append(f"sectorMap = {sectormap_rel}")
        else:
            lines.append(line)
    if not any(line.startswith("sectorMap =") for line in lines):
        lines.append(f"sectorMap = {sectormap_rel}")
    write_lines(out_config, lines)


def line_count(path):
    if not path.exists():
        return None
    with path.open() as f:
        return sum(1 for _ in f)


def run_one(root, runner, s_path, cut_dir, force):
    target = cut_dir / "target"
    series = cut_dir / "series"
    log_path = cut_dir / "run.log"

    target_count = line_count(target)
    series_count = line_count(series)
    if (not force) and series.exists() and target_count == series_count:
        return cut_dir.name, 0, 0.0, "skip"

    t0 = time.time()
    cmd = [
        str(runner),
        str(s_path),
        str(cut_dir / "config"),
        str(target),
        str(series),
    ]
    with log_path.open("w") as log:
        rc = subprocess.run(cmd, cwd=root, stdout=log, stderr=subprocess.STDOUT).returncode
    elapsed = time.time() - t0
    (cut_dir / "elapsed_sec").write_text(f"{elapsed:.6f}\n")
    (cut_dir / "return_code").write_text(f"{rc}\n")
    return cut_dir.name, rc, elapsed, "run"


def main():
    ap = argparse.ArgumentParser(
        description="Generate per-cut series using search/<sector>/G as target."
    )
    ap.add_argument("--root", default=".")
    ap.add_argument("--family", required=True)
    ap.add_argument("--deg", type=int, required=True)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--out-name", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    run_root = root / "pipeline" / "runs" / args.family
    out_name = args.out_name or f"search_series_deg{args.deg}"
    out_root = run_root / out_name

    cut_sectors = [
        line.strip()
        for line in (run_root / "metadata" / "cut_sectors").read_text().splitlines()
        if line.strip()
    ]
    base_config = run_root / "input" / "base_config"
    s_path = run_root / "input" / "S"
    runner = root / "expand" / "tools" / "fi_pipeline_runner" / "fi_pipeline_runner"

    summary = ["sector\tG_integral_count\tseries_ready"]
    cut_dirs = []
    for i, sec in enumerate(cut_sectors):
        name = sector_name(sec)
        cut_dir = out_root / name
        cut_dir.mkdir(parents=True, exist_ok=True)
        g_path = run_root / "search" / name / "G"
        target_text = g_path.read_text()
        (cut_dir / "target").write_text(target_text)
        copy_config(
            base_config,
            cut_dir / "config",
            i,
            len(cut_sectors),
            args.deg,
            "../../symmetry/sectormap",
        )
        g_count = sum(1 for line in target_text.splitlines() if line.strip())
        ready = int((cut_dir / "series").exists() and line_count(cut_dir / "series") == g_count)
        summary.append(f"{sec}\t{g_count}\t{ready}")
        cut_dirs.append(cut_dir)

    write_lines(run_root / "metadata" / f"{out_name}_integral_counts", summary)

    print(
        f"start family={args.family} deg={args.deg} cuts={len(cut_dirs)} "
        f"jobs={args.jobs} out={out_root}",
        flush=True,
    )

    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futures = [ex.submit(run_one, root, runner, s_path, d, args.force) for d in cut_dirs]
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            name, rc, elapsed, mode = fut.result()
            done += 1
            print(
                f"{mode} {name} rc={rc} elapsed={elapsed:.2f} done={done}/{len(cut_dirs)}",
                flush=True,
            )
            if rc != 0:
                failed.append(name)

    print("finished failed=" + ",".join(failed), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
