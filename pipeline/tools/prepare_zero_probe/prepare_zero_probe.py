#!/usr/bin/env python3
import argparse
import concurrent.futures
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline_driver"))
import pipeline_driver as base


def write_lines(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")


def copy_probe_config(base_config, out_config, bc_index, bc_count, deg, sector_map_rel):
    lines = []
    with open(base_config) as f:
        for line in f:
            if line.startswith("bc = "):
                bc = [0] * bc_count
                bc[bc_index] = 1
                lines.append("bc = " + base.vec_text(bc) + "\n")
            elif line.startswith("deg = "):
                lines.append(f"deg = {deg}\n")
            elif line.startswith("reduceMode = "):
                lines.append("reduceMode = cut\n")
            elif line.startswith("sectorMap = "):
                lines.append(f"sectorMap = {sector_map_rel}\n")
            else:
                lines.append(line)
    if not any(x.startswith("sectorMap = ") for x in lines):
        lines.append(f"sectorMap = {sector_map_rel}\n")
    out_config.parent.mkdir(parents=True, exist_ok=True)
    with open(out_config, "w") as f:
        f.writelines(lines)


def run_one(root, runner, s_path, cut_dir):
    log_path = cut_dir / "run.log"
    cmd = [
        str(runner),
        str(s_path),
        str(cut_dir / "config"),
        str(cut_dir / "target"),
        str(cut_dir / "series"),
    ]
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, cwd=root, stdout=log, stderr=subprocess.STDOUT)
    return cut_dir.name, proc.returncode


def parse_series(target_path, series_path):
    out = {}
    with open(target_path) as tf, open(series_path) as sf:
        for t, s in zip(tf, sf):
            t = t.strip()
            if t:
                out[t] = s.rstrip("\n")
    return out


def is_zero_series(text):
    return all(int(x) == 0 for x in re.findall(r"[+-]?\d+", text))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--family", required=True)
    ap.add_argument("--S", required=True)
    ap.add_argument("--base-config", required=True)
    ap.add_argument("--sector-table", required=True)
    ap.add_argument("--sectormap", required=True)
    ap.add_argument("--max-dot", type=int, default=2)
    ap.add_argument("--deg", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    run_root = root / "pipeline" / "runs" / args.family
    rows = base.read_sector_table(root / args.sector_table)
    maps, _, _ = base.read_sectormap(root / args.sectormap)
    cut_rows = [r for r in rows if r["has_cut"]]
    cut_sectors = [r["sector"] for r in cut_rows]
    cut_set = set(cut_sectors)

    meta_dir = run_root / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)
    write_lines(meta_dir / "cut_sectors", [base.vec_text(s) for s in cut_sectors])

    noncase_routes = base.build_noncase_routes(rows, maps)
    write_lines(
        meta_dir / "noncase_routes",
        [
            f"{base.vec_text(k)} -> " + " ".join(base.vec_text(x) for x in v)
            for k, v in sorted(noncase_routes.items(), key=lambda kv: (-sum(kv[0]), kv[0]))
        ],
    )

    probe_root = run_root / "zero_probe"
    runner = root / "expand" / "tools" / "fi_pipeline_runner" / "fi_pipeline_runner"
    s_path = root / args.S
    base_config = root / args.base_config

    summary_lines = ["sector\tnu_count\tkey_count\tseries_ready"]
    work = []
    for ci, sec in enumerate(cut_sectors):
        name = "".join(map(str, sec))
        cut_dir = probe_root / name
        nus = list(base.dot_nus_for_sector(sec, args.max_dot))
        target = []
        for nu in nus:
            target.extend(base.all_tags(nu))
        write_lines(cut_dir / "target", target)
        copy_probe_config(
            base_config,
            cut_dir / "config",
            ci,
            len(cut_sectors),
            args.deg,
            "../../symmetry/sectormap",
        )
        ready = (cut_dir / "series").exists()
        summary_lines.append(f"{base.vec_text(sec)}\t{len(nus)}\t{len(target)}\t{1 if ready else 0}")
        if args.run and not ready:
            work.append(cut_dir)
    write_lines(meta_dir / "zero_probe_key_counts", summary_lines)

    if work:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(run_one, root, runner, s_path, d) for d in work]
            failed = []
            for fut in concurrent.futures.as_completed(futs):
                name, rc = fut.result()
                print(f"zero_probe {name} rc={rc}", flush=True)
                if rc != 0:
                    failed.append(name)
            if failed:
                raise RuntimeError("zero probe failed: " + ", ".join(failed))

    zero_routes = []
    zero_summary = ["sector\tnu_count\tkey_count\tzero_nu_count\tzero_key_count"]
    for sec in cut_sectors:
        name = "".join(map(str, sec))
        cut_dir = probe_root / name
        series_path = cut_dir / "series"
        if not series_path.exists():
            continue
        series_map = parse_series(cut_dir / "target", series_path)
        nus = list(base.dot_nus_for_sector(sec, args.max_dot))
        zero_nus = []
        for nu in nus:
            labels = base.all_tags(nu)
            if all(label in series_map and is_zero_series(series_map[label]) for label in labels):
                routes = base.first_lower_cut_routes(base.sector_of(nu), cut_set, maps)
                if routes:
                    zero_nus.append(nu)
                    zero_routes.append(
                        f"cut={base.vec_text(sec)} nu={base.vec_text(nu)} -> "
                        + " ".join(base.vec_text(x) for x in routes)
                    )
        zero_summary.append(
            f"{base.vec_text(sec)}\t{len(nus)}\t{len(nus) * len(base.HEADS)}"
            f"\t{len(zero_nus)}\t{len(zero_nus) * len(base.HEADS)}"
        )
    write_lines(meta_dir / "zero_routes", zero_routes)
    write_lines(meta_dir / "zero_probe_summary", zero_summary)

    print(f"cut_sectors={len(cut_sectors)}")
    print(f"noncase_routes={len(noncase_routes)}")
    print(f"zero_routes={len(zero_routes)}")


if __name__ == "__main__":
    main()
