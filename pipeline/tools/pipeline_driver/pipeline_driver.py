#!/usr/bin/env python3
import argparse
import itertools
import os
import re
import shutil
import subprocess
from pathlib import Path

HEADS = [
    "FI",
    "BFI[XU]",
    "BFI[XD]",
    "BFI[YU]",
    "BFI[YD]",
    "BBFI[XU,YU]",
    "BBFI[XU,YD]",
    "BBFI[XD,YU]",
    "BBFI[XD,YD]",
]


def parse_vec(text):
    m = re.search(r"\{([^{}]*)\}", text)
    if not m:
        raise ValueError(f"no vector in {text!r}")
    body = m.group(1).strip()
    if not body:
        return []
    return tuple(int(x.strip()) for x in body.split(","))


def vec_text(v):
    return "{" + ",".join(str(x) for x in v) + "}"


def sector_of(nu):
    return tuple(1 if x > 0 else 0 for x in nu)


def tag(head, nu):
    return f"{head}{vec_text(nu)}"


def all_tags(nu):
    return [tag(h, nu) for h in HEADS]


def idx(sector):
    out = 0
    for b in sector:
        out = (out << 1) | b
    return out


def contains(big, small):
    return all((not s) or b for b, s in zip(big, small))


def dot_nus_for_sector(sector, max_dot):
    active = [i for i, b in enumerate(sector) if b]
    base = [1 if b else 0 for b in sector]
    out = []
    for total in range(max_dot + 1):
        for dots in dot_compositions(len(active), total):
            nu = base[:]
            for pos, d in zip(active, dots):
                nu[pos] += d
            out.append(tuple(nu))
    return out


def dot_compositions(n, total):
    if n == 0:
        if total == 0:
            yield ()
        return
    if n == 1:
        yield (total,)
        return
    for first in range(total, -1, -1):
        for rest in dot_compositions(n - 1, total - first):
            yield (first,) + rest


def read_sector_table(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            row["sector"] = parse_vec(row["sector"])
            row["representative"] = parse_vec(row["representative"])
            row["case"] = int(row["case"])
            row["is_representative"] = row["is_representative"] == "1"
            row["has_cut"] = row["has_cut"] == "1"
            row["prop_count"] = int(row["prop_count"])
            rows.append(row)
    return rows


def read_sectormap(path):
    maps = {}
    reverse = {}
    prop_maps = {}
    pat = re.compile(r"\{\s*(\{[^{}]*\})\s*->\s*(\{[^{}]*\})\s*,")
    pair_pat = re.compile(r"(\d+)\s*->\s*(\d+)")
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = pat.match(line)
            if not m:
                raise ValueError(f"bad sectormap line: {line}")
            src = parse_vec(m.group(1))
            dst = parse_vec(m.group(2))
            maps[src] = dst
            reverse.setdefault(dst, []).append(src)
            pairs = {}
            for a, b in pair_pat.findall(line):
                pairs[int(a) - 1] = int(b) - 1
            prop_maps[src] = pairs
    return maps, reverse, prop_maps


def canonical_sector(sector, maps):
    seen = set()
    cur = sector
    while cur in maps:
        if cur in seen:
            raise ValueError(f"cycle in sectormap at {cur}")
        seen.add(cur)
        cur = maps[cur]
    return cur


def canonical_nu(nu, prop_maps):
    cur = tuple(nu)
    seen = set()
    while sector_of(cur) in prop_maps:
        sec = sector_of(cur)
        if sec in seen:
            raise ValueError(f"cycle in sectormap at {sec}")
        seen.add(sec)
        mapped = [0] * len(cur)
        for src, dst in prop_maps[sec].items():
            mapped[dst] = cur[src]
        cur = tuple(mapped)
    return cur


def read_series_map(series_path, target_path):
    data = {}
    with open(target_path) as tf, open(series_path) as sf:
        for t, s in zip(tf, sf):
            t = t.strip()
            if not t:
                continue
            data[t] = s.rstrip("\n")
    return data


def write_lines(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")


def copy_config_template(base_config, out_config, bc_index, bc_count, sectormap_rel):
    lines = []
    with open(base_config) as f:
        for line in f:
            if line.startswith("bc = "):
                bc = [0] * bc_count
                bc[bc_index] = 1
                lines.append("bc = " + vec_text(bc) + "\n")
            elif line.startswith("sectorMap = "):
                lines.append(f"sectorMap = {sectormap_rel}\n")
            elif line.startswith("reduceMode = "):
                lines.append("reduceMode = cut\n")
            else:
                lines.append(line)
    if not any(x.startswith("sectorMap = ") for x in lines):
        lines.append(f"sectorMap = {sectormap_rel}\n")
    out_config.parent.mkdir(parents=True, exist_ok=True)
    with open(out_config, "w") as f:
        f.writelines(lines)


def run(cmd, cwd):
    print("+", " ".join(str(x) for x in cmd), flush=True)
    subprocess.run([str(x) for x in cmd], cwd=cwd, check=True)


def build_noncase_routes(rows, maps):
    valid = {r["sector"]: r for r in rows}
    cut_sectors = {r["sector"] for r in rows if r["has_cut"]}
    routes = {}
    for r in rows:
        if r["case"] == 0:
            continue
        source = r["sector"]
        active = [i for i, b in enumerate(source) if b]
        found = []
        for depth in range(1, len(active) + 1):
            layer = []
            for remove in itertools.combinations(active, depth):
                sub = list(source)
                for i in remove:
                    sub[i] = 0
                sub = tuple(sub)
                rep = canonical_sector(sub, maps)
                if rep in cut_sectors:
                    layer.append(rep)
            if layer:
                found = sorted(set(layer), key=lambda x: (-sum(x), x))
                break
        if found:
            routes[source] = found
    return routes


def first_lower_cut_routes(source, cut_sectors, maps):
    active = [i for i, b in enumerate(source) if b]
    for depth in range(1, len(active) + 1):
        layer = []
        for remove in itertools.combinations(active, depth):
            sub = list(source)
            for i in remove:
                sub[i] = 0
            rep = canonical_sector(tuple(sub), maps)
            if rep in cut_sectors:
                layer.append(rep)
        if layer:
            return sorted(set(layer), key=lambda x: (-sum(x), x))
    return []


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
    ap.add_argument("--old-series-dir", required=True)
    ap.add_argument("--old-base-target", required=True)
    ap.add_argument("--max-dot", type=int, default=2)
    ap.add_argument("--search-deg", type=int, default=1000)
    ap.add_argument("--search-m", type=int, default=10)
    ap.add_argument("--ncheck", type=int, default=1)
    ap.add_argument("--delta", default="571")
    ap.add_argument("--run-expand", action="store_true")
    ap.add_argument("--run-search", action="store_true")
    ap.add_argument("--run-bl", action="store_true")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    run_root = root / "pipeline" / "runs" / args.family
    rows = read_sector_table(root / args.sector_table)
    maps, reverse_maps, prop_maps = read_sectormap(root / args.sectormap)
    cut_rows = [r for r in rows if r["has_cut"]]
    cut_sectors = [r["sector"] for r in cut_rows]
    cut_set = set(cut_sectors)
    all_valid_sectors = [r["sector"] for r in rows]

    meta_dir = run_root / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)
    write_lines(meta_dir / "cut_sectors", [vec_text(s) for s in cut_sectors])

    noncase_routes = build_noncase_routes(rows, maps)
    write_lines(
        meta_dir / "noncase_routes",
        [f"{vec_text(k)} -> " + " ".join(vec_text(x) for x in v)
         for k, v in sorted(noncase_routes.items(), key=lambda kv: (-sum(kv[0]), kv[0]))],
    )

    # Generate the canonical full representative target used by old full-series files.
    rep_sectors = cut_sectors
    rep_lines = []
    for sec in sorted(rep_sectors, key=idx, reverse=True):
        for nu in dot_nus_for_sector(sec, args.max_dot):
            rep_lines.extend(all_tags(nu))
    generated_base_target = run_root / "expand" / "target_dot2_representatives"
    write_lines(generated_base_target, rep_lines)

    old_base_target = root / args.old_base_target
    with open(generated_base_target) as a, open(old_base_target) as b:
        if a.read() != b.read():
            raise RuntimeError(
                f"generated base target differs from old target: {generated_base_target} vs {old_base_target}"
            )

    base_series_by_cut = {}
    for i, sec in enumerate(cut_sectors, start=1):
        series_path = root / args.old_series_dir / f"symmetry_bc{i}" / "series_dot2deg1000"
        if not series_path.exists():
            raise FileNotFoundError(series_path)
        base_series_by_cut[sec] = read_series_map(series_path, old_base_target)

    # Include pre-existing appended series only as series cache, not as pipeline decisions.
    for i, sec in enumerate(cut_sectors, start=1):
        d = root / args.old_series_dir / f"symmetry_bc{i}"
        t = d / "target_appended"
        s = d / "series_appended"
        if t.exists() and s.exists():
            base_series_by_cut[sec].update(read_series_map(s, t))

    # Zero probe from available full series.
    zero_routes = {}
    zero_report = []
    for sec in cut_sectors:
        series_map = base_series_by_cut[sec]
        for nu in dot_nus_for_sector(sec, args.max_dot):
            labels = all_tags(nu)
            if all(label in series_map and is_zero_series(series_map[label]) for label in labels):
                routes = first_lower_cut_routes(sector_of(nu), cut_set, maps)
                if routes:
                    zero_routes[(sec, nu)] = routes
                    zero_report.append(
                        f"cut={vec_text(sec)} nu={vec_text(nu)} -> "
                        + " ".join(vec_text(x) for x in routes)
                    )
    write_lines(meta_dir / "zero_routes", zero_report)

    # Build per-cut targets and series, filling missing series on demand.
    fi_runner = root / "expand/tools/fi_pipeline_runner/fi_pipeline_runner"
    for ci, sec in enumerate(cut_sectors):
        out_dir = run_root / "expand" / "".join(map(str, sec))
        out_dir.mkdir(parents=True, exist_ok=True)
        copy_config_template(
            root / args.base_config,
            out_dir / "config",
            ci,
            len(cut_sectors),
            "../../symmetry/sectormap",
        )

        equivalents = [sec] + reverse_maps.get(sec, [])
        target_sectors = []
        for candidate in all_valid_sectors:
            if any(contains(candidate, eq) for eq in equivalents):
                target_sectors.append(candidate)
        target_lines = []
        for tsec in sorted(set(target_sectors), key=idx, reverse=True):
            for nu in dot_nus_for_sector(tsec, args.max_dot):
                target_lines.extend(all_tags(nu))
        target_lines = list(dict.fromkeys(target_lines))
        write_lines(out_dir / "target", target_lines)

        series_map = dict(base_series_by_cut[sec])
        if (out_dir / "target_missing").exists() and (out_dir / "series_missing").exists():
            series_map.update(read_series_map(out_dir / "series_missing", out_dir / "target_missing"))
        missing = [x for x in target_lines if x not in series_map]
        if missing:
            write_lines(out_dir / "target_missing", missing)
            if args.run_expand:
                run([fi_runner, root / args.S, out_dir / "config",
                     out_dir / "target_missing", out_dir / "series_missing"], root)
                series_map.update(read_series_map(out_dir / "series_missing", out_dir / "target_missing"))
            else:
                print(f"missing series for {vec_text(sec)}: {len(missing)}")
        unresolved = [x for x in target_lines if x not in series_map]
        if unresolved:
            raise RuntimeError(f"unresolved series for {vec_text(sec)}: {len(unresolved)}")
        write_lines(out_dir / "series", [series_map[x] for x in target_lines])

    # Build search inputs.
    p_value = None
    n_value = None
    with open(root / args.base_config) as f:
        for line in f:
            if line.startswith("p ="):
                p_value = line.split("=", 1)[1].strip()
            if line.startswith("N ="):
                n_value = line.split("=", 1)[1].strip()
    for sec in cut_sectors:
        sdir = run_root / "search" / "".join(map(str, sec))
        sdir.mkdir(parents=True, exist_ok=True)
        g_lines = []
        for nu in dot_nus_for_sector(sec, args.max_dot):
            g_lines.extend(all_tags(nu))
        for source, routes in noncase_routes.items():
            if sec in routes:
                for nu in dot_nus_for_sector(source, args.max_dot):
                    g_lines.extend(all_tags(nu))
        for (cut, nu), routes in zero_routes.items():
            if sec in routes:
                g_lines.extend(all_tags(nu))
        g_lines = list(dict.fromkeys(g_lines))
        write_lines(sdir / "G", g_lines)
        with open(sdir / "config", "w") as f:
            f.write(f"N = {n_value}\n")
            f.write(f"deg = {args.search_deg}\n")
            f.write(f"m = {args.search_m}\n")
            f.write(f"ncheck = {args.ncheck}\n")
            f.write(f"p = {p_value}\n")
        series_abs = run_root / "expand" / "".join(map(str, sec)) / "series"
        target_abs = run_root / "expand" / "".join(map(str, sec)) / "target"
        write_lines(sdir / "series_list", [f"{series_abs} {target_abs}"])

    if args.run_search:
        poly = root / "search/tools/poly_relation_searcher/poly_relation_searcher"
        solver = root / "search/tools/integral_solver/integral_solver"
        for sec in cut_sectors:
            sdir = run_root / "search" / "".join(map(str, sec))
            run([poly, sdir / "config", sdir / "G", sdir / "series_list", sdir / "poly_relation"], root)
            run([solver, sdir / "G", sdir / "poly_relation", args.delta, sdir / "integral_solution"], root)
            extract_masters(sdir / "integral_solution", sdir / "masters")

    # Build BL inputs.  The current reducer subtracts all more-complex sector
    # contributions, so each lower sector needs series for masters from all
    # more-complex sectors as well.
    bdir = run_root / "bl_sector_reduce"
    bdir.mkdir(parents=True, exist_ok=True)
    with open(bdir / "config", "w") as f:
        f.write(f"N = {n_value}\n")
        f.write(f"deg = {args.search_deg}\n")
        f.write("m = 1000\n")
        f.write(f"p = {p_value}\n")
        f.write("K_safety = 10\n")
        f.write("K_cert = 10\n")
    master_by_sector = {}
    for sec in cut_sectors:
        mpath = run_root / "search" / "".join(map(str, sec)) / "masters"
        if mpath.exists():
            with open(mpath) as f:
                master_by_sector[sec] = [line.strip() for line in f if line.strip() and not line.startswith("#")]
        else:
            master_by_sector[sec] = []

    fi_runner = root / "expand/tools/fi_pipeline_runner/fi_pipeline_runner"
    lines = []
    for sec in cut_sectors:
        name = "".join(map(str, sec))
        exp_dir = run_root / "expand" / name
        search_dir = run_root / "search" / name
        bl_exp_dir = run_root / "bl_series" / name
        bl_exp_dir.mkdir(parents=True, exist_ok=True)

        series_map = read_series_map(exp_dir / "series", exp_dir / "target")
        if (bl_exp_dir / "target_missing").exists() and (bl_exp_dir / "series_missing").exists():
            series_map.update(read_series_map(bl_exp_dir / "series_missing",
                                              bl_exp_dir / "target_missing"))
        target_lines = list(series_map.keys())
        for higher, masters in master_by_sector.items():
            if sum(higher) <= sum(sec):
                continue
            for m in masters:
                if m not in series_map:
                    target_lines.append(m)
        target_lines = list(dict.fromkeys(target_lines))
        missing = [x for x in target_lines if x not in series_map]
        if missing:
            write_lines(bl_exp_dir / "target_missing", missing)
            run([fi_runner, root / args.S, exp_dir / "config",
                 bl_exp_dir / "target_missing", bl_exp_dir / "series_missing"], root)
            series_map.update(read_series_map(bl_exp_dir / "series_missing", bl_exp_dir / "target_missing"))
        write_lines(bl_exp_dir / "target", target_lines)
        write_lines(bl_exp_dir / "series", [series_map[x] for x in target_lines])

        lines.append(
            f"sector={vec_text(sec)} "
            f"{os.path.relpath(bl_exp_dir / 'series', bdir)} "
            f"{os.path.relpath(bl_exp_dir / 'target', bdir)} "
            f"{os.path.relpath(search_dir / 'masters', bdir)}"
        )
    write_lines(bdir / "sector_series_list", lines)
    dot1_objects = []
    for sec in sorted(all_valid_sectors, key=idx, reverse=True):
        for nu in dot_nus_for_sector(sec, 1):
            dot1_objects.append(tag("FI", canonical_nu(nu, prop_maps)))
    dot1_objects = list(dict.fromkeys(dot1_objects))
    write_lines(bdir / "object_list_dot1", dot1_objects)

    if args.run_bl:
        reducer = root / "bl_sector_reduce/tools/bl_sector_reducer/bl_sector_reducer"
        run([reducer, bdir / "config", bdir / "sector_series_list",
             bdir / "reduction_dot1", bdir / "object_list_dot1"], root)


def extract_masters(solution_path, masters_path):
    masters = []
    in_mis = False
    with open(solution_path) as f:
        for line in f:
            s = line.strip()
            if s == "#MIs":
                in_mis = True
                continue
            if in_mis and s.startswith("["):
                break
            if in_mis and s and not s.startswith("#"):
                masters.extend(re.findall(r"(?:FI|BFI\[[^\]]+\]|BBFI\[[^\]]+\])\{[^}]+\}", s))
    write_lines(masters_path, masters)


if __name__ == "__main__":
    main()
