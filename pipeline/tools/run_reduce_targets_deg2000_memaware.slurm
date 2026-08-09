#!/bin/bash
#SBATCH --job-name=dp_red_t2000
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --partition=cp6
#SBATCH --get-user-env

set -u
ROOT=/fs2/home/mdm/Projects/series2d
RUN=$ROOT/pipeline/runs/dp
TARGET_NAME=reduce_targets_21111111_31111111
RUNNER=$ROOT/expand/tools/fi_pipeline_runner/fi_pipeline_runner
SFILE=$RUN/input/S
ESTIMATE=$RUN/metadata/reduce_targets_deg2000_estimate.tsv
JOB_ID=${SLURM_JOB_ID:-manual}
BATCH_LOG=$RUN/logs/reduce_targets_deg2000_${JOB_ID}.log
MEM_CAP_GIB=${MEM_CAP_GIB:-200}
MAX_PAR=${SLURM_NTASKS:-56}
DEG=2000

mkdir -p "$RUN/logs"
echo "job_id=$JOB_ID host=$(hostname) start=$(date -Is) max_parallel=$MAX_PAR mem_cap_gib=$MEM_CAP_GIB" > "$BATCH_LOG"
echo "root=$ROOT" >> "$BATCH_LOG"
echo "runner=$RUNNER" >> "$BATCH_LOG"
echo "target_name=$TARGET_NAME" >> "$BATCH_LOG"
echo "estimate=$ESTIMATE" >> "$BATCH_LOG"

export ROOT RUN TARGET_NAME RUNNER SFILE BATCH_LOG MEM_CAP_GIB MAX_PAR DEG ESTIMATE
python3 - <<'PY'
import csv
import os
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(os.environ["ROOT"])
RUN = pathlib.Path(os.environ["RUN"])
TARGET_NAME = os.environ["TARGET_NAME"]
RUNNER = os.environ["RUNNER"]
SFILE = os.environ["SFILE"]
BATCH_LOG = pathlib.Path(os.environ["BATCH_LOG"])
MEM_CAP = float(os.environ.get("MEM_CAP_GIB", "200"))
MAX_PAR = int(os.environ.get("MAX_PAR", "56"))
DEG = int(os.environ.get("DEG", "2000"))
ESTIMATE = pathlib.Path(os.environ["ESTIMATE"])


def log(msg):
    with BATCH_LOG.open("a") as f:
        f.write(msg + "\n")
        f.flush()


items = []
with ESTIMATE.open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        items.append({
            "cut": row["cut"],
            "mem": float(row["mem_gib_deg2000"]),
            "time": float(row["time_s_deg2000"]),
        })

items.sort(key=lambda item: item["time"], reverse=True)
waiting = items[:]
running = []
used_mem = 0.0
failed = 0
started = 0
finished = 0


def launch(item):
    global used_mem, started
    cut = item["cut"]
    cut_dir = RUN / "cuts" / cut / TARGET_NAME
    rcfile = cut_dir / f"series_deg{DEG}.return_code"
    outfile = cut_dir / f"series_deg{DEG}"
    if rcfile.exists() and rcfile.read_text().strip() == "0" and outfile.exists() and outfile.stat().st_size > 0:
        log(f"cut={cut} skip_existing_success")
        return None

    inner = f"""
set +e
ROOT={str(ROOT)!r}
RUNNER={RUNNER!r}
SFILE={SFILE!r}
dir={str(cut_dir)!r}
cd "$ROOT" || exit 99
awk 'BEGIN{{done=0}} /^deg[[:space:]]*=/ {{print "deg = {DEG}"; done=1; next}} {{print}} END{{if(!done) print "deg = {DEG}"}}' "$dir/config" > "$dir/config_deg{DEG}"
t0=$(date +%s)
echo "cut={cut} start=$(date -Is) host=$(hostname) deg={DEG}" >> {str(BATCH_LOG)!r}
"$RUNNER" "$SFILE" "$dir/config_deg{DEG}" "$dir/target" "$dir/series_deg{DEG}" > "$dir/series_deg{DEG}.run.log" 2>&1
rc=$?
t1=$(date +%s)
echo "$rc" > "$dir/series_deg{DEG}.return_code"
echo "$((t1-t0))" > "$dir/series_deg{DEG}.elapsed_sec"
echo "cut={cut} finish=$(date -Is) rc=$rc elapsed=$((t1-t0))" >> {str(BATCH_LOG)!r}
exit $rc
"""
    proc = subprocess.Popen(
        ["srun", "--exclusive", "-N", "1", "-n", "1", "-c", "1", "bash", "-lc", inner],
        stdin=subprocess.DEVNULL,
    )
    used_mem += item["mem"]
    started += 1
    log(
        f"cut={cut} launched pid={proc.pid} est_mem_gib={item['mem']:.3f} "
        f"est_time_s={item['time']:.1f} used_mem_gib={used_mem:.1f} running={len(running) + 1}"
    )
    return proc


while waiting or running:
    made_progress = True
    while made_progress and len(running) < MAX_PAR:
        made_progress = False
        for idx, item in enumerate(waiting):
            if used_mem + item["mem"] <= MEM_CAP:
                waiting.pop(idx)
                proc = launch(item)
                if proc is not None:
                    running.append((proc, item))
                made_progress = True
                break

    if not running:
        if waiting:
            log("ERROR no runnable item under memory cap")
            sys.exit(2)
        break

    time.sleep(5)
    still_running = []
    for proc, item in running:
        rc = proc.poll()
        if rc is None:
            still_running.append((proc, item))
            continue
        used_mem -= item["mem"]
        finished += 1
        if rc:
            failed += 1
            log(f"cut={item['cut']} process_failed rc={rc} used_mem_gib={used_mem:.1f} finished={finished}")
        else:
            log(f"cut={item['cut']} process_done rc=0 used_mem_gib={used_mem:.1f} finished={finished}")
    running = still_running

log(f"job_finish failed={failed} started={started} finished={finished} time={time.strftime('%Y-%m-%dT%H:%M:%S%z')}")
sys.exit(1 if failed else 0)
PY
rc=$?
echo "job_id=$JOB_ID finish=$(date -Is) rc=$rc" >> "$BATCH_LOG"
exit "$rc"
