# Tianhe Server Notes

This records the job-system assumptions for running `series2d` on server B
(`tianhe`, project path `/fs2/home/mdm/Projects/series2d`).

## Cluster Facts

- Scheduler: Slurm.
- Submit command: `sbatch`.
- Job-step command: `srun`.
- Queue/status commands: `squeue`, `sinfo`.
- Common partition in existing scripts: `cp6`.
- Debug partition seen in existing scripts: `debug6`.
- Max simultaneously submitted jobs for this account/workflow: 30.
- Charging is by allocated node. Each submitted Slurm job occupies and is charged
  as a node allocation.
- One node has 56 CPU cores and 256 GB memory.

## Existing Script Style

Reference files:

- `/fs2/home/mdm/Projects/test/template_TH.sh`
- `/fs2/home/mdm/Projects/Higgs+jet@NLO_EW/gg-gh/families/gg-gh-10/script/*.sh`

Observed patterns:

- Scripts use `#!/bin/bash` and `#SBATCH` directives.
- They normally include `#SBATCH --get-user-env`.
- `template_TH.sh` uses:
  - `#SBATCH --nodes=2`
  - `#SBATCH --ntasks-per-node=54`
  - `#SBATCH --cpus-per-task=1`
  - `#SBATCH --partition=cp6`
  - one `srun bash -c '...'` block, with task index from `SLURM_PROCID`.
- Other project scripts often allocate several nodes and launch many independent
  one-core subtasks inside one Slurm job:
  - `srun --exclusive -N 1 -n 1 -c 1 bash -c "..." &`
  - short `sleep` between launches
  - final `wait`
- Existing Mathematica jobs call:
  - `bash /fs2/home/s522pkuphy/hep_packages/crack.sh`
  - `/fs2/home/s522pkuphy/Mathematica/bin/math -script ...`

## Implications For `series2d`

- Do not submit one Slurm job per small cut if many cuts are involved; this wastes
  node allocations and can hit the 30-job limit.
- Prefer batching multiple independent cuts inside one Slurm job and use up to 56
  cores on the node.
- For C++ `expand`/`search`/`reduce` tasks, one practical pattern is:
  - one node per batch,
  - up to 56 concurrent single-core `srun --exclusive -N 1 -n 1 -c 1` steps,
  - per-cut logs under the cut directory,
  - one batch-level log recording start/finish times.
- For memory-heavy expand jobs, choose batch size from the estimated memory:
  `cache_keys * deg^2 / 2 * 8 bytes`, plus overhead. Keep total concurrent memory
  comfortably below 256 GB per node.
- Because the cluster limit is 30 submitted jobs, batch scripts should split the
  115 dp cuts into at most 30 Slurm jobs, preferably fewer when each job can
  efficiently fill a node.
