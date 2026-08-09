# SSH Remote Operation Notes

This note records the working rules for operating the remote server `mdm@tianhe`
from local host A. The goal is to avoid repeated SSH failures caused by quoting,
long blocking commands, and unclear job-state checks.

## Fixed Current Method

Use only:

```bash
pipeline/tools/remote_tianhe.py ...
```

The helper uses the local SSH alias `tianhe` with `BatchMode=yes`.
This is the currently verified route. Do not use ad hoc password/pexpect SSH for
normal work.

Verified modes:

```bash
pipeline/tools/remote_tianhe.py run 'echo __OK__; hostname; date'
pipeline/tools/remote_tianhe.py put LOCAL REMOTE
pipeline/tools/remote_tianhe.py get REMOTE LOCAL
pipeline/tools/remote_tianhe.py run 'squeue -u mdm | head'
```

All four succeeded after the server-side transient issue cleared.

Operational rule: if the helper fails, diagnose that helper. Do not switch to a
new SSH method unless the failure is understood and recorded here.

The password/pexpect observations below were made during a server-side SSH
fluctuation and should not be used as the default path.

## Known Failure Modes

1. Complex inline SSH commands are fragile.
   Nested quotes, awk/python one-liners, here-docs, and shell substitutions have
   repeatedly been mangled by local shell parsing before reaching the remote.

2. Long `sleep ...; ssh ...` commands are unreliable in the local tool session.
   The local command session may disappear before the sleep finishes, even though
   the intended remote check has not run.
   If this happens, check for leftover local SSH processes with:
   `pgrep -af 'ssh .*mdm@tianhe'`.
   Kill stale checks before retrying; an old hung SSH process was observed after
   a failed status check.

3. SSH sometimes fails at handshake with:
   `kex_exchange_identification: Connection closed by remote host`.
   This is a transient connection/login-node issue. It should be retried later,
   not treated as proof that the remote command failed.
   Observed again while checking DP reduction status; the correct response is
   to wait and retry with the same `timeout 30s ssh -o ConnectTimeout=20 ...`
   pattern, not to change the remote command.

4. A file existing is not proof that a remote computation completed.
   For series generation, the reliable completion marker is the return-code file,
   not the presence of `series_prod`.

## Required Practice

### Simple Checks

For simple read-only checks, use short SSH commands with an explicit timeout:

```bash
timeout 30s ssh -o ConnectTimeout=20 mdm@tianhe 'hostname; date'
```

Keep these commands short. Do not embed large scripts or complicated quoting.
Use the outer `timeout`; `ConnectTimeout` alone does not cover every hang mode
after the TCP connection is established.

### Complex Remote Logic

For anything involving Python, loops, awk, grep pipelines with many quotes, or
multi-step logic:

1. Write a local script under `/tmp`.
2. Copy it to the remote with `scp`.
3. Run the script remotely with a short SSH command.

Example:

```bash
cat > /tmp/check_job.py <<'PY'
from pathlib import Path
print(Path("/fs2/home/mdm/Projects/series2d").exists())
PY
scp /tmp/check_job.py mdm@tianhe:/tmp/check_job.py
timeout 30s ssh -o ConnectTimeout=20 mdm@tianhe 'python3 /tmp/check_job.py'
```

This is the default pattern. Inline Python in SSH commands should be avoided.

### Authentication Hang

If a short SSH command hangs even with an outer `timeout`, diagnose once with
`ssh -vvv`. A repeated failure was traced to authentication, not to the remote
command: the log reached

```text
Offering public key: /home/mudaoming/.ssh/id_rsa ... agent
we sent a publickey packet, wait for reply
```

In this state:

1. Kill stale local SSH processes:

   ```bash
   pgrep -af 'ssh .*mdm@tianhe|ssh .*192\.168\.10\.50'
   kill <stale-pids>
   ```

2. Avoid the default noninteractive publickey/GSSAPI path.

3. If `sshpass` is unavailable, use a PTY session and explicit password
   authentication:

   ```bash
   ssh -tt -o GSSAPIAuthentication=no -o PreferredAuthentications=password \
     -o PubkeyAuthentication=no mdm@tianhe
   ```

   Then run remote commands inside that session.

### Slurm Job Checks

Do not use long local sleeps combined with SSH as the only state check. If a
periodic check is needed, run one short check per turn:

```bash
timeout 30s ssh -o ConnectTimeout=20 mdm@tianhe 'squeue -j JOBID'
```

If the job is absent from `squeue`, inspect the run directory:

```bash
timeout 30s ssh -o ConnectTimeout=20 mdm@tianhe \
  'cd /path/to/run && cat return_code 2>/dev/null; tail -80 run.err 2>/dev/null'
```

### Completion Markers

For `expand` series jobs, a cut is complete only if:

```text
series_prod.return_code exists and contains 0
```

If `series_prod.return_code` is missing, the job did not finish cleanly or is
still running. A nonempty `series_prod` file may be a partial write.

For `bl_sector_reduce` jobs, check:

```text
return_code
run.err
output
```

Do not report success until `return_code=0` and `output` has no `status=failed`.

## Things To Avoid

- Do not put large Python snippets directly inside `ssh 'python3 - <<PY ...'`.
- Do not rely on `series_prod` file existence as completion.
- Do not use a single path of sectors as a "closed chain" unless all higher
  sectors containing the current sector, in the sectormap sense, are included.
- Do not repeatedly retry a failing SSH handshake in a tight loop; wait and retry
  once with `ConnectTimeout=20`.
- Do not run plain `ssh ...` checks without an outer `timeout`; one check hung
  for more than the requested tool wait even with `ConnectTimeout=20`.
- Do not keep retrying the default SSH authentication path after it hangs at
  publickey authentication. Switch to explicit password authentication in a PTY.

### PTY Login Without Prompt

A PTY password login can authenticate and print the NSCC welcome banner but still
not provide a usable shell prompt. This was observed with:

```bash
ssh -tt -o GSSAPIAuthentication=no -o PreferredAuthentications=password \
  -o PubkeyAuthentication=no mdm@tianhe
```

If this happens, kill the local SSH process. Do not keep sending commands into
that dead PTY. Since `sshpass` is not installed locally but Python `pexpect` is
available, use a small `pexpect` wrapper for single remote commands with explicit
password authentication.

The same explicit authentication options must be used for `scp`. A failed upload
attempt hung because `scp` was called without these options and therefore fell
back to the same broken publickey/GSSAPI path. Use:

```bash
scp -o GSSAPIAuthentication=no -o PreferredAuthentications=password \
  -o PubkeyAuthentication=no LOCAL mdm@tianhe:REMOTE
```

### Process Cleanup

Do not use a broad `pkill -f '...ssh pattern...'` from a shell command whose own
command line contains the same pattern. It can kill the cleanup shell itself.
Use `pgrep -af` first, inspect the PIDs, then kill only the stale child PIDs.

## Current Diagnosis

The SSH problem must be treated separately from the reduction workflow. The
remote job submission should not be attempted until a minimal command such as
`echo __OK__` succeeds repeatedly.

Observed facts:

1. Default noninteractive SSH can hang at publickey authentication:

   ```text
   Offering public key: /home/mudaoming/.ssh/id_rsa ... agent
   we sent a publickey packet, wait for reply
   ```

2. `pexpect` with explicit password can reach the password prompt, send the
   password, and then time out without command output.

3. `SSH_ASKPASS` can authenticate and reach `exec request accepted`, but the
   trivial command `echo __OK__` still produced no output before timeout in one
   test. This points to a remote command/session startup issue, not just password
   entry.

4. After repeated attempts, the login node often rejects connections immediately:

   ```text
   kex_exchange_identification: Connection closed by remote host
   Connection closed by 192.168.10.50 port 22
   ```

5. Because of these states, continuing to retry business commands makes the
   situation worse. The next action after KEX rejection is to stop retrying,
   clear local stale processes, and wait for a clean minimal SSH diagnostic.

Required recovery test before any further remote work:

```bash
echo __OK__
```

must return successfully through the chosen wrapper. Only after that should
`put`, `submit`, or job-status commands be run.
