#!/usr/bin/env bash
set -u
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
write_status() { python3 - "$1" "$2" <<'PY_STATUS'
import json,sys,datetime
from pathlib import Path
Path('run_status.json').write_text(json.dumps({'phase':sys.argv[1],'exit_code':int(sys.argv[2]),'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n')
PY_STATUS
}
trap 'write_status interrupted 143; exit 143' TERM INT
write_status running 0
wolframscript -file "$PWD/run.wl" > run.log 2>&1
rc=$?
if [[ $rc -ne 0 || ! -f completion.json ]]; then
  write_status reduction_failed "$rc"
  exit 1
fi
write_status collecting 0
wolframscript -file "$PWD/collect.wl" > collect.log 2>&1
rc=$?
if [[ $rc -eq 0 && -f system_sizes.json ]]; then
  write_status complete 0
else
  write_status collection_failed "$rc"
  exit 1
fi
