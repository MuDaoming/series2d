import json,subprocess,time,datetime
from pathlib import Path
out=Path(__file__).resolve().parent
root=out.parents[4]
started=time.monotonic()
with (out/'run.log').open('w') as log:
    result=subprocess.run([str(root/'j_integral/src/IBP/build/j-ibp'),'--input',str(out.parent/'input/j_common_prime.json'),'--output',str(out/'result.json')],stdout=log,stderr=log)
status={'exit_code':result.returncode,'seconds':time.monotonic()-started,'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
if result.returncode==0:
    data=json.loads((out/'result.json').read_text())
    status['diagnostics']=data.get('diagnostics',{})
    status['result_keys']=list(data)
(out/'statistics.json').write_text(json.dumps(status,indent=2)+'\n')
raise SystemExit(result.returncode)
