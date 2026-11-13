"""Read-only snapshot of the two local background services."""
import datetime,json,subprocess
from pathlib import Path
base=Path(__file__).resolve().parent
now=datetime.datetime.now(datetime.timezone.utc)
record={'checked_utc':now.isoformat(),'services':{}}
for name,folder in [('dp-fi-native-dot2','fi/blade_native'),('dp-j-ibp-size','j')]:
    result=subprocess.run(['systemctl','--user','show',name+'.service','-p','ActiveState','-p','SubState','-p','Result','-p','MainPID','-p','MemoryCurrent','-p','MemoryPeak','-p','CPUUsageNSec'],capture_output=True,text=True)
    entry=dict(line.split('=',1) for line in result.stdout.splitlines() if '=' in line)
    log=base/folder/'run.log'
    if log.exists():
        with log.open('rb') as handle:
            handle.seek(max(0,log.stat().st_size-2000))
            entry['log_tail']=handle.read().decode(errors='replace')
    record['services'][name]=entry
out=base/'status_checks';out.mkdir(exist_ok=True)
path=out/(now.strftime('%Y%m%dT%H%M%SZ')+'.json')
path.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(record,indent=2,ensure_ascii=False))
