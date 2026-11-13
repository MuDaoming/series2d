"""AB JSON-lines I/O only: saved anchor/diagonals or the public J-IBP executable.

No AB learning, ansatz construction, sector traversal, fitting or sample choice.
"""
import argparse,json,sys,hashlib,fcntl,subprocess,time,os,resource,gzip
from fractions import Fraction
from pathlib import Path

def load(path):
    path=Path(path)
    if path.suffix=='.gz':
        with gzip.open(path,'rt') as file:return json.load(file)
    return json.loads(path.read_text())
def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def send(data):print(json.dumps(data,separators=(',',':')),flush=True)
def write(path,data):
    temp=path.with_suffix('.tmp-'+str(os.getpid()));temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(path)
def event(directory,data):
    with (directory/'events.jsonl').open('a') as file:
        fcntl.flock(file,fcntl.LOCK_EX);file.write(json.dumps(dict(time=time.time(),pid=os.getpid(),**data))+'\n')
p=argparse.ArgumentParser();p.add_argument('capability',choices=['anchor','cut','numeric']);p.add_argument('--input',required=True);p.add_argument('--normalization');p.add_argument('--binary');p.add_argument('--cache');a=p.parse_args()
source=load(a.input);request=None
if a.capability=='numeric':
    cache=Path(a.cache);cache.mkdir(parents=True,exist_ok=True);(cache/'points').mkdir(exist_ok=True);(cache/'logs').mkdir(exist_ok=True);(cache/'receipts').mkdir(exist_ok=True)
    norm=load(a.normalization);q=Fraction(norm['t0']);basis=source['request']['basis'];prime=source['point']['prime'];center=q.numerator*pow(q.denominator,-1,prime)%prime;gauge=norm['gauge_power']
    identity=dict(request_sha256=digest(source),binary_sha256=hashlib.sha256(Path(a.binary).read_bytes()).hexdigest())
    # All sessions at the same setting share this full-point cache.
    with (cache/'identity.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (cache/'identity.json').exists():
            if load(cache/'identity.json')!=identity:raise ValueError('numeric cache identity differs')
        else:write(cache/'identity.json',identity)
for line in sys.stdin:
    message=json.loads(line)
    if a.capability=='anchor':
        if message['operation']!='anchor' or message['prime']!=source['prime'] or message['assignments']!=source['assignments'] or message['masters']!=source['masters'] or set(message['reductions'])!=set(source['reductions']):raise ValueError('anchor input identity/order differs')
        send(dict(request_sha256=digest(message),de=source['de'],reductions=source['reductions']))
    elif a.capability=='cut':
        if message['operation']!='maximal-cut' or message['prime']!=source['prime'] or message['assignments']!=source['assignments']:raise ValueError('diagonal setting differs')
        block=source['blocks'][message['sector']]
        if message['masters']!=block['masters']:raise ValueError('diagonal basis order differs')
        send(dict(request_sha256=digest(message),matrix=block['matrix']))
    elif request is None:
        if message['operation']!='numeric-prepare' or message['mode'] not in ('de','reduction') or message['variable']!=norm['variable'] or message['prime']!=prime or message['assignments']!={} or set(message['masters'])!=set(basis):raise ValueError('unsupported numerical request/setting')
        for target in message['targets']:
            if message['mode']=='de':
                if target['kind']!='derivative' or target['master'] not in basis:raise ValueError('unsupported DE target')
            elif target not in ['J{'+','.join(map(str,v))+'}' for v in source['request']['targets']]:raise ValueError('unsupported numerical target')
        request=message;key=digest(request);event(cache,dict(operation='prepare',mode=request['mode'],task=request['task'],target_count=len(request['targets'])));send(dict(request_sha256=key,ready=True))
    else:
        if message['operation']!='evaluate' or message['request_sha256']!=key:raise ValueError('numeric session identity differs')
        x=message['x'];t=(x+center)%prime
        if t in (0,1):
            event(cache,dict(operation='bad-point',task=request['task'],x=x,t=t));send(dict(request_sha256=key,x=x,bad_sample_point=True));continue
        output=cache/'points'/('point'+str(t)+'.json');receipt=cache/'receipts'/('point'+str(t)+'.json');hit=False
        with (cache/'points'/('point'+str(t)+'.lock')).open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if output.exists() and receipt.exists():
                saved=load(receipt)
                if not saved.get('validated') or saved['output_sha256']!=hashlib.sha256(output.read_bytes()).hexdigest():raise ValueError('cached point receipt differs')
                hit=True
            else:
                data=json.loads(json.dumps(source));data['point']['t']=str(t);data['request']['derivatives']='basis' if request['mode']=='de' else [];inp=cache/'points'/('input'+str(t)+'.json');write(inp,data)
                start=time.monotonic();cpu=resource.getrusage(resource.RUSAGE_CHILDREN)
                with (cache/'logs'/('point'+str(t)+'.log')).open('w') as log:code=subprocess.run([a.binary,'--input',str(inp),'--output',str(output)],stdout=log,stderr=log).returncode
                after=resource.getrusage(resource.RUSAGE_CHILDREN)
                record=dict(t=t,x=x,returncode=code,wall_seconds=time.monotonic()-start,cpu_seconds=after.ru_utime+after.ru_stime-cpu.ru_utime-cpu.ru_stime,validated=False,worker_pid=os.getpid(),input_sha256=digest(data))
                if code:write(receipt,record);event(cache,dict(operation='IBP-failed',**record));raise RuntimeError('J IBP failed; see '+str(cache/'logs'/('point'+str(t)+'.log')))
                values=load(output)
                if values['status']!='complete' or values['basis']!=basis or values['unclosed_derivatives'] :raise ValueError('IBP basis/closure/coverage differs')
                record.update(validated=True,output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),diagnostics=values['diagnostics']);write(receipt,record);event(cache,dict(operation='IBP-complete',mode=request['mode'],task=request['task'],**record))
            values=load(output)
        if request['mode']=='de':
            shift=gauge*(pow(t,-1,prime)+pow((t-1)%prime,-1,prime))%prime
            rows=[[(values['derivatives'][target['master']].get(m,0)+(shift if m==target['master'] else 0))%prime for m in request['masters']] for target in request['targets']]
        else:
            rows=[[values['reductions'][target].get(m,0) for m in request['masters']] for target in request['targets']]
        event(cache,dict(operation='evaluate',mode=request['mode'],task=request['task'],x=x,t=t,cache_hit=hit,target_count=len(request['targets'])))
        send(dict(request_sha256=key,x=x,values=rows))
