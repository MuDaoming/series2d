#!/usr/bin/env python3
from pathlib import Path
import importlib.util
import threading
import time
spec=importlib.util.spec_from_file_location('runner',Path(__file__).with_name('run_reconstruction.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
lock=threading.Lock();events=[];active=0;peak=0

def task(n):
    global active,peak
    with lock:
        active+=1;peak=max(peak,active);events.append(('start',n,time.monotonic()))
    time.sleep(0.25 if n==0 else 0.025)
    with lock:
        active-=1;events.append(('end',n,time.monotonic()))
    return n

assert sorted(module.dispatch(list(range(7)),3,task))==list(range(7))
assert peak==3
assert len([e for e in events if e[0]=='start'])==7
assert next(t for k,n,t in events if (k,n)==('start',3))<next(t for k,n,t in events if (k,n)==('end',0))
print('Dynamic dispatcher: capacity, exactly once, short jobs proceed while long job runs PASS')

# The actual worker input path is a frozen copy, not the subsequently edited source.
import tempfile
with tempfile.TemporaryDirectory() as directory:
    work=Path(directory);source=work/'source.json';source.write_text('{"x":1}')
    frozen=module.freeze_input(work,source,'family',module.sha(source))
    source.write_text('{"x":2}')
    assert frozen.read_text()=='{"x":1}'
print('Frozen runtime input survives source edits PASS')
