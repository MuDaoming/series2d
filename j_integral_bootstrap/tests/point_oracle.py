"""Fresh existing C++ J-IBP samples with explicit field and basis checks."""
import hashlib
import json
from pathlib import Path
import subprocess
from time import perf_counter


class NativeOracle:
    def __init__(self, binary, config, masters, prime, t0, work):
        self.work = Path(work)
        self.work.mkdir(parents=True, exist_ok=True)
        self.binary, self.config = str(Path(binary).resolve()), str(Path(config).resolve())
        self.masters, self.prime, self.t0 = masters, prime, t0
        self.seconds, self.calls = 0., 0
        self.cache = {}
        self.basis = self.work/'basis.json'
        identity = dict(binary_sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
                        config_sha256=hashlib.sha256(Path(config).read_bytes()).hexdigest(),
                        masters=masters, prime=prime, t0=t0)
        manifest = self.work/'identity.json'
        if manifest.exists() and json.loads(manifest.read_text()) != identity:
            raise ValueError('native oracle work directory identity changed')
        manifest.write_text(json.dumps(identity, indent=2)+'\n')
        self.basis.write_text(json.dumps(masters)+'\n')

    def at(self, tau):
        t = (tau+self.t0) % self.prime
        if t in self.cache:
            return self.cache[t]
        output = self.work/f'point_{t}.json'
        start = perf_counter()
        if not output.exists():
            command = [self.binary, self.config, str(output), '--point', str(t),
                       '--basis-file', str(self.basis), '--values', '1', '--compact', '1']
            with (self.work/f'point_{t}.log').open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=600)
            self.calls += 1
            self.seconds += perf_counter()-start
        data = json.loads(output.read_text())
        if (data['prime'] != self.prime or int(data['t_point']) != t
                or data['masters_at_point'] != self.masters or data['unclosed_derivatives']):
            raise ValueError('native sample field, point, basis or closure mismatch')
        self.cache[t] = data
        return data
