"""Ordinary-point jets of known FBI equations; never consumes a J equation."""
from __future__ import annotations

from math import comb
import time
import sympy as s
from sympy.polys.domains import GF
from periods import FixedPeriodReducer


class Model:
    def __init__(self, data):
        self.data = data
        self.p = int(data['prime'])
        if not s.isprime(self.p):
            raise ValueError('prime is not prime')
        self.t, self.z = s.symbols('t z')
        self.L = GF(self.p).frac_field(self.t, self.z)
        self.R = self.L.field.ring
        self.K = GF(self.p).frac_field(self.z)
        self.zz = self.K.gens[0]
        self.t0 = self.scalar(data['t0'])
        self.At = [[self.rational(x) for x in row] for row in data['At']]
        self.Az = [[self.rational(x) for x in row] for row in data['Az']]
        self.rank = len(self.At)
        if not self.rank or any(len(row) != self.rank for row in self.At+self.Az):
            raise ValueError('At and Az must be square matrices of equal size')
        if len(self.Az) != self.rank:
            raise ValueError('At/Az rank mismatch')
        self.targets = {label: [self.rational(x) for x in row]
                        for label, row in data['targets'].items()}
        if any(len(row) != self.rank for row in self.targets.values()):
            raise ValueError('target row length differs from FBI rank')
        endpoints = {int(k): v for k, v in data['endpoints'].items()}
        if set(endpoints) != {0, 1}:
            raise ValueError('Both endpoints must have explicit modes')
        q = self.R.one
        for x in [x for row in self.At for x in row]:
            q = q.lcm(x.denom)
        self.qt = self.poly_jet(q)
        if not self.qt[0]:
            raise ValueError('t0 is not ordinary for the FBI t recurrence')
        self.nt = [[self.poly_jet(self.polynomial(x*self.L.field.new(q)))
                    for x in row] for row in self.At]
        self.object_data = {}
        extra = [1/self.qt[0]]
        for label, row in self.targets.items():
            den = self.R.one
            for x in row:
                den = den.lcm(x.denom)
            dj = self.poly_jet(den)
            if not dj[0]:
                raise ValueError(f't0 is not ordinary for {label}')
            numer = [self.poly_jet(self.polynomial(x*self.L.field.new(den))) for x in row]
            self.object_data[label] = dj, numer
            extra.append(1/dj[0])
        az0 = [[self.specialize(x) for x in row] for row in self.Az]
        self.reducer = FixedPeriodReducer(self.p, az0, extra, endpoints)

    def scalar(self, value):
        q = s.Rational(str(value))
        return int(q.p) * pow(int(q.q), -1, self.p) % self.p

    def rational(self, value):
        if hasattr(value, 'field') and value.field == self.L.field:
            return value
        expr = s.sympify(value, locals={'t': self.t, 'z': self.z})
        if expr.free_symbols - {self.t, self.z}:
            raise ValueError(f'Unspecified symbols in input: {expr.free_symbols}')
        n, d = s.together(expr).as_numer_denom()
        return self.L.field.new(self.R.from_expr(n), self.R.from_expr(d))

    def polynomial(self, value):
        if any(any(monomial) for monomial in value.denom):
            raise ValueError('Failed to clear a rational denominator')
        return value.numer / value.denom[(0, 0)]

    def poly_jet(self, polynomial):
        degree = max((i for i, j in polynomial), default=0)
        out = [self.K.zero for _ in range(degree+1)]
        for (i, j), value in polynomial.items():
            for k in range(i+1):
                c = int(value)*comb(i, k)*pow(self.t0, i-k, self.p) % self.p
                out[k] += c*self.zz**j
        return out

    def specialize(self, value):
        return self.poly_jet(value.numer)[0]/self.poly_jet(value.denom)[0]

    def flatness_check(self):
        t, z = self.L.gens
        for i in range(self.rank):
            for j in range(self.rank):
                x = self.At[i][j].diff(z)-self.Az[i][j].diff(t)
                x += sum(self.At[i][k]*self.Az[k][j]-self.Az[i][k]*self.At[k][j]
                         for k in range(self.rank))
                if x:
                    raise ValueError(f'FBI connection is not flat at ({i},{j})')
        return True


def generate_series(model, order, progress=None):
    if order < 0 or order+1 >= model.p:
        raise ValueError('Require 0 <= order+1 < prime')
    k = model.K
    r = model.rank
    identity = [[k.one if i == j else k.zero for j in range(r)] for i in range(r)]
    f = [identity]
    object_rows = {label: [] for label in model.targets}
    integrated = {label: [] for label in model.targets}
    certificates = {label: [] for label in model.targets}
    timers = {'fbi_jet_seconds': 0., 'target_jet_seconds': 0., 'integration_seconds': 0.}
    for n in range(order+1):
        before = time.perf_counter()
        if n:
            next_row = [[k.zero for _ in range(r)] for _ in range(r)]
            for i in range(r):
                for j in range(r):
                    v = k.zero
                    for a in range(r):
                        for d, coefficient in enumerate(model.nt[i][a]):
                            if d < n and coefficient:
                                v += coefficient*f[n-1-d][a][j]
                    for d in range(1, min(len(model.qt), n+1)):
                        v -= model.qt[d]*(n-d)*f[n-d][i][j]
                    next_row[i][j] = v/model.qt[0]/n
            f.append(next_row)
        timers['fbi_jet_seconds'] += time.perf_counter()-before
        for label, (den, numerator) in model.object_data.items():
            before = time.perf_counter()
            row = []
            for j in range(r):
                v = k.zero
                for i in range(r):
                    for d, coefficient in enumerate(numerator[i]):
                        if d <= n and coefficient:
                            v += coefficient*f[n-d][i][j]
                for d in range(1, min(len(den), n+1)):
                    v -= den[d]*object_rows[label][n-d][j]
                row.append(v/den[0])
            object_rows[label].append(row)
            timers['target_jet_seconds'] += time.perf_counter()-before
            before = time.perf_counter()
            coordinates, cert = model.reducer.reduce(row)
            integrated[label].append(coordinates)
            certificates[label].append([str(x.as_expr()) for x in cert])
            timers['integration_seconds'] += time.perf_counter()-before
        if progress and (n < 3 or n % 10 == 0 or n == order):
            progress(n)
    keys = sorted({key for series in integrated.values() for row in series for key in row}, key=str)
    return {'coordinates': integrated, 'certificates': certificates,
            'period_keys': keys, 'timings': timers,
            'operations': dict(model.reducer.operations), 'raw_rows': object_rows}
