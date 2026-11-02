"""Exact coefficient integration at a fixed ordinary t point.

The input connection belongs to the *FBI integrand*, not to J.  All objects
are over F_p(z).  A reduction records r = remainder + g' + g A_z, so that
integrating requires only fixed-point periods and explicit endpoint terms.
"""
from __future__ import annotations

from collections import defaultdict
import sympy as s
from sympy.polys.domains import GF


def decompose(matrix, p):
    """Return a deterministic partial inverse. Residuals are never discarded."""
    n = len(matrix)
    a = [[int(x) % p for x in row] + [int(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    pivots = []
    for col in range(n):
        row = len(pivots)
        pivot = next((i for i in range(row, n) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        inv = pow(a[row][col], -1, p)
        a[row] = [x * inv % p for x in a[row]]
        for i in range(n):
            if i != row and a[i][col]:
                c = a[i][col]
                a[i] = [(x - c*y) % p for x, y in zip(a[i], a[row])]
        pivots.append(col)
    inverse = [[0] * n for _ in range(n)]
    for i, col in enumerate(pivots):
        inverse[col] = a[i][n:]
    return inverse, len(pivots)


class FixedPeriodReducer:
    def __init__(self, p, connection, extra_denominators=(), endpoints=None):
        self.p = p
        self.z = s.Symbol('z')
        self.K = GF(p).frac_field(self.z)
        self.R = self.K.field.ring
        self.zz = self.K.gens[0]
        self.A = [[self.coerce(x) for x in row] for row in connection]
        self.r = len(self.A)
        if not self.r or any(len(row) != self.r for row in self.A):
            raise ValueError('FBI connection must be nonempty and square')
        self.zero = self.K.zero
        self.one = self.K.one
        q = self.R.one
        for x in [x for row in self.A for x in row] + list(extra_denominators):
            q = q.lcm(self.coerce(x).denom)
        self.Q = q.exquo(q.gcd(q.diff(self.R.gens[0]))).monic()
        self.q = self.Q.degree()
        if self.q < 1:
            raise ValueError('Polynomial-only connections need a separate infinity treatment')
        self.Qrat = self.K.field.new(self.Q)
        self.N = [[self.polynomial(self.Qrat*x) for x in row] for row in self.A]
        if any(x.degree() > self.q-1 for row in self.N for x in row):
            raise ValueError('Connection has a polynomial part at infinity; not supported yet')
        self.endpoints = endpoints or {0: 'regular', 1: 'regular'}
        if any(mode not in ('regular', 'regulated_zero') for mode in self.endpoints.values()):
            raise ValueError('Endpoint mode must be regular or regulated_zero')
        for point, mode in self.endpoints.items():
            if mode == 'regular':
                for row in self.A:
                    for x in row:
                        self.evaluate(x, point)
        self._maps = {}
        self._infinity = {}
        self.operations = {'pole_steps': 0, 'infinity_steps': 0,
                           'coefficient_integrations': 0, 'certificate_checks': 0}

    def coerce(self, value):
        if isinstance(value, str):
            value = s.sympify(value, locals={'z': self.z})
        if isinstance(value, s.Basic):
            num, den = s.together(value).as_numer_denom()
            return self.K.field.new(self.R.from_expr(num), self.R.from_expr(den))
        return self.K.convert(value)

    def polynomial(self, x):
        x = self.K.convert(x)
        if x.denom.degree() != 0:
            raise ValueError('Expected a polynomial after clearing denominators')
        return x.numer / x.denom[(0,)]

    def evaluate(self, x, point):
        x = self.K.convert(x)
        num = x.numer.evaluate(self.R.gens[0], point)
        den = x.denom.evaluate(self.R.gens[0], point)
        if not den:
            raise ValueError(f'Nonregular endpoint {point}; need explicit endpoint analysis')
        return int(num) * pow(int(den), -1, self.p) % self.p

    def _pole_map(self, k):
        if k not in self._maps:
            size = self.r * self.q
            matrix = [[0]*size for _ in range(size)]
            z = self.R.gens[0]
            dq = self.Q.diff(z)
            for source in range(self.r):
                for j in range(self.q):
                    a = z**j
                    for target in range(self.r):
                        poly = a*self.N[source][target]
                        if source == target:
                            poly -= (k-1)*dq*a
                        poly = poly.rem(self.Q)
                        for (degree,), coefficient in poly.items():
                            matrix[target*self.q+degree][source*self.q+j] = int(coefficient)
            inverse, rank = decompose(matrix, self.p)
            self._maps[k] = inverse, rank
        return self._maps[k]

    def _infinity_map(self, power):
        if power not in self._infinity:
            lcq = self.Q[(self.q,)]
            matrix = [[int(self.N[j][i].get((self.q-1,), 0) +
                           (power*lcq if i == j else 0)) % self.p
                       for j in range(self.r)] for i in range(self.r)]
            self._infinity[power] = decompose(matrix, self.p)
        return self._infinity[power]

    def _solve(self, inverse, rhs):
        return [sum(x*y for x, y in zip(row, rhs)) % self.p for row in inverse]

    def covariant(self, g):
        return [g[j].diff(self.zz) + sum(g[i]*self.A[i][j] for i in range(self.r))
                for j in range(self.r)]

    def reduce(self, row, verify=True):
        """Integrate a rational row modulo exact z derivatives.

        Keys ('period', k, component, power) mean integral z^power/Q^k*f_component.
        Keys ('endpoint', point, component) mean f_component(point).
        Obstructions of singular reduction maps are retained as extra periods.
        """
        row = [self.coerce(x) for x in row]
        power = 1
        qp = self.Q
        while any(qp.rem(x.denom) for x in row):
            power += 1
            if power > 10000:
                raise ValueError('Coefficient has a denominator outside the fixed singular locus')
            qp *= self.Q
        polys = [self.polynomial(x*self.K.field.new(qp)) for x in row]
        cert = [self.zero]*self.r
        classes = defaultdict(int)
        z = self.R.gens[0]
        while power > 1:
            rhs = [int(polys[i].rem(self.Q).get((j,), 0))
                   for i in range(self.r) for j in range(self.q)]
            inverse, _ = self._pole_map(power)
            av = self._solve(inverse, rhs)
            a = [sum((z**j*self.R.domain.convert(av[i*self.q+j])
                      for j in range(self.q)), self.R.zero) for i in range(self.r)]
            numer = [self.Q*a[j].diff(z) - (power-1)*self.Q.diff(z)*a[j] +
                     sum((a[i]*self.N[i][j] for i in range(self.r)), self.R.zero)
                     for j in range(self.r)]
            for j in range(self.r):
                left = polys[j] - numer[j]
                obstruction = left.rem(self.Q)
                for (degree,), value in obstruction.items():
                    classes[('period', power, j, degree)] += int(value)
                polys[j] = (left-obstruction).exquo(self.Q)
                cert[j] += self.K.field.new(a[j], self.Q**(power-1))
            power -= 1
            self.operations['pole_steps'] += 1
        while any(poly.degree() >= self.q-1 for poly in polys):
            degree = max(poly.degree() for poly in polys)
            n = degree-self.q+1
            rhs = [int(poly.get((degree,), 0)) for poly in polys]
            inverse, _ = self._infinity_map(n)
            av = self._solve(inverse, rhs)
            a = [z**n*self.R.domain.convert(c) for c in av]
            numer = [self.Q*a[j].diff(z) +
                     sum((a[i]*self.N[i][j] for i in range(self.r)), self.R.zero)
                     for j in range(self.r)]
            for j in range(self.r):
                polys[j] -= numer[j]
                leftover = int(polys[j].get((degree,), 0))
                if leftover:
                    classes[('period', 1, j, degree)] += leftover
                    polys[j] -= z**degree*self.R.domain.convert(leftover)
                cert[j] += self.K.field.new(a[j])
            self.operations['infinity_steps'] += 1
        for j, poly in enumerate(polys):
            for (degree,), value in poly.items():
                classes[('period', 1, j, degree)] += int(value)
        classes = {key: value % self.p for key, value in classes.items() if value % self.p}
        if verify:
            remainder = [self.zero]*self.r
            for (_, k, component, degree), value in classes.items():
                remainder[component] += value*self.zz**degree/self.Qrat**k
            derivative = self.covariant(cert)
            if any(row[i] - remainder[i] - derivative[i] for i in range(self.r)):
                raise ArithmeticError('Exact coefficient-integration certificate failed')
            self.operations['certificate_checks'] += 1
        for point, mode in self.endpoints.items():
            if mode == 'regulated_zero':
                continue
            sign = -1 if point == 0 else 1
            for j, g in enumerate(cert):
                c = sign*self.evaluate(g, point) % self.p
                if c:
                    classes[('endpoint', point, j)] = c
        self.operations['coefficient_integrations'] += 1
        return classes, cert
