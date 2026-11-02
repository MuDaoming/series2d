"""Input-driven two-loop, three-branch FBI front end (nonnegative indices).

No topology, kinematics, test mass, or preferred master is compiled here.
The four FBI cases are the same identities as expand/src/series_solver.tpp.
"""
from __future__ import annotations
from functools import lru_cache
from itertools import product
import sympy as s
from sympy.polys.domains import GF
from sympy.polys.matrices import DomainMatrix


class TwoLoopFamily:
    def __init__(self, cfg, log=None):
        self.cfg = cfg
        self.log = log or (lambda text: None)
        self.p = int(cfg['prime'])
        self.t, self.z = s.symbols('t z')
        self.K = GF(self.p).frac_field(self.t, self.z)
        self.P = GF(self.p).poly_ring(self.t, self.z)
        self.kzero = self.K.zero
        self.kone = self.K.one
        self.tg, self.zg = self.K.gens
        self.dimension = self.value(cfg['dimension'])
        self.props = cfg['propagators']
        self.n = len(self.props)
        self.branches = [int(x['branch']) for x in self.props]
        if set(self.branches) != {0, 1, 2}:
            raise ValueError('Exactly three nonempty branches numbered 0,1,2 are required')
        self.cut = tuple(cfg['cut_sector']) if cfg.get('cut_sector') else None
        if self.cut and (len(self.cut) != self.n or any(x not in (0, 1) for x in self.cut)):
            raise ValueError('cut_sector must be a binary propagator support vector')
        self._build_symanzik()
        self.sectors = {}
        for nu in product((0, 1), repeat=self.n):
            if not self.valid(nu) or not self.keep(nu):
                continue
            self.sectors[nu] = self.analyze_sector(nu)
        self.masters = sorted([nu for nu, info in self.sectors.items() if info['case'] == 0],
                              key=lambda x: (-sum(x), tuple(-a for a in x)))
        self.master_index = {nu: i for i, nu in enumerate(self.masters)}
        self.rank = len(self.masters)
        if not self.rank:
            raise ValueError('No retained MFBI in the supplied family/cut')
        self.log(f'FBI sectors={len(self.sectors)}, retained masters={self.rank}')

    def value(self, value):
        expr = s.sympify(value, locals={'t': self.t, 'z': self.z})
        num, den = s.together(expr).as_numer_denom()
        return self.K.field.new(self.K.field.ring.from_expr(num), self.K.field.ring.from_expr(den))

    def _build_symanzik(self):
        X = [self.t*(1-self.z), (1-self.t)*(1-self.z), self.z]
        y = s.symbols(f'y0:{self.n}')
        gram = s.Matrix([[s.Rational(str(a)) for a in row] for row in self.cfg.get('gram', [])])
        e = gram.rows
        if gram.shape != (e, e) or gram != gram.T:
            raise ValueError('External Gram matrix must be square and symmetric')
        branch_loop = {}
        for pr in self.props:
            c = s.Matrix(pr['loop'])
            if c.shape != (2, 1) or len(pr.get('shift', [])) != e:
                raise ValueError('Every propagator needs two loop coefficients and one shift per external momentum')
            outer = c*c.T
            b = pr['branch']
            if b in branch_loop and branch_loop[b] != outer:
                raise ValueError('Propagators in a branch must share their quadratic loop part')
            branch_loop[b] = outer
        A = sum((X[b]*branch_loop[b] for b in range(3)), s.zeros(2))
        U = s.expand(A.det())
        if not U:
            raise ValueError('Degenerate two-loop quadratic form')
        adj = s.Matrix([[A[1, 1], -A[0, 1]], [-A[1, 0], A[0, 0]]])
        B = s.zeros(2, e)
        C = s.Integer(0)
        for i, pr in enumerate(self.props):
            c = s.Matrix(pr['loop'])
            q = s.Matrix([s.Rational(str(a)) for a in pr.get('shift', [])])
            b = pr['branch']
            if e:
                B += X[b]*y[i]*(c*q.T)
                const = (q.T*gram*q)[0]-s.Rational(str(pr['mass2']))
            else:
                const = -s.Rational(str(pr['mass2']))
            C += X[b]*y[i]*const*sum(y[j] for j in range(self.n) if self.branches[j] == b)
        F = s.expand(sum(adj[i, j]*(B[i, :]*gram*B[j, :].T)[0]
                         for i in range(2) for j in range(2)) - U*C) if e else s.expand(-U*C)
        self.U = self.value(U)
        self.R = [[self.value(s.diff(F, y[i], y[j])) for j in range(self.n)]
                  for i in range(self.n)]
        self.symanzik = {'U': str(U), 'F': str(F), 'variables': [str(a) for a in y]}

    def valid(self, nu):
        return len(nu) == self.n and all(x >= 0 for x in nu) and \
            {self.branches[i] for i, x in enumerate(nu) if x > 0} == {0, 1, 2}

    def keep(self, nu):
        return not self.cut or all(not b or nu[i] > 0 for i, b in enumerate(self.cut))

    def analyze_sector(self, nu):
        active = [i for i, x in enumerate(nu) if x]
        n = len(active)+3
        S = [[self.kzero for _ in range(n)] for _ in range(n)]
        for j, prop in enumerate(active):
            b = self.branches[prop]
            S[b][j+3] = S[j+3][b] = self.kone
            for k, other in enumerate(active):
                S[j+3][k+3] = self.R[prop][other]
        polynomial = all(all(all(not any(m) for m in x.denom) for x in row) for row in S)
        if not polynomial:
            raise ValueError('Internal S was expected to be polynomial')
        pdata = [[self.P.ring.from_dict(dict(x.numer/x.denom[(0, 0)])) for x in row] for row in S]
        dm = DomainMatrix(pdata, (n, n), self.P)
        # Eliminate the three affine branch constraints before inversion.
        # This is a general KKT block identity, not a topology-specific shortcut.
        pivots = [next(j for j, prop in enumerate(active) if self.branches[prop] == b)
                  for b in range(3)]
        free = [j for j in range(len(active)) if j not in pivots]
        rr = [[S[3+i][3+j] for j in range(len(active))] for i in range(len(active))]
        v = [[self.value(int(i == j)-int(i == pivots[self.branches[active[j]]]))
              for j in free] for i in range(len(active))]
        h = [[sum((v[i][a]*rr[i][j]*v[j][b]
                   for i in range(len(active)) for j in range(len(active))), self.kzero)
              for b in range(len(free))] for a in range(len(free))]
        if free:
            hp = [[self.P.ring.from_dict(dict(x.numer/x.denom[(0,0)])) for x in row] for row in h]
            hd = DomainMatrix(hp,(len(free),len(free)),self.P)
            determinant = hd.det()
        else:
            determinant = self.P.one
        inv = None
        if determinant:
            if free:
                numerator, denominator = hd.inv_den(method='rref')
                de = self.value(denominator.as_expr())
                hi = [[self.value(x.as_expr())/de for x in row] for row in numerator.to_list()]
            else:
                hi = []
            nn = len(active)
            k = [[sum((v[i][a]*hi[a][b]*v[j][b]
                       for a in range(len(free)) for b in range(len(free))), self.kzero)
                  for j in range(nn)] for i in range(nn)]
            bmat = [[self.value(int(i == pivots[b]))-
                     sum((k[i][j]*rr[j][pivots[b]] for j in range(nn)),self.kzero)
                     for b in range(3)] for i in range(nn)]
            inv = [[self.kzero for _ in range(n)] for _ in range(n)]
            for a in range(3):
                for b in range(3):
                    inv[a][b] = -sum((rr[pivots[a]][j]*bmat[j][b] for j in range(nn)),self.kzero)
                for j in range(nn):
                    inv[a][3+j] = inv[3+j][a] = bmat[j][a]
            for i in range(nn):
                for j in range(nn):
                    inv[3+i][3+j] = k[i][j]
            cz = [sum(row[:3], self.kzero) for row in inv]
            C = sum(cz[:3], self.kzero)
            case = 0 if C else 1
        else:
            null = dm.nullspace().to_list()
            candidates = [[self.value(x.as_expr()) for x in row] for row in null]
            cz = next((v for v in candidates if sum(v[:3], self.kzero)), candidates[0])
            C = sum(cz[:3], self.kzero)
            case = 2 if C else 3
        self.log(f'sector {"".join(map(str,nu))}: case={case}')
        return {'case': case, 'active': active, 'inverse': inv, 'C': C, 'z': cz[3:]}

    def scaled_sum(self, terms):
        out = [self.kzero]*self.rank
        for coefficient, row in terms:
            if coefficient:
                out = [a+coefficient*b for a, b in zip(out, row)]
        return tuple(out)

    @lru_cache(None)
    def reduce(self, nu, shift):
        if not self.valid(nu) or not self.keep(nu):
            return (self.kzero,)*self.rank
        sec = tuple(int(x > 0) for x in nu)
        info = self.sectors[sec]
        active, C, zs = info['active'], info['C'], info['z']
        delta = self.dimension+shift
        total = sum(nu)
        lower = lambda i: tuple(x-int(j == i) for j, x in enumerate(nu))
        if info['case'] == 0:
            if max(nu) <= 1:
                if shift == 0:
                    out = [self.kzero]*self.rank
                    out[self.master_index[nu]] = self.kone
                    return tuple(out)
                if shift > 0:
                    fac = 2*delta-total-3
                    return self.scaled_sum([(C/fac, self.reduce(nu, shift-1))] +
                        [(-zs[j]/fac, self.reduce(lower(i), shift-1)) for j, i in enumerate(active)])
                fac = 2*(delta+1)-total-3
                return self.scaled_sum([(fac/C, self.reduce(nu, shift+1))] +
                    [(zs[j]/C, self.reduce(lower(i), shift)) for j, i in enumerate(active)])
            i = max(active, key=lambda j: nu[j])
            prev = lower(i)
            row = info['inverse'][3+active.index(i)]
            terms = [(-sum(row[:3], self.kzero)/(nu[i]-1), self.reduce(prev, shift-1))]
            for j, prop in enumerate(active):
                src = tuple(x-int(a == prop) for a, x in enumerate(prev))
                terms.append((row[3+j]/(nu[i]-1), self.reduce(src, shift-1)))
            return self.scaled_sum(terms)
        if info['case'] == 1:
            return self.scaled_sum([(-zs[j]/(2*delta-total-3), self.reduce(lower(i), shift-1))
                                    for j, i in enumerate(active)])
        if info['case'] == 2:
            return self.scaled_sum([(zs[j]/C, self.reduce(lower(i), shift)) for j, i in enumerate(active)])
        j = next(j for j, x in enumerate(zs) if x)
        pivot = active[j]
        terms = []
        for a, i in enumerate(active):
            if a == j:
                continue
            src = tuple(x+int(k == pivot)-int(k == i) for k, x in enumerate(nu))
            terms.append((-zs[a]/zs[j], self.reduce(src, shift)))
        return self.scaled_sum(terms)

    def derivative(self, nu, variable):
        terms = []
        active = [i for i, x in enumerate(nu) if x > 0]
        for i in active:
            for j in active:
                coeff = self.R[i][j].diff(variable)
                if not coeff:
                    continue
                factor = nu[i]*(nu[j]+int(i == j))
                src = tuple(x+int(a == i)+int(a == j) for a, x in enumerate(nu))
                terms.append((-coeff*factor/2, self.reduce(src, 1)))
        return self.scaled_sum(terms)

    def integer_prefactor(self, nu):
        X = [self.tg*(1-self.zg), (1-self.tg)*(1-self.zg), self.zg]
        W = self.kone
        for b in range(3):
            W *= X[b]**(sum(nu[i] for i in range(self.n) if self.branches[i] == b)-1)
        return (1-self.zg)*W*self.U**sum(nu)

    def connection_model(self):
        pref = [self.integer_prefactor(nu) for nu in self.masters]
        matrices = []
        for variable in (self.tg, self.zg):
            matrix = []
            for i, nu in enumerate(self.masters):
                raw = self.derivative(nu, variable)
                row = [pref[i]*raw[j]/pref[j] for j in range(self.rank)]
                row[i] += pref[i].diff(variable)/pref[i] - \
                    3*self.dimension*self.U.diff(variable)/(2*self.U)
                matrix.append([str(x.as_expr()) for x in row])
            matrices.append(matrix)
        targets = {}
        requested = self.cfg['targets']
        for values in requested:
            nu = tuple(values)
            if not self.valid(nu):
                raise ValueError(f'Target is outside the nonnegative three-branch scope: {nu}')
            raw = self.reduce(nu, 0)
            fac = self.integer_prefactor(nu)
            targets['J{'+','.join(map(str,nu))+'}'] = [str((fac*raw[j]/pref[j]).as_expr())
                                                    for j in range(self.rank)]
        return {'prime': self.p, 't0': self.cfg['t0'], 'At': matrices[0], 'Az': matrices[1],
                'targets': targets, 'endpoints': self.cfg['endpoints'],
                'metadata': {'mode': 'cut' if self.cut else 'full', 'cut_sector': self.cut,
                             'fbi_masters': self.masters, 'symanzik': self.symanzik,
                             'dimension': self.cfg['dimension']}}
