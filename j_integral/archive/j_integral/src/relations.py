"""Relations among exact fixed-period coordinate jets, with held-out orders."""
from __future__ import annotations
import sympy as s
from sympy.polys.domains import GF


def rref(a, p):
    if not a:
        return [], []
    a = [[x % p for x in row] for row in a]
    pivots = []
    for col in range(len(a[0])):
        row = len(pivots)
        if row == len(a):
            break
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        inv = pow(a[row][col], -1, p)
        a[row] = [x*inv % p for x in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                c = a[i][col]
                a[i] = [(x-c*y) % p for x, y in zip(a[i], a[row])]
        pivots.append(col)
    return a[:len(pivots)], pivots


def search_relations(series, p, degree, train_order):
    labels = list(series)
    order = min(len(v) for v in series.values())-1
    if not 0 <= train_order < order:
        raise ValueError('Require a nonempty held-out coefficient window')
    keys = sorted({key for v in series.values() for row in v for key in row}, key=str)
    variables = [(label, d) for label in labels for d in range(degree+1)]
    rows = [[series[label][n-d].get(key, 0) if n >= d else 0
             for label, d in variables] for key in keys for n in range(train_order+1)]
    rr, pivots = rref(rows, p)
    free = [c for c in range(len(variables)) if c not in pivots]
    accepted, rejected = [], 0
    for fc in free:
        v = [0]*len(variables)
        v[fc] = 1
        for row, pc in enumerate(pivots):
            v[pc] = -rr[row][fc] % p
        ok = True
        for key in keys:
            for n in range(train_order+1, order+1):
                residue = sum(c*series[label][n-d].get(key, 0)
                              for c, (label, d) in zip(v, variables) if n >= d)
                if residue % p:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            accepted.append({label: v[i*(degree+1):(i+1)*(degree+1)]
                             for i, label in enumerate(labels)})
        else:
            rejected += 1
    return {'degree': degree, 'train_order': train_order, 'verification_order': order,
            'projection_count': len(keys), 'train_rank': len(pivots),
            'candidate_count': len(free), 'rejected_count': rejected,
            'relations': accepted}


def solve_rational_relations(labels, relations, p):
    tau = s.Symbol('tau')
    K = GF(p).frac_field(tau)
    x = K.gens[0]
    rows = [[sum((c*x**i for i, c in enumerate(rel[label])), K.zero)
             for label in labels] for rel in relations]
    pivot_columns = []
    for col in range(len(labels)):
        row = len(pivot_columns)
        pivot = next((i for i in range(row, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[row], rows[pivot] = rows[pivot], rows[row]
        lead = rows[row][col]
        rows[row] = [v/lead for v in rows[row]]
        for i in range(len(rows)):
            if i != row and rows[i][col]:
                c = rows[i][col]
                rows[i] = [a-c*b for a, b in zip(rows[i], rows[row])]
        pivot_columns.append(col)
    free = [c for c in range(len(labels)) if c not in pivot_columns]
    reductions = {labels[pc]: {labels[fc]: str((-rows[i][fc]).as_expr())
                               for fc in free if rows[i][fc]}
                  for i, pc in enumerate(pivot_columns)}
    return {'masters': [labels[fc] for fc in free], 'reductions': reductions,
            'master_count': len(free)}


def search_master_de(series, masters, p, degree, train_order):
    output = {}
    for label in masters:
        derivative = [{key: (n+1)*v % p for key, v in series[label][n+1].items()
                       if (n+1)*v % p} for n in range(len(series[label])-1)]
        data = {'DERIVATIVE': derivative}
        data.update({m: series[m][:-1] for m in masters})
        found = None
        for d in range(degree+1):
            result = search_relations(data, p, d, min(train_order, len(derivative)-2))
            for rel in result['relations']:
                if any(rel['DERIVATIVE']):
                    found = {'degree': d, 'relation': rel,
                             'verification_order': result['verification_order']}
                    break
            if found:
                break
        output[label] = found or {'status': 'unresolved; supplied J set may not be derivative-closed'}
    return output
