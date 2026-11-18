"""Cross-check of the identification with prior work (Losev-Manin 2000, Dotsenko-Shadrin-Tamaroff 2022).

(1) Build the multilinear part of U(k) = k<G_T>/(r_{i,j,S}) directly from the relators r_{i,j,S} of the paper
    (eq. (r)), restricted to content [n]; compute its dimension by D-degree (= n - number of letters).
    The same space is, by definition, arity n of the twisted associative algebra tHyperCom of
    [DST, Prop. 2.4.2] (generators m_T <-> G_T, relators R_{i,j} <-> r_{i,j,S} with S = [n]-{i,j}).
    We check: the DST relator R_{i,j} literally coincides with the paper's r_{i,j,S}; dimension = n!;
    graded dimensions = Eulerian numbers = Poincare polynomial p_n(q) of L_n [LM, Thm 2.3, (2.6)].
(2) sympy: EGF(U(g)) * e^{-xy} at z = 1 equals LM's (q-1)/(q - e^{(q-1)x}) with q = y.
Exact arithmetic over Q (fractions)."""
from fractions import Fraction
from itertools import combinations
import sympy as sp

def ordered_set_partitions(s):
    s = tuple(s)
    if not s:
        yield ()
        return
    n = len(s)
    for r in range(1, n + 1):
        for first in combinations(s, r):
            rest = tuple(x for x in s if x not in first)
            for p in ordered_set_partitions(rest):
                yield (frozenset(first),) + p

def relator_paper(i, j, S):
    """r_{i,j,S} = sum_{P subset S} [G_{P+i}, G_{S-P+j}] as dict word->coeff (associative expansion)."""
    S = tuple(S); out = {}
    for r in range(len(S) + 1):
        for P in combinations(S, r):
            A = frozenset(P) | {i}; B = frozenset(S) - frozenset(P) | {j}
            out[(A, B)] = out.get((A, B), 0) + 1
            out[(B, A)] = out.get((B, A), 0) - 1
    return {k: v for k, v in out.items() if v}

def relator_DST(i, j, n):
    """R_{i,j} = sum_{n=I u J, i in I, j in J} m_I m_J - m_J m_I."""
    N = frozenset(range(1, n + 1)); out = {}
    others = sorted(N - {i, j})
    for r in range(len(others) + 1):
        for P in combinations(others, r):
            I = frozenset(P) | {i}; J = N - I
            out[(I, J)] = out.get((I, J), 0) + 1
            out[(J, I)] = out.get((J, I), 0) - 1
    return {k: v for k, v in out.items() if v}

def rank_Q(rows, ncols):
    rows = [dict(r) for r in rows if r]
    piv = {}; rank = 0
    for r in rows:
        r = {k: Fraction(v) for k, v in r.items() if v}
        while r:
            c = min(r)
            if c in piv:
                pr = piv[c]; f = r[c] / pr[c]
                for k, v in pr.items():
                    nv = r.get(k, 0) - f * v
                    if nv: r[k] = nv
                    else: r.pop(k, None)
            else:
                piv[c] = r; rank += 1; break
    return rank

def multilinear_dims(n):
    N = list(range(1, n + 1))
    words = list(ordered_set_partitions(N))
    idx = {w: t for t, w in enumerate(words)}
    rows = []
    # ideal component of content [n]: u * r * v with u, v words, contents disjoint, union [n]
    # enumerate: choose ordered set partition of [n] into (U-blocks..., X, V-blocks...) where X is the support of r
    for p in ordered_set_partitions(N):
        for pos, X in enumerate(p):
            if len(X) < 2: continue
            u = p[:pos]; v = p[pos + 1:]
            xs = sorted(X)
            i = xs[0]
            for j in xs[1:]:   # r_{i,j}, j ranges: spans all r_{a,b} on X since R_ab = R_ai + R_ib
                S = [s for s in xs if s not in (i, j)]
                rel = relator_paper(i, j, S)
                row = {}
                for (A, B), c in rel.items():
                    row[idx[u + (A, B) + v]] = c
                rows.append(row)
    # graded by number of letters p
    dims = {}
    bydeg = {}
    for w in words: bydeg.setdefault(len(w), []).append(idx[w])
    for k, cols in bydeg.items():
        cs = set(cols)
        rws = [r for r in rows if r and next(iter(r)) in cs]   # rows are homogeneous in the number of letters
        dims[k] = len(cols) - rank_Q(rws, len(cols))
    return dims

if __name__ == "__main__":
    # (0) literal coincidence of relators
    for n in range(2, 7):
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                if i == j: continue
                S = [s for s in range(1, n + 1) if s not in (i, j)]
                assert relator_paper(i, j, S) == relator_DST(i, j, n)
    print("relators r_{i,j,S} (paper) == R_{i,j} (DST Prop 2.4.2) for all n<=6, i!=j: OK")
    q, x, y, z = sp.symbols('q x y z')
    for n in range(1, 7):
        d = multilinear_dims(n)
        # D-degree a = n - (#letters)
        print(f"n={n}: dim U(k)_[n] by D-degree a=0..{n-1}:", [d[n - a] for a in range(n)], " total", sum(d.values()), " n! =", sp.factorial(n))
    # LM (2.6)
    LM = (q - 1) / (q - sp.exp((q - 1) * x))
    ser = sp.series(LM, x, 0, 8).removeO()
    for n in range(1, 8):
        c = sp.expand(sp.simplify(ser.coeff(x, n) * sp.factorial(n)))
        print(f"LM p_{n}(q) =", sp.Poly(c, q).all_coeffs()[::-1])
    # (2) EGF of U(g) from the paper: 1/H(-x) with (y-z)H(x) = y e^{xy} - z e^{xz}; at z=1
    H = lambda t: (y * sp.exp(t * y) - z * sp.exp(t * z)) / (y - z)
    Ug = 1 / H(-x)
    expr = sp.simplify((Ug * sp.exp(-x * y)).subs(z, 1).subs(y, q) - LM)
    print("EGF(U(g)) * e^{-xy} |_{z=1,y=q} - (q-1)/(q-e^{(q-1)x}) simplifies to:", expr)
    diff = sp.series((Ug * sp.exp(-x * y)).subs(z, 1).subs(y, q) - LM, x, 0, 10).removeO()
    print("series difference through x^9:", sp.simplify(diff))
