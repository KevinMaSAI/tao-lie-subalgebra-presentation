"""Rigorous dims of multilinear components of U(g) via its PBW (Groebner) basis of normal words
(diamond lemma; PBW property verified exactly in pbw_check.py), refined by D-degree, n <= NMAX.
Then dims of g via the PBW theorem U(g) = Sym(g) (EGF: log), compared with the formula
x(y+z) - log(1-Y) (= K^top ⋊ D, Lemma B) and hence with the Koszul prediction."""
import sys
from collections import defaultdict
import sympy as sp
import pbw_check as pc
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
lk = lambda l: (l[1], 'DE'.index(l[0]))          # letter order: index first, then D<E (a PBW order)
key = lambda w: tuple(lk(l) for l in w)
_, _, leads2, _ = pc.check(pc.Ug_relations(2), 2, key)
_, _, leads3, _ = pc.check(pc.Ug_relations(3), 3, key)
pat = set((a[0], b[0], a[1] > b[1]) for (a, b) in leads2 if a[1] != b[1])
pat3 = set((a[0], b[0], a[1] > b[1]) for (a, b) in leads3 if a[1] != b[1])
assert pat == pat3
print("leading words (multilinear patterns typeA,typeB,idxA>idxB):", sorted(pat))
def count(n):
    dp = defaultdict(lambda: [0]*(n+1))
    for i in range(n):
        for t in 'DE':
            dp[(1 << i, i, t)][1 if t == 'D' else 0] += 1
    for _ in range(n-1):
        nd = defaultdict(lambda: [0]*(n+1))
        for (S, i, t), c in dp.items():
            for j in range(n):
                if S >> j & 1: continue
                for u in 'DE':
                    if (t, u, i > j) in pat: continue
                    tgt = nd[(S | 1 << j, j, u)]
                    sh = 1 if u == 'D' else 0
                    for a in range(n+1-sh):
                        if c[a]: tgt[a+sh] += c[a]
        dp = nd
    tot = [0]*(n+1)
    for c in dp.values():
        for a in range(n+1): tot[a] += c[a]
    return tot
x, y, z = sp.symbols('x y z')
U = 1
for n in range(1, NMAX+1):
    c = count(n)
    print("U(g) multilinear n=%d by D-degree:" % n, c, "total", sum(c), flush=True)
    U += sum(c[a]*y**a*z**(n-a) for a in range(n+1))*x**n/sp.factorial(n)
G = sp.expand(sp.series(sp.log(U), x, 0, NMAX+1).removeO())
Y = sum(x**w/sp.factorial(w)*sum(y**(s-1)*z**(w-s+1) for s in range(2, w+1)) for w in range(2, NMAX+1))
Gpred = sp.expand(x*(y+z) + sp.series(-sp.log(1-Y), x, 0, NMAX+1).removeO())
print("log(PBW count EGF) - (x(y+z) - log(1-Y)) =", sp.expand(G - Gpred))
for n in range(1, NMAX+1):
    cc = sp.Poly(sp.expand(G.coeff(x, n)*sp.factorial(n)), y, z)
    print("g multilinear n=%d by D-degree:" % n, [cc.coeff_monomial(y**a*z**(n-a)) for a in range(n+1)])
