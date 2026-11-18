import itertools, sys, ast
from fractions import Fraction
from math import factorial
from collections import Counter
def partitions(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for p in partitions(n-k, k): yield (k,)+p
def mn_char(lam, mu):
    # Murnaghan-Nakayama via beta-sets
    if not mu: return 1 if sum(lam) == 0 else 0
    r = mu[0]; rest = mu[1:]
    L = len(lam)
    beta = [lam[i] + (L-1-i) for i in range(L)]
    tot = 0
    for i, b in enumerate(beta):
        nb = b - r
        if nb < 0 or nb in beta: continue
        newbeta = beta[:i] + [nb] + beta[i+1:]
        sign = (-1)**sum(1 for x in beta if nb < x < b)
        newbeta.sort(reverse=True)
        Ln = len(newbeta)
        newlam = tuple(newbeta[j] - (Ln-1-j) for j in range(Ln))
        newlam = tuple(x for x in newlam if x > 0)
        tot += sign * mn_char(newlam, rest)
    return tot
def classsize(n, mu):
    c = Counter(mu); d = 1
    for k, m in c.items(): d *= k**m * factorial(m)
    return factorial(n)//d
def decompose(n, chars):
    res = {}
    for lam in partitions(n):
        s = sum(classsize(n, mu)*chars[mu]*mn_char(lam, mu) for mu in chars)
        s = Fraction(s, factorial(n))
        if s: res[lam] = s
    return res
for line in sys.stdin:
    if not line.startswith('k'): continue
    parts = line.split(' ', 4)
    k = int(parts[1]); data = ast.literal_eval(parts[4])
    n = sum(data[0][0])
    chars = {tuple(ct): v for ct, v in data}
    print("n", n, "k", k, {''.join(map(str,l)): int(m) for l, m in decompose(n, chars).items()})
