"""K^top = FreeLie(G_T) / ( [G_{S+i},G_j] - [G_{S+j},G_i] ), the excess-leading-form algebra.
Compare multilinear dims with K (kmodel.py).  If equal, the excess filtration is flat and
l -> K is injective in that multidegree."""
import sys, itertools
import kmodel
from kmodel import *

def relations_top(S, maxS0=None):
    idx = [1<<t for t in range(S.bit_length()) if S>>t & 1]
    res = []
    for i, j in itertools.combinations(idx, 2):
        S0 = S & ~i & ~j
        v = {}
        for w, c in br2(S0|i, j).items(): v[w] = (v.get(w,0)+c) % P
        for w, c in br2(S0|j, i).items(): v[w] = (v.get(w,0)-c) % P
        v = {w: c for w, c in v.items() if c}
        if v: res.append(v)
    return res
kmodel.relations = relations_top
N = int(sys.argv[1])
from math import factorial
I = kmodel.compute(N)
for n in range(2, N+1):
    S = (1<<n)-1
    print(n, [stirling2(n,k)*factorial(k-1) - len(I[S].get(k, [])) for k in range(1, n+1)], flush=True)
