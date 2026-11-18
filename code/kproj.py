"""Direct injectivity check in multilinear degree n (mod p), memory-lean.
pi: F -> F_big = coordinates (words) containing at least one block of size >= 3.
dim(K/h) = dim F/(Fs + I) = dim F_big - rank pi(I).
Recursion: pi(I(S)) = pi(R(S)) + sum_{y small} [y, pi(I(S-y))] + sum_{y big} [y, I(S-y)]
(for y small, pi[y,u] = [y,pi u]; for y big, [y,u] is already in F_big).
Full I(S') is needed only for |S'| <= n-3.
Injective in degree (n,k)  <=>  dim(K/h)(n,k) == dim K(n,k) - dim l(n,k)."""
import sys, time, resource
import kmodel as km
from kmodel import Echelon, popcount, relations, blocks_of, bracket_letter, I_get, is_small, stirling2, free_small_dim
from math import factorial, comb

def big(w): return not is_small(w)

def run(n, kmin):
    full = {}   # canonical size m -> {k: rows} full ideal, m <= n-3
    proj = {}   # canonical size m -> {k: rows} projected ideal
    for m in range(2, n+1):
        S = (1 << m) - 1
        eF = {}; eP = {}
        ks = range(1, m+1) if m < n else range(kmin, m+1)
        for k in ks:
            if m <= n-3:
                e = Echelon()
                for r in relations(S):
                    if len(next(iter(r))) == k: e.add(r)
                for y in blocks_of(S):
                    for row in I_get(full, S & ~y).get(k-1, []):
                        e.add(bracket_letter(y, row, S & ~y))
                eF[k] = e.rows()
            e = Echelon()
            for r in relations(S):
                if len(next(iter(r))) == k:
                    v = {w: c for w, c in r.items() if big(w)}
                    if v: e.add(v)
            for y in blocks_of(S):
                Sp = S & ~y
                src = I_get(proj, Sp) if popcount(y) <= 2 else I_get(full, Sp)
                for row in src.get(k-1, []):
                    v = bracket_letter(y, row, Sp)
                    if popcount(y) <= 2: v = {w: c for w, c in v.items() if big(w)}  # (already big)
                    if v: e.add(v)
            eP[k] = e.rows()
            if m == n:
                fb = stirling2(n, k)*factorial(k-1) - free_small_dim(n, k)
                q = fb - len(eP[k])
                print("n=%d k=%d  dim F_big=%d  rank pi(I)=%d  dim K/h=%d  maxRSS=%.0fMB" %
                      (n, k, fb, len(eP[k]), q, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024), flush=True)
                eP[k] = []  # free
        if m <= n-3: full[S] = eF
        # store canonical
        proj[S] = eP
        # relabel helper expects dict keyed by canonical S
    return

def I_get_wrap(store, S):
    return I_get(store, S)
if __name__ == '__main__':
    n = int(sys.argv[1]); kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    t = time.time(); run(n, kmin); print("time %.0fs" % (time.time()-t))
