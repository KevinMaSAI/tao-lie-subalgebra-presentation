# S_n character of multilinear K(n,k) (method A coordinates)
import sys, itertools
from kmodel import *
from math import factorial
def full_reduce(e, v):
    v = dict(v); out = {}
    while v:
        c = max(v); r = e.piv.get(c)
        if r is None:
            out[c] = v.pop(c); continue
        f = v[c]
        for cc, x in r.items():
            y = (v.get(cc, 0) - f*x) % P
            if y: v[cc] = y
            else: v.pop(cc, None)
    return out

def to_pivot(word_coeffs):
    """given Lie element as dict over arbitrary-ordered words that is (relabeled) pivot coords
    of an element w.r.t. old pivot, convert to new pivot coords"""
    pass

def relabel_elem(Pv, perm_bits):
    # Pv in pivot coords; apply index permutation; result: Lie elem whose coords are on words starting with image of old pivot
    # convert: new coeff on word starting with new pivot q: coef_full of old element on preimage word
    inv = {v: k for k, v in perm_bits.items()}
    def mapb(b, f):
        r = 0
        for t, u in f.items():
            if b & t: r |= u
        return r
    out = {}
    parts = set(frozenset(w) for w in Pv)
    for part in parts:
        newletters = [mapb(b, perm_bits) for b in part]
        allm = 0
        for b in newletters: allm |= b
        q = [b for b in newletters if b & lowbit(allm)][0]
        rest = [b for b in newletters if b != q]
        for perm in itertools.permutations(rest):
            nw = (q,)+perm
            ow = tuple(mapb(b, inv) for b in nw)
            c = coef_full(Pv, ow)
            if c: out[nw] = c
    return out

N = int(sys.argv[1])
IK = compute(N)
S = (1<<N)-1
def cycle_types(n):
    def parts(n, m):
        if n == 0: yield (); return
        for k in range(min(n, m), 0, -1):
            for p in parts(n-k, k): yield (k,)+p
    return list(parts(n, n))
for k in range(1, N+1):
    rows = IK[S].get(k, [])
    e = Echelon()
    for r in rows: e.add(r)
    # free coords for k
    coords = []
    def setparts(elems, k):
        if not elems:
            if k == 0: yield []
            return
        first, rest = elems[0], elems[1:]
        for p in setparts(rest, k-1): yield [1<<first] + p
        for p in setparts(rest, k):
            for i in range(len(p)):
                yield p[:i] + [p[i] | (1<<first)] + p[i+1:]
    for part in setparts(list(range(N)), k):
        piv = [b for b in part if b & 1][0]
        rest = [b for b in part if b != piv]
        for perm in itertools.permutations(rest): coords.append((piv,)+perm)
    comp = [c for c in coords if c not in e.piv]
    chars = []
    for ct in cycle_types(N):
        # permutation with cycle type ct on bits 0..N-1
        perm = {}; pos = 0
        for L in ct:
            for t in range(L): perm[1<<(pos+t)] = 1<<(pos+(t+1) % L)
            pos += L
        tr = 0
        for c in comp:
            v = relabel_elem({c: 1}, perm)
            w = full_reduce(e, v)
            tr += w.get(c, 0)
        tr %= P
        if tr > P//2: tr -= P
        chars.append((ct, tr))
    print("k", k, "dim", len(comp), chars, flush=True)
