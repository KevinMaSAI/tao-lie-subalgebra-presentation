"""
Method A: the ideal K = <E> of g, presented as
  K = FreeLie(G_T : T nonempty set of indices) / ( r_{i,j,S} ),
  r_{i,j,S} = sum_{A u B = S} [G_{A+i}, G_{B+j}]   (i != j, i,j notin S)
(multilinear part only: every index appears at most once).
G_T = ad(D_{T-t}) E_t (symmetric in T).  h = subalgebra gen. by G_T, |T|<=2.
l = FreeLie(G_T, |T|<=2) / (r with |S|<=1).
Multilinear free Lie elements on distinct letters are stored by their
coefficients on associative words that start with the pivot letter
(= block containing min index); left-normed [[p,a1],...] basis.
Arithmetic mod prime P.
"""
# EXACT RATIONAL VERSION (auto-generated from kmodel.py)
import sys, itertools, time
from functools import lru_cache
from fractions import Fraction
class _Q:
    def __rmod__(self, x): return x
P = _Q()

def popcount(x): return bin(x).count('1')
def lowbit(x): return x & -x

def subsets(m):
    s = m
    while True:
        yield s
        if s == 0: break
        s = (s-1) & m

class Echelon:
    """row echelon basis mod P; rows are dicts col->val"""
    def __init__(self, key=None):
        self.piv = {}   # pivot col -> row (dict), normalized pivot=1
        self.key = key  # column ordering key (pivot = max key)
    def reduce(self, v):
        v = dict(v)
        key = self.key
        while v:
            c = max(v, key=key) if key else max(v)
            r = self.piv.get(c)
            if r is None: return v, c
            f = v[c]
            for cc, x in r.items():
                y = (v.get(cc, 0) - f*x) % P
                if y: v[cc] = y
                else: v.pop(cc, None)
        return v, None
    def add(self, v):
        v = {c: x % P for c, x in v.items() if x % P}
        w, c = self.reduce(v)
        if not w: return False
        inv = Fraction(1, 1) / w[c]
        self.piv[c] = {cc: x*inv % P for cc, x in w.items()}
        return True
    def rank(self): return len(self.piv)
    def rows(self): return list(self.piv.values())

def coef_full(Pv, word):
    """coefficient of P (pivot-coords dict word->c, words start with pivot) on arbitrary word"""
    # pivot of P = letter containing min index among word letters
    allm = 0
    for b in word: allm |= b
    m = lowbit(allm)
    pos = [i for i, b in enumerate(word) if b & m][0]
    p = word[pos]; u = word[:pos]; v = word[pos+1:]
    ru = tuple(reversed(u))
    s = 0
    # shuffles of ru and v
    L = len(ru) + len(v)
    for posset in itertools.combinations(range(L), len(ru)):
        a = []; iu = iv = 0; ps = set(posset)
        for t in range(L):
            if t in ps: a.append(ru[iu]); iu += 1
            else: a.append(v[iv]); iv += 1
        s += Pv.get((p,)+tuple(a), 0)
    return (-1)**len(u) * s

def bracket_letter(y, Pv, Sset):
    """[y, P] where P lives on index set Sset (disjoint from y)"""
    out = {}
    S = Sset | y
    if lowbit(S) & y:
        # pivot is y: coefficient on y w = coef_P(w), for all words w in letters of P
        # letters of P per partition: enumerate partitions present in Pv
        parts = set(frozenset(w) for w in Pv)
        for part in parts:
            letters = list(part)
            for perm in itertools.permutations(letters):
                c = coef_full(Pv, perm)
                if c: out[(y,)+perm] = c
    else:
        for w, c in Pv.items():
            out[w+(y,)] = -c
    return out

def br2(a, b):
    """[G_a, G_b] for blocks a,b in pivot coords"""
    if lowbit(a|b) & a: return {(a,b): 1}
    return {(b,a): -1}

def relations(S, maxS0=None):
    """all r_{i,j,S0} with {i,j} u S0 = S, i<j (r_{j,i}=-r_{i,j})"""
    idx = [1<<t for t in range(S.bit_length()) if S>>t & 1]
    res = []
    for i, j in itertools.combinations(idx, 2):
        S0 = S & ~i & ~j
        if maxS0 is not None and popcount(S0) > maxS0: continue
        v = {}
        for A in subsets(S0):
            B = S0 & ~A
            for w, c in br2(A|i, B|j).items():
                v[w] = (v.get(w, 0) + c) % P
        v = {w: c for w, c in v.items() if c}
        if v: res.append(v)
    return res

def blocks_of(S, small=False):
    return [b for b in subsets(S) if b and b != S and (not small or popcount(b) <= 2)]

def compute(n, small=False):
    """ideal I(S) for all S subsets of [n]; graded by k = number of letters.
       returns dict S -> {k: list of basis rows}"""
    I = {}
    full = (1<<n) - 1
    for size in range(2, n+1):
        for S in range(1, full+1):
            if popcount(S) != size: continue
            # only need canonical S = [size] for speed? compute all (needed as subsets)
            if S != (1<<size)-1:
                continue
            ech = {}
            def addv(v):
                if not v: return
                k = len(next(iter(v)))
                e = ech.setdefault(k, Echelon())
                e.add(v)
            for r in relations(S, 1 if small else None):
                addv(r)
            for y in blocks_of(S, small):
                Sp = S & ~y
                for k, rows in I_get(I, Sp).items():
                    for row in rows:
                        addv(bracket_letter(y, row, Sp))
            I[S] = {k: e.rows() for k, e in ech.items()}
    return I

def relabel_map(Scan, S):
    """order preserving map from bits of Scan=[m] to bits of S"""
    tgt = [1<<t for t in range(S.bit_length()) if S>>t & 1]
    def f(b):
        r = 0; t = 0
        while b:
            if b & 1: r |= tgt[t]
            b >>= 1; t += 1
        return r
    return f

def I_get(I, S):
    m = popcount(S)
    if m < 2: return {}
    can = (1<<m)-1
    base = I[can]
    if S == can: return base
    f = relabel_map(can, S)
    return {k: [{tuple(f(b) for b in w): c for w, c in row.items()} for row in rows] for k, rows in base.items()}

def is_small(w): return all(popcount(b) <= 2 for b in w)

def free_small_dim(n, k):
    # partitions of [n] into k blocks of size 1 or 2, times (k-1)!
    # number with t pairs: n = k + t
    t = n - k
    if t < 0 or k < 1: return 0
    from math import comb, factorial
    npart = comb(n, 2*t) * factorial(2*t) // (factorial(t) * 2**t)
    return npart * factorial(k-1)

def stirling2(n, k):
    from functools import lru_cache
    @lru_cache(None)
    def S2(n, k):
        if n == k: return 1
        if k == 0 or k > n: return 0
        return k*S2(n-1, k) + S2(n-1, k-1)
    return S2(n, k)

if __name__ == '__main__':
    from math import factorial
    N = int(sys.argv[1])
    t0 = time.time()
    IK = compute(N, small=False)
    IL = compute(N, small=True)
    print("n k | dimFree dimI dimK | dimFreeSmall dim(I cap small) dimI_l  dim l  dim h  kernel")
    for n in range(2, N+1):
        S = (1<<n)-1
        for k in range(1, n+1):
            rowsK = IK[S].get(k, [])
            dI = len(rowsK)
            # rank of projection of I onto big coordinates
            e = Echelon()
            for row in rowsK:
                e.add({w: c for w, c in row.items() if not is_small(w)})
            dcap = dI - e.rank()
            dIl = len(IL[S].get(k, []))
            fs = free_small_dim(n, k)
            fr = stirling2(n, k)*factorial(k-1)
            if fs == 0 and fr == 0: continue
            print(n, k, '|', fr, dI, fr-dI, '|', fs, dcap, dIl, fs-dIl, fs-dcap, dcap-dIl, flush=True)
    print("time", time.time()-t0)
