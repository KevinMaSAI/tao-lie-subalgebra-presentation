"""
Method B (independent check): work directly in g = FreeLie(D_i,E_i)/([D_i,D_j],[E_i,E_j],[D_i,E_j]-[D_j,E_i]),
multilinear in indices (each index used once, as D or as E).
Letters: (bit, t) with t=0 for D, 1 for E.  Pivot = letter with smallest index.
Free Lie elements stored by coefficients on words starting with pivot.
Computes: dim g_S by D-degree, and dim of image of h = <E_i, F_ij=[D_i,E_j]> in g_S
graded by (number of E's and F's).
"""
import sys, itertools, time
P = 2147483647
def popcount(x): return bin(x).count('1')
def lowbit(x): return x & -x

class Echelon:
    def __init__(self): self.piv = {}
    def add(self, v):
        v = {c: x % P for c, x in v.items() if x % P}
        while v:
            c = max(v); r = self.piv.get(c)
            if r is None:
                inv = pow(v[c], P-2, P)
                self.piv[c] = {cc: x*inv % P for cc, x in v.items()}
                return True
            f = v[c]
            for cc, x in r.items():
                y = (v.get(cc, 0) - f*x) % P
                if y: v[cc] = y
                else: v.pop(cc, None)
        return False
    def rank(self): return len(self.piv)
    def rows(self): return list(self.piv.values())

def coef_full(Pv, word):
    m = min(l[0] for l in word)
    pos = [i for i, l in enumerate(word) if l[0] == m][0]
    p = word[pos]; u = word[:pos]; v = word[pos+1:]
    ru = tuple(reversed(u)); L = len(ru)+len(v); s = 0
    for posset in itertools.combinations(range(L), len(ru)):
        a = []; iu = iv = 0; ps = set(posset)
        for t in range(L):
            if t in ps: a.append(ru[iu]); iu += 1
            else: a.append(v[iv]); iv += 1
        s += Pv.get((p,)+tuple(a), 0)
    return (-1)**len(u) * s % P

def add_into(out, v, f=1):
    for w, c in v.items():
        y = (out.get(w, 0) + f*c) % P
        if y: out[w] = y
        else: out.pop(w, None)

def br_letter(y, Pv):
    """[y, P], y a letter (bit,t) whose index is not in P"""
    out = {}
    if not Pv: return out
    m = min(l[0] for l in next(iter(Pv)))
    if y[0] < m:
        letters_sets = set(frozenset(w) for w in Pv)
        for ls in letters_sets:
            for perm in itertools.permutations(ls):
                c = coef_full(Pv, perm)
                if c: out[(y,)+perm] = c
    else:
        for w, c in Pv.items(): out[w+(y,)] = (-c) % P
    return out

def letter(b, t): return (b, t)
def gen_vec(l): return {(l,): 1}
def br_gen(a, b):  # [a,b] letters
    return br_letter(a, gen_vec(b))

def bits(S): return [1<<t for t in range(S.bit_length()) if S>>t & 1]

def relations(S):
    i, j = bits(S)
    D, E = 0, 1
    r1 = br_gen((i, D), (j, D)); r2 = br_gen((i, E), (j, E))
    r3 = dict(br_gen((i, D), (j, E))); add_into(r3, br_gen((j, D), (i, E)), -1)
    return [r1, r2, r3]

def Fvec(i, j):  # F_ij = [D_i, E_j]
    return br_gen((i, 0), (j, 1))

def brF(i, j, Pv):  # [[D_i,E_j], P] = [D_i,[E_j,P]] - [E_j,[D_i,P]]
    out = dict(br_letter((i, 0), br_letter((j, 1), Pv)))
    add_into(out, br_letter((j, 1), br_letter((i, 0), Pv)), -1)
    return out

def relabel(v, f):
    return {tuple((f[l[0]], l[1]) for l in w): c for w, c in v.items()}

def main(N):
    I = {}; H = {}
    for n in range(1, N+1):
        S = (1<<n)-1
        eI = {}; eH = {}
        def addI(v):
            if not v: return
            a = sum(1 for l in next(iter(v)) if l[1] == 0)
            eI.setdefault(a, Echelon()).add(v)
        def addH(v, kk):
            if not v: return
            eH.setdefault(kk, Echelon()).add(v)
        if n == 1:
            addH(gen_vec((1, 1)), 1)
        if n == 2:
            for r in relations(S): addI(r)
            addH(Fvec(1, 2), 1)
        # ideal: [y, I(S\y)]
        for b in bits(S):
            Sp = S & ~b
            if Sp == 0: continue
            tgt = bits(Sp); m = len(tgt)
            f = {1<<t: tgt[t] for t in range(m)}
            for a, rows in I.get(m, {}).items():
                for row in rows:
                    rr = relabel(row, f)
                    for t in (0, 1): addI(br_letter((b, t), rr))
            # h: [E_b, h(S\b)]
            for kk, rows in H.get(m, {}).items():
                for row in rows: addH(br_letter((b, 1), relabel(row, f)), kk+1)
        # h: [F_{b c}, h(S\{b,c})]
        bl = bits(S)
        for b, c in itertools.combinations(bl, 2):
            Sp = S & ~b & ~c
            if Sp == 0: continue
            tgt = bits(Sp); m = len(tgt)
            f = {1<<t: tgt[t] for t in range(m)}
            for kk, rows in H.get(m, {}).items():
                for row in rows: addH(brF(b, c, relabel(row, f)), kk+1)
        I[n] = {a: e.rows() for a, e in eI.items()}
        H[n] = {kk: e.rows() for kk, e in eH.items()}
        # report
        from math import comb, factorial
        for a in range(0, n+1):
            fr = comb(n, a)*factorial(n-1)
            dI = len(I[n].get(a, []))
            print("g", n, "Ddeg", a, "dim g =", fr-dI, flush=True)
        for kk, rows in sorted(H[n].items()):
            # image dim = rank(h + I) - rank(I) in D-degree a = n - kk
            a = n - kk
            e = Echelon()
            for r in I[n].get(a, []): e.add(r)
            r0 = e.rank()
            for r in rows: e.add(r)
            print("h", n, "k", kk, "dim h image =", e.rank()-r0, flush=True)
        # keep H rows reduced modulo I to keep small: (not needed for correctness)
main(int(sys.argv[1]))
