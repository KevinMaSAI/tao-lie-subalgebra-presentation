"""Search for a letter order making the quadratic relations of U(g) (or of its dual A^!) a
Groebner basis (PBW algebra).  Quadratic algebra T(W)/(R) is PBW for deglex order iff
dim(W R + R W)_3 == #(degree-3 words containing a leading word of R) (diamond lemma; all
ambiguities are of degree 3).  Indices 1..3 suffice: everything depends only on the
order pattern of indices (<=3 distinct indices in a degree-3 word).  Exact over Q."""
import itertools, sys
from fractions import Fraction

def rank_and_leads(rows, key):
    """rows: list of dict word->coef. returns reduced echelon with leading words (max by key)"""
    piv = {}
    for r in rows:
        v = {w: Fraction(c) for w, c in r.items() if c}
        while v:
            lw = max(v, key=key)
            if lw in piv:
                f = v[lw]
                for w, c in piv[lw].items():
                    v[w] = v.get(w, 0) - f*c
                    if v[w] == 0: del v[w]
            else:
                f = v[lw]
                piv[lw] = {w: c/f for w, c in v.items()}
                break
    return piv

def Ug_relations(N):
    R = []
    idx = range(1, N+1)
    for i, j in itertools.combinations(idx, 2):
        R.append({(('D',i),('D',j)): 1, (('D',j),('D',i)): -1})
        R.append({(('E',i),('E',j)): 1, (('E',j),('E',i)): -1})
        R.append({(('D',i),('E',j)): 1, (('E',j),('D',i)): -1, (('D',j),('E',i)): -1, (('E',i),('D',j)): 1})
    return R

def dual_relations(N):
    """A^! = T(W*)/(R^perp); R^perp = annihilator of R in W*⊗W* (pairing word-by-word)."""
    letters = [(t, i) for t in 'DE' for i in range(1, N+1)]
    words = [(a, b) for a in letters for b in letters]
    R = Ug_relations(N)
    # nullspace of R as functionals: vectors v on words with sum_w r[w] v[w] = 0 for all r
    # compute via Fraction gaussian elimination
    cols = {w: k for k, w in enumerate(words)}
    M = [[Fraction(r.get(w, 0)) for w in words] for r in R]
    # rref
    piv_cols = []; rI = 0
    M = [row[:] for row in M]
    for c in range(len(words)):
        p = next((r for r in range(rI, len(M)) if M[r][c] != 0), None)
        if p is None: continue
        M[rI], M[p] = M[p], M[rI]
        f = M[rI][c]; M[rI] = [x/f for x in M[rI]]
        for r in range(len(M)):
            if r != rI and M[r][c] != 0:
                g = M[r][c]; M[r] = [a - g*b for a, b in zip(M[r], M[rI])]
        piv_cols.append(c); rI += 1
        if rI == len(M): break
    free = [c for c in range(len(words)) if c not in piv_cols]
    basis = []
    for fc in free:
        v = {words[fc]: Fraction(1)}
        for r, pc in enumerate(piv_cols):
            if M[r][fc] != 0: v[words[pc]] = -M[r][fc]
        basis.append(v)
    return basis

def check(R, N, key):
    letters = [(t, i) for t in 'DE' for i in range(1, N+1)]
    piv = rank_and_leads(R, key)
    leads = set(piv)
    rows3 = []
    for r in piv.values():
        for x in letters:
            rows3.append({(x,)+w: c for w, c in r.items()})
            rows3.append({w+(x,): c for w, c in r.items()})
    piv3 = rank_and_leads(rows3, key)
    nonnormal = sum(1 for w in itertools.product(letters, repeat=3)
                    if (w[0], w[1]) in leads or (w[1], w[2]) in leads)
    return len(piv3), nonnormal, leads, piv

def letter_keys():
    out = {}
    for prim in ('type', 'index'):
        for tord in ('DE', 'ED'):
            for isign in (1, -1):
                def lk(l, prim=prim, tord=tord, isign=isign):
                    t = tord.index(l[0]); i = isign*l[1]
                    return (t, i) if prim == 'type' else (i, t)
                out[(prim, tord, isign)] = lk
    return out

if __name__ == '__main__':
    N = 3
    for name, R in (('U(g)', Ug_relations(N)), ('A^!', dual_relations(N))):
        print("==", name, "number of quadratic relations:", len(R))
        for kname, lk in letter_keys().items():
            key = lambda w, lk=lk: tuple(lk(l) for l in w)   # deglex (all words same degree)
            d3, nn, leads, piv = check(R, N, key)
            print(kname, "dim I_3 =", d3, " non-normal words =", nn, "PBW" if d3 == nn else "")
