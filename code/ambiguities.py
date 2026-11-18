"""Explicit resolution of the four overlap ambiguities of the reduction system for U(g)
(Bergman's diamond lemma).  Letters (X,i), X in 'DE'; order: index first, then D<E.
Rules (i<j):  D_j D_i -> D_i D_j ;  E_j E_i -> E_i E_j ;  E_j D_i -> D_i E_j - D_j E_i + E_i D_j.
For each overlap X_k Y_j Z_i (k>j>i) we reduce the left pair first, resp. the right pair
first, then compute full normal forms (leftmost reduction) and compare. Exact integers."""
from collections import defaultdict

def rule(a, b):
    """if ab is a leading word return its replacement (dict word->coef), else None"""
    (X, j), (Y, i) = a, b
    if j <= i: return None
    if (X, Y) == ('D', 'D'): return {(('D', i), ('D', j)): 1}
    if (X, Y) == ('E', 'E'): return {(('E', i), ('E', j)): 1}
    if (X, Y) == ('E', 'D'): return {(('D', i), ('E', j)): 1, (('D', j), ('E', i)): -1, (('E', i), ('D', j)): 1}
    return None

def add(acc, w, c):
    acc[w] = acc.get(w, 0) + c
    if acc[w] == 0: del acc[w]

def reduce_at(word, p):
    rep = rule(word[p], word[p+1])
    out = {}
    for w, c in rep.items(): add(out, word[:p] + w + word[p+2:], c)
    return out

def nf(poly):
    res = {}
    todo = dict(poly)
    while todo:
        w, c = todo.popitem()
        for p in range(len(w)-1):
            if rule(w[p], w[p+1]) is not None:
                for w2, c2 in reduce_at(w, p).items():
                    todo[w2] = todo.get(w2, 0) + c*c2
                    if todo[w2] == 0: del todo[w2]
                break
        else:
            add(res, w, c)
    return res

def tex(poly):
    def wt(w): return ''.join('%s_{%d}' % (X, i) for X, i in w)
    items = sorted(poly.items())
    s = ''
    for w, c in items:
        sg = '-' if c < 0 else '+'
        a = abs(c)
        s += ' %s %s%s' % (sg, '' if a == 1 else str(a), wt(w))
    return s.strip().lstrip('+').strip()

ok = True
for XYZ in ['DDD', 'EDD', 'EED', 'EEE']:
    w = ((XYZ[0], 3), (XYZ[1], 2), (XYZ[2], 1))
    left = reduce_at(w, 0); right = reduce_at(w, 1)
    nl, nr = nf(left), nf(right)
    print('overlap', tex({w: 1}))
    print('  left  step :', tex(left))
    print('  right step :', tex(right))
    print('  NF(left)   :', tex(nl))
    print('  NF(right)  :', tex(nr))
    print('  equal      :', nl == nr)
    ok &= (nl == nr)
print('ALL RESOLVABLE' if ok else 'FAILURE')
