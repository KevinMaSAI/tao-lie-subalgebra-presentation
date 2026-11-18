"""Sanity check that the degree-3 test of pbw_check.py is not vacuous: for random quadratic
relations (same number as for U(g), indices 1..3, deglex with the order of Section 7) the
test should typically FAIL (dim I_3 > number of non-normal words; note dim I_3 >= #non-normal
always holds, because the normal words span modulo I_3)."""
import random, itertools
import pbw_check as pc
random.seed(20261115)
N = 3
letters = [(t, i) for t in 'DE' for i in range(1, N+1)]
words = [(a, b) for a in letters for b in letters]
lk = lambda l: (l[1], 'DE'.index(l[0]))
key = lambda w: tuple(lk(l) for l in w)
for trial in range(3):
    R = [{w: random.randint(-3, 3) for w in random.sample(words, 6)} for _ in range(9)]
    d3, nn, leads, piv = pc.check(R, N, key)
    print("trial", trial, "dim I_3 =", d3, " non-normal words =", nn, "PBW" if d3 == nn else "NOT PBW")
