# Koszul prediction vs K^top formula vs Eulerian numbers
import sympy as sp
x, y, z = sp.symbols('x y z')
M = 10
Ashrek = sum((-x)**n/sp.factorial(n)*sum(y**a*z**(n-a) for a in range(n+1)) for n in range(M+1))
# G = -log A^!(-x)  (series in x)
G = sp.series(-sp.log(Ashrek), x, 0, M+1).removeO()
Y = sum(x**w/sp.factorial(w)*sum(y**(s-1)*z**(w-s+1) for s in range(2, w+1)) for w in range(2, M+1))
Ktop = x*(y+z) + sp.series(-sp.log(1-Y), x, 0, M+1).removeO()
d = sp.expand(G - Ktop)
print("G - (x(y+z) - log(1-Y)) =", d)
def eulerian(n, k):
    from math import comb
    return sum((-1)**j*comb(n+1, j)*(k+1-j)**n for j in range(k+2))
for n in range(2, M+1):
    c = sp.Poly(sp.expand(G.coeff(x, n)*sp.factorial(n)), y, z)
    # coefficient by E-degree k (z power)
    row = [c.coeff_monomial(y**(n-k)*z**k) for k in range(1, n)]
    print(n, row, [eulerian(n-1, k-1) for k in range(1, n)])
