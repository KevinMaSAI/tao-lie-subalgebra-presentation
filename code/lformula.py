# dim l(n,k) and dim K(n,k) (proven formulas) for comparison
import sympy as sp, sys
x, y, z = sp.symbols('x y z'); M = int(sys.argv[1])
Y = sum(x**w/sp.factorial(w)*sum(y**(s-1)*z**(w-s+1) for s in range(2, w+1)) for w in range(2, M+1))
Y2 = sum(x**w/sp.factorial(w)*y*z**(w-1) for w in range(2, M+1))
K = sp.expand(x*z + sp.series(-sp.log(1-Y), x, 0, M+1).removeO())
L = sp.expand(x*z + sp.series(-sp.log(1-Y2), x, 0, M+1).removeO())
for n in range(2, M+1):
    kc = sp.Poly(sp.expand(K.coeff(x, n)*sp.factorial(n)), y, z); lc = sp.Poly(sp.expand(L.coeff(x, n)*sp.factorial(n)), y, z)
    print(n, "K:", [kc.coeff_monomial(y**(n-k)*z**k) for k in range(1, n+1)], " l:", [lc.coeff_monomial(y**(n-k)*z**k) for k in range(1, n+1)],
          " K-l:", [kc.coeff_monomial(y**(n-k)*z**k)-lc.coeff_monomial(y**(n-k)*z**k) for k in range(1, n+1)])
