import pickle
import sympy as sym
import numpy as np
from sympy import oo
from sympy import hermite
from sympy import factorial


def f(z, sigma=1):
    return sym.exp(-z**2 / (2 * sigma)) / sym.sqrt(2 * sym.pi * sigma)

N = 10

z = sym.symbols('z')


def hermitehe(n, z):
    return 2**(-float(n)/2) \
        * hermite(n, z / sym.sqrt(2)) / sym.sqrt(factorial(n))

# Triple correlations of Hermite polynomials
e = np.zeros((N, N, N))
for i in range(N):
    for j in range(N):
        for k in range(N):
            e[i, j, k] = sym.integrate(
                hermitehe(i, z) * hermitehe(j, z) * hermitehe(k, z) * f(z),
                (z, -oo, oo))

for k in range(N):
    print "Eqn {0}".format(k)
    print e[:, :, k]

with open('e_ijk.pickle', 'w') as f:
    pickle.dump(e, f)
