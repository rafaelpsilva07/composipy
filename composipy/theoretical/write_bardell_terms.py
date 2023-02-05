import sympy as sp
import numpy as np


n_terms = 10
xi, r1, r2, t1, t2  = sp.symbols(['xi', 'r1', 'r2', 't1', 't2'])
f1 = 1/2 - 3/4*xi + 1/4*xi**3
f2 = 1/8 - 1/8*xi - 1/8*xi**2 + 1/8*xi**3
f3 = 1/2 + 3/4*xi - 1/4*xi**3
f4 = -1/8 - 1/8*xi + 1/8*xi**2 + 1/8*xi**3

f_xi = {
    0: t1 * f1,
    1: r1 * f2,
    2: t2 * f3,
    3: r2 * f4
} # first terms of set of functions multiplied by constraint flag


def _orth_polynomials(r):   
    '''
    r term from the set of functions used by Bardell (1991)

    Parameters
    ----------
    r : int

    Returns
    -------
    f_xi_r : sympy equation
    r term from the set of functions.
    '''
    global xi
    nmax = r//2 + 1 #denotes its own integer part (Bardel, 1991).
    f_xi_r = 0

    for n in range(0, nmax):
        f_xi_r += (
            (-1)**n * sp.factorial2(2*r-2*n-7)
            / (2**n * sp.factorial(n) * sp.factorial(r-2*n-1))
            * xi**(r-2*n-1)
            )
    f_xi_r = sp.sympify(f_xi_r)
    return f_xi_r


# Complete set of bardell functions
for i in range(4, n_terms):
    f_xi[i] = _orth_polynomials(i+1)


print(f_xi)