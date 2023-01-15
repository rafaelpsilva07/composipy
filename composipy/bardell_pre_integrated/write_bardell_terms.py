import sympy as sp
import numpy as np
import pickle

n_terms = 10
xi = sp.symbols('xi')
f1 = 1/2 - 3/4*xi + 1/4*xi**3
f2 = 1/8 - 1/8*xi - 1/8*xi**2 + 1/8*xi**3
f3 = 1/2 + 3/4*xi - 1/4*xi**3
f4 = -1/8 - 1/8*xi + 1/8*xi**2 + 1/8*xi**3

f_xi = {
    0: f1,
    1: f2,
    2: f3,
    3: f4
} # first terms of set of functions

d_f_xi = {
    0: sp.diff(f_xi[0], xi), 
    1: sp.diff(f_xi[1], xi),
    2: sp.diff(f_xi[2], xi),
    3: sp.diff(f_xi[3], xi)
}

dd_f_xi = {
    0: sp.diff(f_xi[0], xi, xi), 
    1: sp.diff(f_xi[1], xi, xi),
    2: sp.diff(f_xi[2], xi, xi),
    3: sp.diff(f_xi[3], xi, xi)
}

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


if __name__ == '__main__':
    
    n_terms = 30

    # Complete set of bardell functions
    for i in range(4, n_terms):
        f_xi[i] = _orth_polynomials(i+1)
        d_f_xi[i] = sp.diff(f_xi[i], xi)
        dd_f_xi[i] = sp.diff(f_xi[i], xi, xi)

    # Converting simbols into strings
    txt_fxi = 'FXI = [\n'
    txt_d_fxi = 'D_FXI = [\n'
    txt_dd_fxi = 'DD_FXI = [\n'
    for k, v in f_xi.items():
        txt_fxi += '    ' + '\'' + str(v) + '\'' + ',\n'
    for k, v in d_f_xi.items():
        txt_d_fxi += '    ' + '\'' + str(v)  + '\''+ ',\n'
    for k, v in d_f_xi.items():
        txt_dd_fxi += '    ' + '\'' + str(v)  + '\'' + ',\n'

    txt_fxi += ']'
    txt_d_fxi += ']'
    txt_dd_fxi += ']'

    # Building python file
    txt = f'# This file contains the {n_terms} first bardel terms\n'
    txt += f'# The results herein are produced by {__file__}\n\n\n'
    txt += '__all__ = [\'FXI\', \'D_FXI\', \'DD_FXI\', \'I_D_FXI\', \'I_DD_FXI\']\n\n\n'
    txt += txt_fxi + '\n'*3 + txt_d_fxi + '\n'*3 + txt_dd_fxi + '\n'*3 

    
    # integration
    i_d_fxi = {}
    i_dd_fxi = {}
    for i in range(n_terms):
        for j in range(n_terms):
            i_d_fxi[(i, j)] = sp.integrate(
                                d_f_xi[i] * d_f_xi[j])
            i_dd_fxi[(i, j)] = sp.integrate(
                                dd_f_xi[i] * dd_f_xi[j])

    # Converting simbols into strings
    txt_i_d_fxi = 'I_D_FXI = {\n'
    txt_i_dd_fxi = 'I_DD_FXI = {\n'
    for k, v in i_d_fxi.items():
        txt_i_d_fxi += '    ' + str(k) + ': \'' + str(v) + '\'' + ',\n'
    for k, v in i_dd_fxi.items():
        txt_i_dd_fxi += '    ' + str(k) + ': \'' + str(v)  + '\''+ ',\n'

    txt_i_d_fxi += '}'
    txt_i_dd_fxi += '}'
    txt += txt_i_d_fxi + '\n'*3 + txt_i_dd_fxi

    with open('_fxi.py', 'w') as f:
        f.write(txt)


