import sympy as sp
import numpy as np
import pickle

n_terms = 10
xi, eta, a, b = sp.symbols(['xi', 'eta', 'a', 'b'])
f1 = 1/2 - 3/4*xi + 1/4*xi**3
f2 = 1/8 - 1/8*xi - 1/8*xi**2 + 1/8*xi**3
f3 = 1/2 + 3/4*xi - 1/4*xi**3
f4 = -1/8 - 1/8*xi + 1/8*xi**2 + 1/8*xi**3

g1 = 1/2 - 3/4*eta + 1/4*eta**3
g2 = 1/8 - 1/8*eta - 1/8*eta**2 + 1/8*eta**3
g3 = 1/2 + 3/4*eta - 1/4*eta**3
g4 = -1/8 - 1/8*eta + 1/8*eta**2 + 1/8*eta**3

f_xi = {0: f1, 1: f2, 2: f3, 3: f4} # first terms of set of functions
g_eta = {0: g1, 1: g2, 2: g3, 3: g4}


def _orth_polynomials(r, sym):   
    '''
    r term from the set of functions used by Bardell (1991)

    Parameters
    ----------
    r : int
    sym : sympy.symbol
        symbol to be used in equation

    Returns
    -------
    f_sym_r : sympy equation
    r term from the set of functions.
    '''
    nmax = r//2 + 1 #denotes its own integer part (Bardel, 1991).
    f_sym_r = 0

    for n in range(0, nmax):
        f_sym_r += (
            (-1)**n * sp.factorial2(2*r-2*n-7)
            / (2**n * sp.factorial(n) * sp.factorial(r-2*n-1))
            * sym**(r-2*n-1)
            )
    f_sym_r = sp.sympify(f_sym_r)
    return f_sym_r


if __name__ == '__main__':
    
    n_terms = 30

    # Complete set of bardell functions
    for i in range(4, n_terms):
        f_xi[i] = _orth_polynomials(i+1, xi)
        g_eta[i] = _orth_polynomials(i+1, eta)

    S = {}
    B11, B22, B31, B32, B43, B53, B63 = {}, {}, {}, {}, {}, {}, {}
    for i, fi in f_xi.items():
        for j, gj in g_eta.items():
            S[(i, j)] = f_xi * g_eta
            B11[(i, j)] = (2/a) * sp.diff(S[(i, j)], xi)
            B22[(i, j)] = (2/b) * sp.diff(S[(i, j)], eta)
            B31[(i, j)] = (2/b) * sp.diff(S[(i, j)], eta)
            B32[(i, j)] = (2/a) * sp.diff(S[(i, j)], xi)
            B43[(i, j)] = -(4/a**2) * sp.diff(S[(i, j)], xi, xi)
            B53[(i, j)] = -(4/b**2) * sp.diff(S[(i, j)], eta, eta)
            B63[(i, j)] = -2 * (4/(a*b)) * sp.diff(S[(i, j)], xi, eta)

   




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


