import sys
import numpy as np
import sympy as sp

sys.path.append('D:/repositories/composipy/composipy/bardell_pre_integrated')
from _sym_fxi_geta import *


m = 15
n = 15

uf_xi = FXI.copy()
vf_xi = FXI.copy()
wf_xi = FXI.copy()
uf_eta = GETA.copy()
vf_eta = GETA.copy()
wf_eta = GETA.copy()


a, b = sp.symbols(['a', 'b'])
xi, eta = sp.symbols(['xi', 'eta'])

A11, A12, A16, A22, A26, A66 = sp.symbols(['A11', 'A12', 'A16', 'A22', 'A26', 'A66']) # jones (5.16)
B11, B12, B16, B22, B26, B66 = sp.symbols(['B11', 'B12', 'B16', 'B22', 'B26', 'B66'])
D11, D12, D16, D22, D26, D66 = sp.symbols(['D11', 'D12', 'D16', 'D22', 'D26', 'D66'])

Nxx, Nyy, Nxy = sp.symbols(['N_xx', 'N_yy', 'N_xy'])
xi_lbound, xi_ubound,  = sp.symbols(['xi_lbound', 'xi_ubound'])
eta_lbound, eta_ubound,  = sp.symbols(['eta_lbound', 'eta_ubound'])

Su, Sv, Sw = [[]], [[]], [[]]
for i in range(m):
    for j in range(n):
        Su[0].append(uf_xi[i]*uf_eta[j])
        Sv[0].append(vf_xi[i]*vf_eta[j])
        Sw[0].append(wf_xi[i]*wf_eta[j])
Su, Sv, Sw = sp.Matrix(Su), sp.Matrix(Sv), sp.Matrix(Sw)


B0_11 = np.array((2/a) * sp.diff(Su, xi))
B0_22 = np.array((2/b) * sp.diff(Sv, eta))
B0_31 = np.array((2/b)*sp.diff(Su, eta))
B0_32 = np.array((2/a)*sp.diff(Sv, xi))
B0_43 = np.array(-(4/a**2) * sp.diff(Sw, xi, xi))
B0_53 = np.array(-(4/b**2) * sp.diff(Sw, eta, eta))
B0_63 = np.array(-2*(4/(a*b)) * sp.diff(Sw, xi, eta))
Z = np.zeros(m*n).reshape(1, m*n)

B0_kappa = sp.Matrix([
    np.hstack([*B0_11, *Z, *Z]),
    np.hstack([*Z, *B0_22, *Z]),
    np.hstack([*B0_31, *B0_32, *Z]),
    np.hstack([*Z, *Z, *B0_43]),
    np.hstack([*Z, *Z, *B0_53]),
    np.hstack([*Z, *Z, *B0_63])
])


F = sp.Matrix([
    [A11, A12, A16, B11, B12, B16],
    [A12, A22, A26, B12, B22, B26],
    [A16, A26, A66, B16, B26, B66],
    [B11, B12, B16, D11, D12, D16],
    [B12, B22, B26, D12, D22, D26],
    [B16, B26, B66, D16, D26, D66]
])


G = sp.Matrix([
    np.hstack([*Z, *Z, *sp.diff(Sw, xi)*(2/a)]),
    np.hstack([*Z, *Z, *sp.diff(Sw, eta)*(2/b)])
 ])

    
G = sp.Matrix(G)

N = sp.Matrix([[Nxx, Nxy],
                [Nxy, Nyy]])



KG = sp.integrate(
     sp.integrate(
        G.T * N * G, (xi, -1, 1)), (eta, -1, 1)
) * (a*b/4)


K = sp.integrate(
    sp.integrate(
        B0_kappa.T * F * B0_kappa, (xi, -1, 1)), (eta, -1, 1)
) * (a*b/4)


txt_kg = 'KG_IJKL = {\n'
txt_k = 'K_IJKL = {\n'

size = 3 * m * n

for i in range(size):
    for j in range(size):
        txt_kg += f'    {(i, j)}: \'' + str(KG[i, j]) + '\',\n'
        txt_k +=  f'    {(i, j)}: \'' + str(K[i, j]) + '\',\n'

txt_kg += '}\n\n\n'
txt_k += '}\n\n\n'
txt = txt_k + txt_kg

with open('_K_KG_IJKL.py', 'w') as f:
    f.write(txt)
