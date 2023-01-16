#from .bardell_pre_integrated import *
import sys
sys.path.append('D:/repositories/composipy/composipy/bardell_pre_integrated')
from _fxi import *
from _K_KG_IJKL import * 


__all__ = ['i_d_fxi', 'i_dd_fxi', 'calc_k']


def fxi(i, xi):
    '''returns the ith approximation function evaluated at xi'''
    fxi = FXI[i]
    return eval(fxi)


def d_fxi(i, xi):
    '''returns the first derivative of  ith approximation function evaluated at xi'''
    d_fxi = D_FXI[i]
    return eval(d_fxi)


def dd_fxi(i, xi):
    '''returns the second derivative of ith approximation function evaluated at xi'''
    dd_fxi = DD_FXI[i]
    return eval(dd_fxi)


def _i_d_fxi(ij, xi):
    '''returns the integral of the first derivative of ith approximation function evaluated at xi'''
    i_d_fxi = I_D_FXI[ij]
    return eval(i_d_fxi)


def _i_dd_fxi(ij, xi):
    '''returns the integral of the second derivative of ith approximation function evaluated at xi'''
    i_dd_fxi = I_DD_FXI[ij]
    return eval(i_dd_fxi)


def i_d_fxi(ij, lbound, ubound):
    '''
    Definite integral of the first derivative of ijth approximation function evaluated at xi

    Parameters
    ----------
    ij : tuple
        index correspondent to the desired function
    lbound : float
        Lower bound of the definite integral
    ubound : float
        Upper bound of the definite integral
    
    Returns
    -------
    definite_integral : float
        integral result   
    '''
    definite_integral = _i_d_fxi(ij, ubound) - _i_d_fxi(ij, lbound)
    return definite_integral


def i_dd_fxi(ij, lbound, ubound):
    '''
    Definite integral of the second derivative of ijth approximation function evaluated at xi

    Parameters
    ----------
    ij : tuple
        index correspondent to the desired function
    lbound : float
        Lower bound of the definite integral
    ubound : float
        Upper bound of the definite integral
    
    Returns
    -------
    definite_integral : float
        integral result   
    '''
    definite_integral = _i_dd_fxi(ij, ubound) - _i_dd_fxi(ij, lbound)
    return definite_integral


def calc_k(F, i, j, a, b):
    A11, A12, A16, A22, A26, A66 = F[0, 0], F[0, 1], F[0, 2], F[1, 1], F[1, 2], F[2, 2]
    B11, B12, B16, B22, B26, B66 = F[0, 3], F[0, 4], F[0, 5], F[1, 4], F[1, 5], F[2, 5]
    D11, D12, D16, D22, D26, D66 = F[3, 3], F[3, 4], F[3, 5], F[4, 4], F[4, 5], F[5, 5]

# F = sp.Matrix([
#     [A11, A12, A16, B11, B12, B16],
#     [A12, A22, A26, B12, B22, B26],
#     [A16, A26, A66, B16, B26, B66],
#     [B11, B12, B16, D11, D12, D16],
#     [B12, B22, B26, D12, D22, D26],
#     [B16, B26, B66, D16, D26, D66]])
    kij = eval(K_IJKL[(i, j)])
    return kij