#from .bardell_pre_integrated import *
import sys
sys.path.append('D:/repositories/composipy/composipy/bardell_pre_integrated')
from _fxi import *


__all__ = ['i_d_fxi', 'i_dd_fxi']


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
