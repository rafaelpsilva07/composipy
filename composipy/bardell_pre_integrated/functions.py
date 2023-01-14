#from .bardell_pre_integrated import *
import sys
sys.path.append('D:/repositories/composipy/composipy/bardell_pre_integrated')
from _fxi_dfxi_ddfxi import *


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

