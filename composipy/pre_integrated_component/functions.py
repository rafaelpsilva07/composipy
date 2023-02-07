import sys

import numpy as np

from composipy.pre_integrated_component._ii_F import *
from composipy.pre_integrated_component._S import *


__all__ = ['ii_ff', 'ii_fxi_f', 'ii_fxi_fxi', 'ii_fxixi_f', 'ii_fxixi_fxi', 'ii_fxixi_fxixi']
__all__.extend(['ii_f_fxi', 'ii_f_fxixi', 'ii_fxi_fxixi'])
__all__.extend(['fxi', 'sxieta'])


def ii_ff(ik):
    return eval(ii_FF[ik])


def ii_fxi_f(ik):
    return eval(ii_FXI_F[ik])


def ii_f_fxi(ik):
    ki = ik[::-1]
    return eval(ii_FXI_F[ki])


def ii_fxi_fxi(ik):
    return eval(ii_FXI_FXI[ik])


def ii_fxixi_f(ik):
    return eval(ii_FXIXI_F[ik])


def ii_f_fxixi(ik):
    ki = ik[::-1]
    return eval(ii_FXIXI_F[ki])


def ii_fxixi_fxi(ik):
    return eval(ii_FXIXI_FXI[ik])


def ii_fxi_fxixi(ik):
    ki = ik[::-1]
    return eval(ii_FXIXI_FXI[ki])


def ii_fxixi_fxixi(ik):
    return eval(ii_FXIXI_FXIXI[ik])
    

def fxi(i, xi):
    return eval(FXI[i])
    

def sxieta_ij(ij, xi, eta):
    return eval(S[ij])


def sxieta(s_idx, xi, eta):
    s = []
    for ij in s_idx:
        s.append(
            sxieta_ij(ij, xi, eta)
            )
    s = np.array(s)

    return s

