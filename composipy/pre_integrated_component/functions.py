import sys

import numpy as np

from composipy.pre_integrated_component._ii_F import *
from composipy.pre_integrated_component._S import *


__all__ = ['integrate_ff', 'integrate_fxi_f', 'integrate_fxi_fxi', 'integrate_fxixi_f', 'integrate_fxixi_fxi', 'integrate_fxixi_fxixi']
__all__.extend(['integrate_f_fxi', 'integrate_f_fxixi', 'integrate_fxi_fxixi'])
__all__.extend(['calc_f', 'calc_fxi', 'sxieta'])


def integrate_ff(ik):
    return eval(ii_FF[ik])


def integrate_fxi_f(ik):
    return eval(ii_FXI_F[ik])


def integrate_f_fxi(ik):
    ki = ik[::-1]
    return eval(ii_FXI_F[ki])


def integrate_fxi_fxi(ik):
    return eval(ii_FXI_FXI[ik])


def integrate_fxixi_f(ik):
    return eval(ii_FXIXI_F[ik])


def integrate_f_fxixi(ik):
    ki = ik[::-1]
    return eval(ii_FXIXI_F[ki])


def integrate_fxixi_fxi(ik):
    return eval(ii_FXIXI_FXI[ik])


def integrate_fxi_fxixi(ik):
    ki = ik[::-1]
    return eval(ii_FXIXI_FXI[ki])


def integrate_fxixi_fxixi(ik):
    return eval(ii_FXIXI_FXIXI[ik])
    

def calc_f(i, xi):
    return eval(F[i])


def calc_fxi(i, xi):
    return eval(F[i])
    

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

