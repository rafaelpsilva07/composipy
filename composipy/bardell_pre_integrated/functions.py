import sys
sys.path.append('D:/repositories/composipy/composipy/bardell_pre_integrated')
from _ii_F import *


__all__ = ['ii_ff', 'ii_fxi_f', 'ii_fxi_fxi', 'ii_fxixi_f', 'ii_fxixi_fxi', 'ii_fxixi_fxixi']
__all__.extend(['ii_f_fxi', 'ii_f_fxixi', 'ii_fxi_fxixi'])


ii_FF, ii_FXI_F, ii_FXI_FXI, ii_FXIXI_F, ii_FXIXI_FXI, ii_FXIXI_FXIXI

def ii_ff(ik):
    return ii_FF[ik]


def ii_fxi_f(ik):
    return ii_FXI_F[ik]


def ii_f_fxi(ik):
    ki = ik[::-1]
    return ii_FXI_F[ki]


def ii_fxi_fxi(ik):
    return ii_FXI_FXI[ik]


def ii_fxixi_f(ik):
    return ii_FXIXI_F[ik]


def ii_f_fxixi(ik):
    ki = ik[::-1]
    return ii_FXIXI_F[ki]


def ii_fxixi_fxi(ik):
    return ii_FXIXI_FXI[ik]


def ii_fxi_fxixi(ik):
    ki = ik[::-1]
    return ii_FXIXI_FXI[ki]


def ii_fxixi_fxixi(ik):
    return ii_FXIXI_FXIXI[ik]