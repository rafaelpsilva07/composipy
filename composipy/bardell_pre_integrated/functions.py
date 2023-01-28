import sys
sys.path.append('D:/repositories/composipy/composipy/bardell_pre_integrated')
from _ii_S import *


__all__ = ['ii_sxi_sxi', 'ii_sxi_seta', 'ii_sxixi_sxi', 'ii_setaeta_sxi',
           'ii_sxieta_sxi', 'ii_sxixi_sxixi','ii_sxixi_setaeta',
            'ii_sxixi_sxieta', 'ii_sxieta_sxieta']


def ii_sxi_sxi(ijkl):
    return II_SXI_SXI[ijkl]

def ii_sxi_seta(ijkl):
    return II_SXI_SETA[ijkl]


def ii_sxixi_sxi(ijkl):
    return II_SXIXI_SXI[ijkl]


def ii_setaeta_sxi(ijkl):
    return II_SETAETA_SXI[ijkl]


def ii_sxieta_sxi(ijkl):
    return II_SXIETA_SXI[ijkl]


def ii_sxixi_sxixi(ijkl):
    return II_SXIXI_SXIXI[ijkl]


def ii_sxixi_setaeta(ijkl):
    return II_SXIXI_SETAETA[ijkl]


def ii_sxixi_sxieta(ijkl):
    return II_SXIXI_SXIETA[ijkl]


def ii_sxieta_sxieta(ijkl):
    return II_SXIETA_SXIETA[ijkl]
