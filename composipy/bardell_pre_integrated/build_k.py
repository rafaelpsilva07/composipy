import sys
sys.path.append('D:/repositories/composipy/composipy/bardell_pre_integrated')

from functions import *


def calc_K11_ijkl(a, b, ui, uj, uk, ul, A11, A16, A66):
    uik = (ui, uk)
    ujl = (uj, ul)
    uxi_ij_uxi_kl = ii_fxi_fxi(uik) * ii_ff(ujl)
    uxi_ij_ueta_kl = ii_fxi_f(uik) * ii_f_fxi(ujl)
    ueta_ij_uxi_kl = ii_f_fxi(uik) * ii_fxi_f(ujl)
    ueta_ij_ueta_kl = ii_ff(uik) * ii_fxi_fxi(ujl)

    k11 = (
        A11 * (b/a) *  uxi_ij_uxi_kl
        + A16 * (uxi_ij_ueta_kl+ueta_ij_uxi_kl)
        + A66 * (b/a) * ueta_ij_ueta_kl)
    return k11


def calc_k12_ijkl(a, b, ijkl, A12, A16, A26, A66):
    klij = (ijkl[2], ijkl[3], ijkl[0],ijkl[1])
    k12 = (
        A12 * ii_sxi_seta(ijkl)
        + A16 * (b/a) * ii_sxi_sxi(ijkl)
        + A26 * (a/b) * ii_sxi_sxi(ijkl)
        + A66 * ii_sxi_seta(klij)
    )
    return k12


def calc_k13_ijkl():#(a, b, ijkl, B11, B12, B16, B26, B66):
    klij = (ijkl[2], ijkl[3], ijkl[0],ijkl[1])
    k13 = (
        -2 * (b/a**2) * B11 * ii_sxixi_sxi(klij)
        -(2/b) * B12 * ii_setaeta_sxi(klij)
        -(2/a) * B16 * (2*ii_sxieta_sxi(klij)+ii_setaeta_sxi(klij))
        -2 * (a/b**2) * B26 * ii_sxixi_sxi(klij)
        -(4/b) * B66 * ii_sxieta_sxi(klij)
    )
    return 0


def calc_k21_ijkl(a, b, uk, ul, vi, vj, A12, A16, A26, A66):
    vij_ukl = (vi, vj, uk, ul)
    k21 = (
        A12 * ii_sxi_seta(vij_ukl)
        + A16 * (a/b) * ii_sxi_sxi(vij_ukl)
        + A26 * (a/b) * ii_sxi_sxi(vij_ukl)
        + A66 * ii_sxi_seta(vij_ukl)
    )
    return k21


def calc_k22_ijkl(a, b, vi, vj, vk, vl, A22, A26, A66):
    vij_vkl = (vi, vj, vk, vl)
    k22 = (
        A22 * (a/b) * ii_sxi_sxi(vij_vkl)
        + A26 * (ii_sxi_seta(vij_vkl)+ii_sxi_seta(vij_vkl))
        + A66 * (b/a) * ii_sxi_sxi(vij_vkl)
    )
    return k22


def calc_k23_ijkl(): #(a, b, vi, vj, wk, wl, B12, B16, B22, B26, B66):
    return 0


def calc_k31_ijkl():
    return


def calc_k32_ijkl():
    return 0


def calc_k33_ijkl(a, b, wi, wj, wk, wl, D11, D12, D22, D16, D26, D66):
    wijkl = (wi, wj, wk, wl)
    wklij = (wk, wl, wi, wj)
    k33 = (
        4*(b/a**3) * D11 * ii_sxixi_sxixi(wijkl)
        + (4/(a*b)) * D12 * (ii_sxixi_setaeta(wijkl)+ii_sxixi_setaeta(wklij))
        + 4*(a/b**3) * D22 * ii_sxixi_sxixi(wijkl)
        + 8/a**2 * D16 * (ii_sxixi_sxieta(wklij))
        + 8/b**2 * D26 * (ii_sxixi_sxieta(wklij)+ii_sxixi_sxieta(wijkl))
        + (16/(a*b)) * D66 * ii_sxieta_sxieta(wijkl)
    )
    return k33





