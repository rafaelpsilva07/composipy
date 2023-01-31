import sys
sys.path.append('D:/repositories/composipy/composipy/pre_integrated_component')

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


def calc_k12_ijkl(a, b, ui, uj, vk, vl, A12, A16, A26, A66):
    uivk = (ui, vk)
    ujvl = (uj, vl)
    uxi_ij_veta_kl = ii_fxi_f(uivk) * ii_f_fxi(ujvl)
    uxi_ij_vxi_kl = ii_fxi_fxi(uivk) * ii_ff(ujvl)
    ueta_ij_veta_kl = ii_ff(uivk) * ii_fxi_fxi(ujvl)
    ueta_ij_vxi_kl = ii_f_fxi(uivk) * ii_fxi_f(ujvl)

    k12 = (
        A12 * uxi_ij_veta_kl
        + A16 * (b/a) * uxi_ij_vxi_kl
        + A26 * (a/b) * ueta_ij_veta_kl
        + A66 * ueta_ij_vxi_kl
    )
    return k12


def calc_k13_ijkl():#(a, b, ijkl, B11, B12, B16, B26, B66):

    k13 = 0
    #-2 * (b/a**2) * B11 #* ii_sxixi_sxi(klij)
    #-(2/b) * B12 #* ii_setaeta_sxi(klij)
    #-(2/a) * B16 #* (2*ii_sxieta_sxi(klij)+ii_setaeta_sxi(klij))
    #-2 * (a/b**2) * B26 #* ii_sxixi_sxi(klij)
    #-(4/b) * B66 #* ii_sxieta_sxi(klij)
    return 0


def calc_k21_ijkl(a, b, vi, vj, uk, ul, A12, A16, A26, A66):
    viuk = (vi, uk)
    vjul = (vj, ul)
    veta_ij_uxi_kl = ii_f_fxi(viuk) * ii_fxi_f(vjul)
    vxi_ij_uxi_kl = ii_fxi_fxi(viuk) * ii_ff(vjul)
    veta_ij_ueta_ij = ii_ff(viuk) * ii_fxi_fxi(vjul)
    vxi_ij_ueta_kl = ii_fxi_f(viuk) * ii_f_fxi(vjul)

    k21 = (
        A12 * veta_ij_uxi_kl
        + A16 * (a/b) * vxi_ij_uxi_kl
        + A26 * (a/b) * veta_ij_ueta_ij
        + A66 * vxi_ij_ueta_kl
    )
    return k21


def calc_k22_ijkl(a, b, vi, vj, vk, vl, A22, A26, A66):
    vivk = (vi, vk)
    vjvl = (vj, vl)
    veta_ij_veta_kl = ii_ff(vivk) * ii_fxi_fxi(vjvl)
    vxi_ij_veta_kl = ii_fxi_f(vivk) * ii_f_fxi(vjvl)
    veta_ij_vxi_kl = ii_f_fxi(vivk) * ii_fxi_f(vjvl)
    vxi_ij_vxi_kl = ii_fxi_fxi(vivk) * ii_ff(vjvl)

    k22 = (
        A22 * (a/b) * veta_ij_veta_kl
        + A26 * (vxi_ij_veta_kl+veta_ij_vxi_kl)
        + A66 * (b/a) * vxi_ij_vxi_kl
    )
    return k22


def calc_k23_ijkl(): #(a, b, vi, vj, wk, wl, B12, B16, B22, B26, B66):
    return 0


def calc_k31_ijkl():
    return 0


def calc_k32_ijkl():
    return 0


def calc_k33_ijkl(a, b, wi, wj, wk, wl, D11, D12, D22, D16, D26, D66):
    wiwk = (wi, wk)
    wjwl = (wj, wl)
    wxixi_ij_wxixi_kl = ii_fxixi_fxixi(wiwk) * ii_ff(wjwl)
    wxixi_ij_wetaeta_kl = ii_fxixi_f(wiwk) * ii_f_fxixi(wjwl)
    wetaeta_ij_wxixi_ij = ii_f_fxixi(wiwk) * ii_fxixi_f(wjwl)
    wetaeta_ij_wetaeta_kl = ii_ff(wiwk) * ii_fxixi_fxixi(wjwl)
    wxixi_ij_wxieta_kl = ii_fxixi_fxi(wiwk) * ii_f_fxi(wjwl)
    wxieta_ij_wxixi_kl = ii_fxi_fxixi(wiwk) * ii_fxi_f(wjwl)
    wxieta_ij_wetaeta_kl = ii_fxi_f(wiwk) * ii_fxi_fxixi(wjwl)
    wetaeta_ij_wxieta_kl = ii_f_fxi(wiwk) * ii_fxixi_fxi(wjwl)
    wxieta_ij_wxieta_kl = ii_fxi_fxi(wiwk) * ii_fxi_fxi(wjwl)

    k33 = (
        4*(b/a**3) * D11 * wxixi_ij_wxixi_kl
        + (4/(a*b)) * D12 * (wxixi_ij_wetaeta_kl+wetaeta_ij_wxixi_ij)
        + 4*(a/b**3) * D22 * wetaeta_ij_wetaeta_kl
        + 8/a**2 * D16 * (wxixi_ij_wxieta_kl+wxieta_ij_wxixi_kl)
        + 8/b**2 * D26 * (wxieta_ij_wetaeta_kl+wetaeta_ij_wxieta_kl)
        + (16/(a*b)) * D66 * wxieta_ij_wxieta_kl
    )
    return k33


def calc_kG33_ijkl(a, b, wi, wj, wk, wl, Nxx, Nyy, Nxy):
    wiwk = (wi, wk)
    wjwl = (wj, wl)
    wxi_ij_wxi_kl = ii_fxi_fxi(wiwk) * ii_ff(wjwl)
    weta_ij_weta_kl = ii_ff(wiwk) * ii_fxi_fxi(wjwl)
    weta_ij_wxi_kl = ii_f_fxi(wiwk) * ii_fxi_f(wjwl)
    wxi_ij_weta_kl = ii_fxi_f(wiwk) * ii_f_fxi(wjwl)

    kG33 = (
        Nxx * (b/a) * wxi_ij_wxi_kl
        + Nyy * (a/b) * weta_ij_weta_kl
        + Nxy * (weta_ij_wxi_kl+wxi_ij_weta_kl)
    )
    return kG33


