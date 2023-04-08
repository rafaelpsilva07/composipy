from composipy.pre_integrated_component.functions import *

# TODO: espalhar para matriz
def calc_u2_xcte(b, k, xicte, ui, uj, uk, ul):
    jl = uj, ul
    uijkl = k * (b/2) * calc_f(ui, xicte) * calc_f(uk, xicte) * integrate_ff(jl)
    return uijkl


def calc_uaub_xcte(b, k, xictea, xicteb, uia, uja, uka, ula, uib, ujb, ukb, ulb):
    jalb = uja, ulb
    jbla = ujb, ula
    uijkl = k * (b/2) * calc_f(uia, xictea) * calc_f(ukb, xicteb) * integrate_ff(jalb)
    uijkl += k * (b/2) * calc_f(uib, xicteb) * calc_f(uka, xictea) * integrate_ff(jbla)
    return uijkl


# TODO: espalhar para matrix
def calc_u2_ycte(a, k, etacte, ui, uj, uk, ul):
    ik = ui, uk
    uijkl = k * (a/2) * integrate_ff(ik) * calc_f(uj, etacte) * calc_f(ul, etacte)
    return uijkl

