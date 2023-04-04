from composipy.pre_integrated_component.functions import *

# TODO: espalhar para matriz
def calc_u2_xcte(b, k, xicte, ui, uj, uk, ul):
    jl = uj, ul
    uijkl = k * (b/2) * calc_f(ui, xicte) * calc_f(uk, xicte) * integrate_ff(jl)
    return uijkl


#TODO: fazer matriz simetrica como no theory
def calc_uaub_xcte(b, k, xicte, ui, uj, uk, ul):
    #jl = uj, ul
    #uijkl = k * (b/2) * calc_f(ui, xicte) * calc_f(uk, xicte) * integrate_ff(jl)
    #return uijkl
    pass


# TODO: espalhar para matrix
def calc_u2_ycte(a, k, etacte, ui, uj, uk, ul):
    ik = ui, uk
    uijkl = k * (a/2) * integrate_ff(ik) * calc_f(uj, etacte) * calc_f(ul, etacte)
    return uijkl

