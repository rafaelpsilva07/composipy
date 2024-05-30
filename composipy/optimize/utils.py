'''Utilities functions to be used in optimization'''

from composipy import OrthotropicMaterial, LaminateProperty, PlateStructure



def Ncr_from_lp(a, b, T, m, n, xi1, xi3, E1, E2, v12, G12, Nxx, Nyy, Nxy, constraints):
    '''
    Calculates critical buckling load based in lamination parameters, geometry and constraints
    '''
    
    stacking = {
    'xiD': [xi1, 0.0, xi3, 0.0],
     'T': T}

    mat1 = OrthotropicMaterial(E1, E2, v12, G12, thickness=0.1)
    laminate1 = LaminateProperty(stacking, mat1)
    plate1 = PlateStructure(dproperty=laminate1, a=a, b=b, constraints=constraints, Nxx=Nxx, Nyy=Nyy, Nxy=Nxy, m=m, n=n)

    return plate1.buckling_analysis()[0][0]


def normalize_critical_load(Nxx, Nyy, Nxy):
    '''Normalizes the critical load to be used in optimization functions'''
    max_abs_load = [abs(Nxx), abs(Nyy), abs(Nxy)]
    max_abs_load.sort()
    max_abs_load = max_abs_load[-1]

    Nxx_norm = Nxx / max_abs_load
    Nyy_norm = Nyy / max_abs_load
    Nxy_norm = Nxy / max_abs_load

    return Nxx_norm, Nyy_norm, Nxy_norm, max_abs_load


def penalty_g1(y0):
    '''Defines the boundary conditions g1 that delimits the triangle.'''
    xi1, xi3 = y0
    g1 = xi3 + 2*xi1 + 1
    return g1


def penalty_g2(y0):
    '''Defines the boundary conditions g2 that delimits the triangle.'''
    xi1, xi3 = y0
    g2 = xi3 - 2*xi1 + 1
    return g2
