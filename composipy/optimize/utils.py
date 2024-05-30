'''Utilities functions to be used in optimization'''

import warnings

import numpy as np
import matplotlib.pyplot as plt

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

    if len(y0) == 2: #unpacking for maximize_buckling
        xi1, xi3 = y0
    elif len(y0) == 3: #unpacking for minimize_panel_weight
        T, xi1, xi3 = y0

    g1 = xi3 + 2*xi1 + 1
    return g1


def penalty_g2(y0):
    '''Defines the boundary conditions g2 that delimits the triangle.'''
    if len(y0) == 2: #unpacking for maximize_buckling
        xi1, xi3 = y0
    elif len(y0) == 3: #unpacking for minimize_panel_weight
        T, xi1, xi3 = y0

    g2 = xi3 - 2*xi1 + 1
    return g2


def natural_constraint_g(y0):
    '''Defines the boundary conditions g1 that delimits the parabola.'''

    if len(y0) == 2: #unpacking for maximize_buckling
        xi1, xi3 = y0
    elif len(y0) == 3: #unpacking for minimize_panel_weight
        T, xi1, xi3 = y0

    g = xi3 - 2*xi1**2 + 1
    return g


def check_loads(Nxx, Nyy, Nxy):
    '''Given the load inputs, check if the function may present a problem'''

    # Positive loads warning   
    if Nxx >= 0 and Nyy >= 0 and Nxy == 0:
        warnings.warn(f'Buckling analysis is supposed to take at least a negative normal load or shear. (Nxx = {Nxx}, Nyy = {Nyy}, Nxy = {Nxy})')

    #Non normalized loads warning
    if abs(Nxx) > 1 or abs(Nyy) > 1 or abs(Nxy) > 1:
        warnings.warn(f'Loads will be normalized, prefer to use loads between -1 and 1 (Nxx = {Nxx}, Nyy = {Nyy}, Nxy = {Nxy})')


def _constraint(xi1, xi3, silent=True, penalty=True):
    g_natural = xi3 - 2*xi1**2 + 1
    g1 = xi3 + 2*xi1 + 1
    g2 = xi3 - 2*xi1 + 1

    if not penalty:
        g = g_natural
    elif xi1 < 0:
        g = g1
    else:
        g = g2

    if g < 0:
        if not silent:
            print(xi1, xi3)
        return False
    else:
        return True



def plot_optimization(a, b, T, m, n, E1, E2, v12,
                       G12, Nxx, Nyy, Nxy, panel_constraint, points_to_plot, res, penalty):

        print('generating plot...')
        x = np.linspace(-1., 1., points_to_plot)
        init_args = [0, 0, -1]
        Nx_arr = []
        xi1_arr, xi3_arr = [], []

        for xi1_ in x:
            for xi3_ in x:
                g_curr = _constraint(xi1_, xi3_, silent=True, penalty=penalty)
                if g_curr:
                    critc_N = Ncr_from_lp(a, b, T, m, n, xi1_, xi3_, E1, E2, v12, G12, Nxx, Nyy, Nxy, panel_constraint)
                    Nx_arr.append(critc_N)
                    xi1_arr.append(xi1_)
                    xi3_arr.append(xi3_)

                    if critc_N > init_args[2]:
                        init_args = [xi1_, xi3_, critc_N]
                else:
                    pass
                    Nx_arr.append(np.nan)
                    xi1_arr.append(xi1_)
                    xi3_arr.append(xi3_)

        g1_xi1 = np.linspace(-1, 0, points_to_plot)
        g2_xi1 = np.linspace(0, 1, points_to_plot)

        if penalty:
            g1_xi3 = -2*g1_xi1 - 1
            g2_xi3 = 2*g2_xi1 - 1
        else:
            print('penalty False')
            g1_xi3 = 2*g1_xi1**2 - 1
            g2_xi3 = 2*g2_xi1**2 - 1


        n_valid_points = int(np.sqrt(len(Nx_arr)))
        Nx_arr = np.array(Nx_arr)
        xi1_arr = np.array(xi1_arr)
        xi3_arr = np.array(xi3_arr)                       
        Nx_arr = Nx_arr.reshape(n_valid_points, n_valid_points)
        xi1_arr = xi1_arr.reshape(n_valid_points, n_valid_points)
        xi3_arr = xi3_arr.reshape(n_valid_points, n_valid_points)

        #plots
        plt.figure()

        # Nxcrit
        cs1 = plt.contour(xi1_arr, xi3_arr, Nx_arr, 20)
        plt.clabel(cs1)

        #g1 and g2
        plt.plot(g1_xi1, g1_xi3, 'k')
        plt.plot(g2_xi1, g2_xi3, 'k')

        #
        plt.plot(*res['x'], 'ro', label='optimum')


        plt.xlabel('xi1')
        plt.ylabel('xi3')
        plt.legend()
        plt.xlim([-1, 1])
        plt.ylim([-1, 1])
        plt.show()