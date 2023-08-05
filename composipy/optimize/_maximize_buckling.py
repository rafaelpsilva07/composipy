import numpy as np
import matplotlib.pyplot as plt

from composipy import OrthotropicMaterial, LaminateProperty, PlateStructure
from scipy.optimize import NonlinearConstraint, Bounds, minimize, LinearConstraint
from scipy.sparse.linalg import ArpackError



def _Ncr(a, b, T, m, n, xi1, xi3, E1, E2, v12, G12, Nxx, Nyy, Nxy, constraints):
    '''
    Calculates critical buckling load based in lamination parameters   
    '''
    
    stacking = {
    'xiD': [xi1, 0.0, xi3, 0.0],
     'T': T}

    mat1 = OrthotropicMaterial(E1, E2, v12, G12, thickness=0.1)
    laminate1 = LaminateProperty(stacking, mat1)
    plate1 = PlateStructure(dproperty=laminate1, a=a, b=b, constraints=constraints, Nxx=Nxx, Nyy=Nyy, Nxy=Nxy, m=m, n=n)

    return plate1.buckling_analysis()[0][0]


def _penalty_g1(y0):
    xi1, xi3 = y0
    g1 = xi3 + 2*xi1 + 1
    return g1


def _penalty_g2(y0):
    xi1, xi3 = y0
    g2 = xi3 - 2*xi1 + 1
    return g2





def _define_critical_load(Nxx, Nyy, Nxy):
    max_abs_load = [abs(Nxx), abs(Nyy), abs(Nxy)]
    max_abs_load.sort()
    max_abs_load = max_abs_load[-1]

    Nxx_norm = Nxx / max_abs_load
    Nyy_norm = Nyy / max_abs_load
    Nxy_norm = Nxy / max_abs_load

    return Nxx_norm, Nyy_norm, Nxy_norm, max_abs_load


def maximize_buckling_load(a, b, T,
                   E1, E2, v12, G12,
                   Nxx, Nyy, Nxy,
                   m=7, n=7,
                   panel_constraint='PINNED',
                   options=None,
                   tol=None,
                   plot=False,
                   points_to_plot=30
                   ):
    '''
    Parameters
    ----------

    Returns
    -------  
    '''

    # Normalizing loads
    Nxx_norm, Nyy_norm, Nxy_norm, max_load = _define_critical_load(Nxx, Nyy, Nxy)

    # Critical load error handling
    def critical_load(y0):
        xi1, xi3 = y0
        try:
            eig = _Ncr(a, b, T, m, n, xi1, xi3, E1, E2, v12, G12, Nxx_norm, Nyy_norm, Nxy_norm, constraints=panel_constraint)
        except ArpackError: #cause if xi1 and xi3 chosen by optimizer is out of the bounds, it may crash
            eig = -1
        return (-1) * eig #negative so the solver can maximize instead of minimize
    
    # Constraints and boundaries
    c2 = NonlinearConstraint(_penalty_g1, -0.0001, 10000)
    c3 = NonlinearConstraint(_penalty_g2, -0.0001, 10000)
    b1 = ([-1.0, 1.0], [-1.0, 1.0])

    x0 = [0.0, 0.0]
   
    res = minimize(critical_load, x0, method='SLSQP', constraints=[c2, c3], bounds=b1, options=options, tol=tol)

    if plot:
        print('generating plot...')
        x = np.linspace(-1., 1., points_to_plot)
        init_args = [0, 0, -1]
        Nx_arr = []
        xi1_arr, xi3_arr = [], []

        def constraint(xi1, xi3, silent=True):
            #g = xi3 - 2*xi1**2 + 1
            g1 = xi3 + 2*xi1 + 1
            g2 = xi3 - 2*xi1 + 1

            if xi1 < 0:
                g = g1
            else:
                g = g2

            if g < 0:
                if not silent:
                    print(xi1, xi3)
                return False
            else:
                return True


        for xi1_ in x:
            for xi3_ in x:
                g_curr = constraint(xi1_, xi3_, silent=True)
                if g_curr:
                    critc_N = _Ncr(a, b, T, m, n, xi1_, xi3_, E1, E2, v12, G12, Nxx, Nyy, Nxy, panel_constraint)
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
        g1_xi3 = -2*g1_xi1 - 1

        g2_xi1 = np.linspace(0, 1, points_to_plot)
        g2_xi3 = 2*g2_xi1 - 1

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
        plt.show()


    return res



if __name__ == '__main__':

    E1 = 128e3
    E2 = 13e3
    G12 = 6.4e3
    nu12 = 0.3
    t = 0.127 # not used

    a = 508
    b = 254
    T = 3.048

    Nxx = -1
    Nyy = -.5
    Nxy = -0.

    m=7
    n = 7


    res = maximize_buckling_load(a, b, T, E1, E2, nu12, G12, Nxx, Nyy, Nxy, m, n, panel_constraint="PINNED", plot=False)

    print(res)

    print(res['x'])

    print(
        _Ncr(a, b, T, m, n, -0.235, -0.529, E1, E2, nu12, G12, Nxx, Nyy, Nxy, constraints="PINNED")
        )