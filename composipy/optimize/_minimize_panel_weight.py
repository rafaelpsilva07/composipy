import numpy as np

from composipy import OrthotropicMaterial, LaminateProperty, PlateStructure
from scipy.optimize import NonlinearConstraint, Bounds, minimize, LinearConstraint
from scipy.sparse.linalg import ArpackError

from .utils import Ncr_from_lp, normalize_critical_load, penalty_g1, penalty_g2, check_loads


__all__ = ['minimize_panel_weight']


def _objective_function(y0):
    T, *_ = y0 #(T, xi1, xi3)
    return T**3


def minimize_panel_weight(a, b,
                   E1, E2, v12, G12,
                   Nxx=0, Nyy=0, Nxy=0,
                   m=7, n=7,
                   x0 = [0.1, 0.0, 0.0],
                   panel_constraint='PINNED',
                   options=None,
                   tol=None,
                   ):
    '''
    Parameters
    ----------
    a : float
        Plate dimension a
    b : float
        Plate dimension b
    E1 : float
        Young Modulus direction 1
    E2 : float
        Young Modulus direction 2
    v12 : float
        Poisson ratio
    G12 : float
        Shear Modulus
    Nxx : float, default 0
        Normal load in x direction
    Nyy : float, default 0
        Normal load in y direction
    Nxy : float, default 0
        Normal load in xy direction
    m : float, dafault 7
        Size of shape function along x axis
    n : float, default 0
        Size of shape function along y axis
    x0 : list, default [0.1, 0.0, 0.0]
        Initial value
    panel_constraint : composipy constraints format
        Plate boundary constraints
    options : scipyminimize optiions
    tol : scipy minimize tol
    plot : bool, default False
        Plot optimization function
    points_to_plot : float, default 30
      Number of points to plot

    Returns
    ------- 
    res : scipy minimize result
    '''

    check_loads(Nxx, Nyy, Nxy)

    # Normalizing loads
    Nxx_norm, Nyy_norm, Nxy_norm, max_load = normalize_critical_load(Nxx, Nyy, Nxy)

    # Critical load error handling
    def critical_load(y0):
        T, xi1, xi3 = y0
        try:
            eig = Ncr_from_lp(a, b, T, m, n, xi1, xi3, E1, E2, v12, G12, Nxx_norm, Nyy_norm, Nxy_norm, constraints=panel_constraint)
        except ArpackError: #cause if xi1 and xi3 chosen by optimizer is out of the bounds, it may crash
            eig = -1
        return eig
    
    # Constraints and boundaries
    c1 = NonlinearConstraint(critical_load, max_load, np.inf)
    c2 = NonlinearConstraint(penalty_g1, -0.0001, 10000)
    c3 = NonlinearConstraint(penalty_g2, -0.0001, 10000)
    b1 = ([0.001, 1000000], [-1.0, 1.0], [-1.0, 1.0])
   
    res = minimize(_objective_function, x0, method='SLSQP',
                   constraints=[c1, c2, c3], bounds=b1, options=options, tol=tol)

    return res
