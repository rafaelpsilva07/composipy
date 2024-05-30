import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import NonlinearConstraint, Bounds, minimize, LinearConstraint
from scipy.sparse.linalg import ArpackError

from .utils import Ncr_from_lp, normalize_critical_load, penalty_g1, penalty_g2, check_loads, natural_constraint_g, plot_optimization


__all__ = ['maximize_buckling_load']


def maximize_buckling_load(a, b, T,
                   E1, E2, v12, G12,
                   Nxx=0, Nyy=0, Nxy=0,
                   m=7, n=7,
                   panel_constraint='PINNED',
                   options=None,
                   tol=None,
                   plot=False,
                   points_to_plot=30,
                   penalty=True
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
    panel_constraint : composipy constraints format
        Plate boundary constraints
    options : scipyminimize optiions
    tol : scipy minimize tol
    plot : bool, default False
        Plot optimization function
    points_to_plot : float, default 30
      Number of points to plot
    penalty : bool, default True
      Applies boundary conditions into the feasible region.
      If penalty is True, the boundary condition is the triangle the delimits 0, 45, 90°
      If penalty is False, the boundary conditions is the parabola that defines the natural constraints.

    Returns
    ------- 
    res : scipy minimize result
    '''

    check_loads(Nxx, Nyy, Nxy)

    # Normalizing loads
    Nxx_norm, Nyy_norm, Nxy_norm, max_load = normalize_critical_load(Nxx, Nyy, Nxy)

    # Critical load error handling
    def critical_load(y0):
        xi1, xi3 = y0
        try:
            eig = Ncr_from_lp(a, b, T, m, n, xi1, xi3, E1, E2, v12, G12, Nxx_norm, Nyy_norm, Nxy_norm, constraints=panel_constraint)
        except ArpackError: #cause if xi1 and xi3 chosen by optimizer is out of the bounds, it may crash
            eig = -1
        return (-1) * eig #negative so the solver can maximize instead of minimize
    
    # Constraints and boundaries   
    b1 = ([-1.0, 1.0], [-1.0, 1.0])
    x0 = [0.0, 0.0]
    if penalty:
        c2 = NonlinearConstraint(penalty_g1, -0.0001, 10000)
        c3 = NonlinearConstraint(penalty_g2, -0.0001, 10000)
        res = minimize(critical_load, x0, method='SLSQP', constraints=[c2, c3], bounds=b1, options=options, tol=tol)
    elif not penalty:
        c = NonlinearConstraint(natural_constraint_g, -0.0001, 10000)
        res = minimize(critical_load, x0, method='SLSQP', constraints=[c], bounds=b1, options=options, tol=tol)
 

    if plot:
        plot_optimization(a, b, T, m, n, E1, E2, v12,
                       G12, Nxx, Nyy, Nxy, panel_constraint, points_to_plot, res, penalty)

    return res
