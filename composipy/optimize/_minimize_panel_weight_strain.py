import numpy as np

from scipy.optimize import NonlinearConstraint, minimize

from .utils import (
    epsilon_min_from_lp,
    penalty_g1,
    penalty_g2,
)


__all__ = ['minimize_panel_weight_strain']


def _objective_function(y0):
    T, *_ = y0  # (T, xi1, xi3)
    return T  # strain is governed by A ~ T (linear), so minimizing T directly is correct


def minimize_panel_weight_strain(
        a, b,
        E1, E2, v12, G12,
        Nxx=0, Nyy=0, Nxy=0,
        epsilon_allow=-3500e-6,
        x0=None,
        panel_constraint='PINNED',
        options=None,
        tol=None,
):
    '''
    Minimizes composite panel weight subject to a minimum principal midplane
    strain constraint.

    The optimization problem is:

        minimize    T
        subject to  epsilon_min(T, xi1, xi3) >= epsilon_allow   (strain)
                    xi3 + 2*xi1 + 1          >= 0              (LP feasibility g1)
                    xi3 - 2*xi1 + 1          >= 0              (LP feasibility g2)
                    T > 0,  -1 <= xi1 <= 1,  -1 <= xi3 <= 1

    where epsilon_min is the minimum (most compressive) principal strain at the
    laminate midplane, computed via Mohr\'s circle from the applied membrane loads.

    The objective is T (not T^3) because the membrane stiffness A scales linearly
    with thickness, making T the natural objective for strain-driven sizing.
    Compare with ``minimize_panel_weight`` which uses T^3 because bending
    stiffness D scales with T^3 for buckling-driven sizing.

    This function addresses the strain failure mode independently from buckling.
    To also check buckling, run ``minimize_panel_weight`` separately and compare
    the two results — the governing (thicker) solution drives the final design.

    Parameters
    ----------
    a : float
        Plate dimension along x axis (mm or consistent length unit).
    b : float
        Plate dimension along y axis.
    E1 : float
        Young modulus in the fibre direction.
    E2 : float
        Young modulus transverse to the fibre.
    v12 : float
        Poisson ratio.
    G12 : float
        In-plane shear modulus.
    Nxx : float, default 0
        Membrane load in x direction (N/mm or consistent force/length unit).
    Nyy : float, default 0
        Membrane load in y direction.
    Nxy : float, default 0
        Membrane shear load.
    epsilon_allow : float, default -3500e-6
        Minimum principal strain allowable (negative for compression).
        The constraint enforces epsilon_min >= epsilon_allow.
        Default is -3500 microstrain, a typical industry value for composite
        laminates under compression-dominated loading.
    x0 : list or None, default None
        Initial guess [T, xi1, xi3]. If None, defaults to [0.1, 0.0, 0.0].
    panel_constraint : str or dict, default \'PINNED\'
        Plate boundary conditions in composipy format. Kept for API consistency
        with other optimization functions; not used in the strain calculation
        (strains depend only on ABD and loads, not plate geometry or BCs).
    options : dict or None
        Options passed directly to ``scipy.optimize.minimize``.
    tol : float or None
        Tolerance passed directly to ``scipy.optimize.minimize``.

    Returns
    -------
    res : scipy.optimize.OptimizeResult
        Optimization result. Key fields:
        - ``res.x``       -> [T_opt, xi1_opt, xi3_opt]
        - ``res.fun``     -> optimized objective T
        - ``res.success`` -> True if converged

    Notes
    -----
    The minimum principal strain is computed at the laminate midplane using
    the membrane loads only (moments = 0). For bending-dominated cases a more
    detailed through-thickness analysis would be needed.

    The strain allowable ``epsilon_allow`` should be provided as a dimensionless
    strain value (e.g. -3500e-6, not -3500). Typical aerospace industry values
    for carbon/epoxy laminates range from -3000e-6 to -4000e-6 for compression.

    To obtain the final panel design, run both this function and
    ``minimize_panel_weight`` (buckling) and take the result with the larger
    thickness T — that is the governing constraint.

    Examples
    --------
    >>> from composipy.optimize import minimize_panel_weight_strain
    >>> res = minimize_panel_weight_strain(
    ...     a=500, b=250,
    ...     E1=128000, E2=13000, v12=0.3, G12=6400,
    ...     Nxx=-200, Nyy=0, Nxy=0,
    ...     epsilon_allow=-3500e-6,
    ... )
    >>> T_opt, xi1_opt, xi3_opt = res.x
    '''

    if x0 is None:
        x0 = [0.1, 0.0, 0.0]

    # ------------------------------------------------------------------ #
    # Constraint: epsilon_min >= epsilon_allow                            #
    # Physical loads used directly — epsilon_allow is in strain units.   #
    # ------------------------------------------------------------------ #
    def strain_constraint(y0):
        T, xi1, xi3 = y0
        eps_min = epsilon_min_from_lp(
            T, xi1, xi3,
            E1, E2, v12, G12,
            Nxx, Nyy, Nxy,
        )
        # eps_min - epsilon_allow >= 0: strain must not be more compressive
        # than the allowable (both negative for compression-dominated loads).
        return eps_min - epsilon_allow

    # ------------------------------------------------------------------ #
    # Assemble constraints and bounds                                     #
    # ------------------------------------------------------------------ #
    c_strain = NonlinearConstraint(strain_constraint, 0.0, np.inf)
    c_lp_g1  = NonlinearConstraint(penalty_g1, -0.0001, 10000)
    c_lp_g2  = NonlinearConstraint(penalty_g2, -0.0001, 10000)

    bounds = ([0.001, 1_000_000], [-1.0, 1.0], [-1.0, 1.0])

    res = minimize(
        _objective_function,
        x0,
        method='SLSQP',
        constraints=[c_strain, c_lp_g1, c_lp_g2],
        bounds=bounds,
        options=options,
        tol=tol,
    )

    return res
