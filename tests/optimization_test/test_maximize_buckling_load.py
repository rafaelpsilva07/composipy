'''
Test of maximize_buckling_load function
========================================

Test using a plate from reference 1.

References
----------
1 - Gurdal, Z. and Haftka, R. T. “Optimization of Composite Laminates,”
 presented at the NATO Advanced Study Institute on Optimization of Large Structural Systems, Berchtesgaden, Germany, Sept. 23 – Oct. 4, 1991.
'''

import pytest
import numpy as np

from composipy.optimize import  maximize_buckling_load


@pytest.fixture
def panel_definition():
    panel = {'E1': 128e3,
             'E2': 13e3,
             'G12': 6.4e3,
             'v12': 0.3,
             'a': 508,
             'b': 254,
             'm': 7,
             'n': 7,
             'panel_constraint' : 'PINNED',
             'plot': False
    }
   
    return panel


def test_minimize_nyy_nxx_ratio0(panel_definition):
    T = 3.048
    Nxx = -1.
    Nyy = -0.
    Nxy = -0.

    res = maximize_buckling_load(T=T,Nxx=Nxx, Nyy=Nyy, Nxy=Nxy, **panel_definition)
    fun = res['fun']
    res = np.array(res['x'])

    fun_reference = np.array([-102.73893])
    res_reference = np.array([ 5.000e-05, -1.000e+00])

    #np.testing.assert_almost_equal(fun, fun_reference)
    assert np.allclose(res, res_reference, rtol=0.1, atol=0.1)
    assert np.isclose(fun, fun_reference, rtol=0.1, atol=0.1)


