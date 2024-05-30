'''
Test of minimize_panel_weight function
========================================

Test using a plate from reference 1.

References
----------
1 - [1]	Liu, B. & Haftka, Raphael & Akgun, Mehmet. (2000).
Two-level composite wing structural optimization using response surfaces. Structural and Multidisciplinary Optimization. 20. 87-96. 10.1007/s001580050140.
'''

import pytest
import numpy as np

from composipy.optimize import  minimize_panel_weight


@pytest.fixture
def panel_definition():
    panel = {'E1': 127560,
             'E2': 13030,
             'G12': 6410,
             'v12': 0.3,
             'a': 1181.1,
             'b': 746.74,
             'm': 7,
             'n': 7,
             'panel_constraint' : 'PINNED',
             #'plot': False # removed from this function
    }
 
    return panel


def test_minimize_panel16_weight(panel_definition):
    Nxx, Nyy, Nxy = -3100.,	-907, 236

    res = minimize_panel_weight(Nxx=Nxx, Nyy=Nyy, Nxy=Nxy, **panel_definition)
    res = np.array(res['x'])

    res_reference = np.array([2.209e+01, -3.949e-01, -2.102e-01])

    assert np.allclose(res, res_reference, rtol=0.1, atol=0.1)

