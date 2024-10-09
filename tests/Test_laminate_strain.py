'''
Test of strains of LaminateStrength
===================================

Test LaminateStrength class using examples from reference 1

References
----------
1 - https://ntrs.nasa.gov/api/citations/19950009349/downloads/19950009349.pdf
'''

import pytest
import numpy as np
import pandas as pd

from composipy import OrthotropicMaterial, LaminateProperty, LaminateStrength


@pytest.fixture
def laminateconfig():
    '''
    Material configuration from Example 1 in page 31.
    Laminate configuration from Example 4 in page 46
    This laminate is used in all tests
    '''
    E1 = 2.01e7 #lb/in2
    E2 = 1.301e6 #lb/in2
    v12 = 0.3
    G12 = 1.001e6 #lb/in2
    t = 0.005 #in
    material = OrthotropicMaterial(E1, E2, v12, G12, t)
    stacking = [0, 45, 45, 0]
    return LaminateProperty(stacking, material)


def test_example_3_midplane(laminateconfig):
    '''
    Example 3 in page 46
    This test check the stains at the midplane
    '''  
    Nxx = 500 #lb/in
    strain_analysis = LaminateStrength(dproperty=laminateconfig, Nxx=Nxx)
    calculated_strains = strain_analysis.epsilon0()
    reference_result = np.array([0.00221, -0.0007, -0.00115, 0, 0, 0])
    np.testing.assert_allclose(calculated_strains, reference_result, rtol=1e-2, atol=1e-4)