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
    Laminate configuration from Example 3 in page 46
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


@pytest.fixture
def example_3_setup(laminateconfig):
    '''
    Example 3 setup in page 46
    '''
    Nxx = 500 #lb/in
    return LaminateStrength(dproperty=laminateconfig, Nxx=Nxx)


def test_example_3_midplane(example_3_setup):
    '''This test check the stains at the midplane'''
    calculated_strains = example_3_setup.epsilon0()
    reference_result = np.array([0.00221, -0.0007, -0.00115, 0, 0, 0])
    np.testing.assert_allclose(calculated_strains, reference_result, rtol=1e-2, atol=1e-4)


def test_example_3_strain_45(example_3_setup):
    '''This test check the stains for the 45° plies in material direction'''
    df_strains = example_3_setup.calculate_strain()
    epsilon1 = df_strains.iloc[2]['epsilon1']
    epsilon2 = df_strains.iloc[2]['epsilon2']
    gamma12 = df_strains.iloc[2]['gamma12']
    assert np.isclose(epsilon1, 0.00018, rtol=1e-3, atol=1e-4)
    assert np.isclose(epsilon2, 0.00133, rtol=1e-3, atol=1e-4)
    assert np.isclose(gamma12, -0.00291, rtol=1e-3, atol=1e-4)


def test_example_3_stress_45(example_3_setup):
    '''This test check the stresses for the 45° plies in material direction'''
    df_stresses = example_3_setup.calculate_stress()
    sigma1 = df_stresses.iloc[2]['sigma1']
    sigma2 = df_stresses.iloc[2]['sigma2']
    tau12 = df_stresses.iloc[2]['tau12']
    assert np.isclose(sigma1, 4146, rtol=1e-3, atol=100)
    assert np.isclose(sigma2, 1812, rtol=1e-3, atol=100)
    assert np.isclose(tau12, -2913, rtol=1e-3, atol=100)