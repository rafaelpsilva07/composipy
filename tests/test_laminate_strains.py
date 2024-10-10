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


@pytest.fixture
def example_6_setup(laminateconfig):
    '''
    Example 6 setup in page 56
    '''
    Mxx = 5 #lb/in
    return LaminateStrength(dproperty=laminateconfig, Mxx=Mxx)


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


def test_example_6_midplane(example_6_setup):
    '''This test check the stains at the midplane eq E6.1'''
    calculated_strains = example_6_setup.epsilon0()
    reference_result = np.array([0, 0, 0, 0.418, -0.165, -0.0975])
    np.testing.assert_allclose(calculated_strains, reference_result, rtol=1e-2, atol=1e-4)


def test_example_6_strain_top_ply1(example_6_setup):
    '''This test check the stains at the top of ply 1 in laminate direction'''
    df_strains = example_6_setup.calculate_strain()
    epsilonx = df_strains.iloc[0]['epsilonx']
    epsilony = df_strains.iloc[0]['epsilony']
    gammaxy = df_strains.iloc[0]['gammaxy']
    assert np.isclose(epsilonx, 0.00418, rtol=1e-3, atol=1e-4)
    assert np.isclose(epsilony, -0.00165, rtol=1e-3, atol=1e-4)
    assert np.isclose(gammaxy, -0.000975, rtol=1e-3, atol=1e-4)


def test_example_6_strain_top_ply2(example_6_setup):
    '''This test check the stains at the top of ply 1 in material direction'''
    df_strains = example_6_setup.calculate_strain()
    epsilon1 = df_strains.iloc[2]['epsilon1']
    epsilon2 = df_strains.iloc[2]['epsilon2']
    gamma12 = df_strains.iloc[2]['gamma12']
    assert np.isclose(epsilon1, 0.00039, rtol=1e-3, atol=1e-4)
    assert np.isclose(epsilon2, 0.00088, rtol=1e-3, atol=1e-4)
    assert np.isclose(gamma12, -0.00292, rtol=1e-3, atol=1e-4)


def test_example_6_stress_bot_ply3(example_6_setup):
    '''This test check the stresses for the ply 3 in material direction E6.12'''
    df_stresses = example_6_setup.calculate_stress()
    sigma1 = df_stresses.iloc[5]['sigma1']
    sigma2 = df_stresses.iloc[5]['sigma2']
    tau12 = df_stresses.iloc[5]['tau12']
    assert np.isclose(sigma1, -8094, rtol=1e-3, atol=100)
    assert np.isclose(sigma2, -1296, rtol=1e-3, atol=100)
    assert np.isclose(tau12, 2923, rtol=1e-3, atol=100)