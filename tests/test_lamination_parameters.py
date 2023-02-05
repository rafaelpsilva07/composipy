import pytest
import numpy as np

from composipy import OrthotropicMaterial, LaminateProperty


@pytest.fixture
def stacking_plies():
    E11 = 120000
    E22 = 6500
    nu12 = 0.07
    G12 = 4550
    t = 0.21

    laminaprop = (E11, E22, nu12, G12, G12, G12)
    stack = [45, -45, 0, 90, 45, 0, 0]
    stack += stack[::-1]

    mat1 = OrthotropicMaterial(E11, E22, nu12, G12, t)
   
    return stack, mat1


def test_xiA(stacking_plies):
    stacking, plies = stacking_plies
    xiA_reference = np.array([
        0.2857142857142857,
        0.14285714285714285,
        0.14285714285714285,
        -1.7494954273533614e-17
    ]) # calculated with composites
    laminate = LaminateProperty(stacking, plies)
    np.testing.assert_almost_equal(laminate.xiA, xiA_reference)


def test_xiD(stacking_plies):
    stacking, plies = stacking_plies
    xiD_reference = np.array([
        0.09329446064139943,
        0.16034985422740528,
        -0.38192419825072876,
        -6.783757779533445e-18
    ]) # calculated with composites
    laminate = LaminateProperty(stacking, plies)
    np.testing.assert_almost_equal(laminate.xiD, xiD_reference)