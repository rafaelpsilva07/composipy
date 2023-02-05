'''
Test of PlateStructure.buckling_analysis
========================================

Mathematical verification performed according compmech.
'''

import pytest
import numpy as np

from math import isclose

from composipy import OrthotropicMaterial, LaminateProperty, PlateStructure


@pytest.fixture
def stacking_plies():
    #Properties (N and mm)
    E1 = 60800
    E2 = 58250
    v12 = 0.07
    G12 = 4550
    t = 0.21

    # Plate Parameters
    a = 360
    b = 360
    m = 10
    n = 10

    stacking = [45,-45,0,90,45,-45,0,90,45,-45]
    stacking += stacking[::-1]

    a = 360.
    b = 360.

    mat1 = OrthotropicMaterial(E1, E2, v12, G12, t)
    laminate1 = LaminateProperty(stacking, mat1)

    plateconfig = {
        'dproperty' : laminate1,
        'a': a, 
        'b': b,
        'Nxx': -1, 
        'Nyy': -1,
        'm': m,
        'n': n
    }
   
    return plateconfig


def test_clamped_plate(plateconfig):
    first_eigen_reference = 119.76079572
    plate1 = PlateStructure(**plateconfig, constraints='Clamped')
    eigenvalues, eigvectors = plate1.buckling_analysis()
    assert isclose(eigenvalues[0], first_eigen_reference)