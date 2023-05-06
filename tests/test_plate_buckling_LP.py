'''
Test of PlateStructure.buckling_analysis using lamination parameters
====================================================================




'''
import pytest
import numpy as np

from math import isclose

from composipy import OrthotropicMaterial, LaminateProperty, PlateStructure


@pytest.fixture
def matconfig():
    #Properties (N and mm)
    matconfig = {
    'e1' : 60000.,
    'e2': 7000., 
    'v12': 0.3,
    'g12': 4550, 
    'thickness': 0.01,
    }
   
    return matconfig


@pytest.fixture
def plateconfig():

    # Plate Parameters (N and mm)
    plateconfig = {
        'a': 200, 
        'b': 200,
        'Nxx': -1,
        'm': 10,
        'n': 10,
        'constraints': 'PINNED'
    }
   
    return plateconfig


def test_LP_plate(matconfig, plateconfig):
    '''
    Test plate build with lamination parameters
    xiD1 = 0 and xiD3 = -1, equivalent to +-45°
    '''
    stack = {'xiD': [0., 0., -1., 0.], 'T': 1.0}
    
    first_eigen_reference = 2.784#50183
    mat1 = OrthotropicMaterial(**matconfig)
    laminate1 = LaminateProperty(stacking=stack, plies=mat1)
    
    plate1 = PlateStructure(laminate1, **plateconfig)
    eigenvalues, eigvectors = plate1.buckling_analysis()
    assert isclose(eigenvalues[0], first_eigen_reference, abs_tol=1e-4, rel_tol=0.1)


def test_equivalent_plate(matconfig, plateconfig):
    '''
    Test plate build with lamination parameters   
    stack sequence of ([45, -45]25)s
    
    '''
    stack = 25 * [45, -45]
    stack += stack[::-1]

    
    first_eigen_reference = 2.784#376578
    mat1 = OrthotropicMaterial(**matconfig)
    laminate1 = LaminateProperty(stacking=stack, plies=mat1)
    
    plate1 = PlateStructure(laminate1, **plateconfig)
    eigenvalues, eigvectors = plate1.buckling_analysis()
    assert isclose(eigenvalues[0], first_eigen_reference, abs_tol=1e-4, rel_tol=0.1)



