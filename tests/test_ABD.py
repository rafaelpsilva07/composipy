'''
Test of LaminateProperty
========================

Test LaminateProperty class using example from page 201 of reference 1

References
----------
1 - Mendonça, Paulo de Tarso R. Materiais Compostos e Estruturas-sanduíche: Projeto e Análise. 2005 Editora Manole Ltda.
'''

import pytest
import numpy as np

from composipy import OrthotropicMaterial, LaminateProperty


@pytest.fixture
def stacking_plies():
    E1 = 54870
    E2 = 18320
    G12 = 8900
    v12 = 0.25
    ply1 = OrthotropicMaterial(E1, E2, v12, G12, 3.0)
    ply2 = OrthotropicMaterial(E1, E2, v12, G12, 6.0)
    stacking = [45, 0, 45]
    plies = [ply1, ply2, ply1]
   
    return stacking, plies


def test_A(stacking_plies):
    stacking, plies = stacking_plies
    A_reference = np.array([[515794.00930665, 100823.25453699, 55993.44346208],
                            [100823.25453699, 291820.23545835, 55993.44346208],
                            [ 55993.44346208, 55993.44346208, 151491.93392275]])
    laminate = LaminateProperty(stacking, plies)
    np.testing.assert_almost_equal(laminate.A, A_reference)


def test_B(stacking_plies):
    stacking, plies = stacking_plies
    B_reference = np.array([[0., 0., 0.],
                            [0., 0., 0.],
                            [0., 0., 0.]])
    laminate = LaminateProperty(stacking, plies)
    np.testing.assert_almost_equal(laminate.B, B_reference)


def test_D(stacking_plies):
    stacking, plies = stacking_plies
    D_reference = np.array([[4779418.7240577, 1612106.45974872, 1175862.31270358],
                            [1612106.45974872, 4107497.4025128, 1175862.31270358],
                            [1175862.31270358, 1175862.31270358, 2220130.61237785]])
    laminate = LaminateProperty(stacking, plies)
    np.testing.assert_almost_equal(laminate.D, D_reference)