'''
Test of Q0 of OrthotropicMaterial
=================================

Test OrthotropicMaterial class using example from page 102 of reference 1

References
----------
1 - Mendonça, Paulo de Tarso R. Materiais Compostos e Estruturas-sanduíche: Projeto e Análise. 2005 Editora Manole Ltda.
'''

import pytest
import numpy as np

from composipy import OrthotropicMaterial


@pytest.fixture
def properties():
    E1 = 39000
    E2 = 8300
    G12 = 4000
    v12 = 0.26
    T = 1.0
    return tuple([E1, E2, v12, G12, T])


def test_orthotropic(properties):
    reference = np.array([[39569.26989624, 2189.49960093, 0.],
                          [2189.49960093, 8421.15231125, 0.],
                          [0., 0., 4000.]])
    ply1 = OrthotropicMaterial(*properties)
    np.testing.assert_array_almost_equal(ply1.Q_0, reference)