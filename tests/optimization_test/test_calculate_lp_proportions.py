'''
Test of calc_lp_proportions function
========================================

Test if the function is calculating the correct proportions of 
0, 45, 90° plies using lamination parameters.
'''

import pytest
import numpy as np

from composipy.optimize import  calc_lp_proportions


def test_100_0_0():   
    '''
    0° : 100%
    90° : 0%
    45° : 0%   
    '''
    proportions = calc_lp_proportions(1, 1)
    p0 = proportions['p0']
    p90 = proportions['p90']
    p45 = proportions['p45']

    assert np.isclose(p0, 100.)
    assert np.isclose(p90, 0.)
    assert np.isclose(p45, 0.)


def test_0_100_0():   
    '''
    0° : 0%
    90° : 100%
    45° : 0%   
    '''
    proportions = calc_lp_proportions(-1, 1)
    p0 = proportions['p0']
    p90 = proportions['p90']
    p45 = proportions['p45']

    assert np.isclose(p0, 0.)
    assert np.isclose(p90, 100.)
    assert np.isclose(p45, 0.)


def test_0_0_100():   
    '''
    0° : 0%
    90° : 0%
    45° : 100%   
    '''
    proportions = calc_lp_proportions(0, -1)
    p0 = proportions['p0']
    p90 = proportions['p90']
    p45 = proportions['p45']

    assert np.isclose(p0, 0.)
    assert np.isclose(p90, 0.)
    assert np.isclose(p45, 100.)


def test_0_60_40():   
    '''
    0° : 0%
    90° : 60%
    45° : 40%   
    '''
    proportions = calc_lp_proportions(-0.936, 0.872)
    p0 = proportions['p0']
    p90 = proportions['p90']
    p45 = proportions['p45']

    assert np.isclose(p0, 0.)
    assert np.isclose(p90, 60.)
    assert np.isclose(p45, 40.)


def test_30_20_50():   
    '''
    0° : 30%
    90° : 20%
    45° : 50%   
    '''
    proportions = calc_lp_proportions(0.439, 0.75)
    p0 = proportions['p0']
    p90 = proportions['p90']
    p45 = proportions['p45']

    assert np.isclose(p0, 30.)
    assert np.isclose(p90, 20.)
    assert np.isclose(p45, 50.)


def test_50_0_50():   
    '''
    0° : 50%
    90° : 0%
    45° : 50%   
    '''
    proportions = calc_lp_proportions(0.875, 0.75)
    p0 = proportions['p0']
    p90 = proportions['p90']
    p45 = proportions['p45']

    assert np.isclose(p0, 50.)
    assert np.isclose(p90, 0.)
    assert np.isclose(p45, 50.)