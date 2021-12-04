# function
# This function calculates buckling load using a pre integration formulation

import sympy as sp
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

from ._plot_surface_function import _plot_surface
from sympy import sin, cos, pi, Piecewise, Ne
from ._pre_integrated_component import integration_KGx_component, integration_K_component

def buckling_load(a, b, D, n = 3, shape_plot = False , eig = False):
    '''
    ===========================================
    The function

    This function calculates the critical buckling load of a composite plate.
    Some formulation characteristics are:
        - The plate is considered simple supported
        - Distributed load entry in the X direction
        - Laminate must be symmetric
        - It uses Rayleigh Ritz formulation and shape function as series of sine waves
        - Pre integration of the KGx and K matrices was added, so the function is much much faster
        - Available after release 0.1.0
    ===========================================

    Use

    buckling_load(a, b, D, n = 3, shape_plot = False , eig = False)
    It returns the Nx load of the plate (float)

        - a ==> is plate lengh, paralell to the X direction and paralel to the load
        - b ==> is the plate width, paralel to the Y direction, where the load is distributed
        - D ==> bending matrix of the laminate
        - n ==> number of shape functions to be used (default = 3)
        - shape_plot ==> display the critical buckling shape of the plate (default = False)
        - eig ==> the function a dictionary with Nx, eigen values and eigen vectors  (default = False)


    ===========================================

    References:
        - JONES, M. Robert. Mechanics of Composite Materials. Taylor & Francis: 2nd ed 1999.
        - Analysis and Design of composite structures. Class notes. ITA 2020.
        - CASTRO, G. P. Saullo. Semi-Analytical Tools for the Analysis of Laminated Composite Cylindrical
    and Conical Imperfect Shells under Various Loading and Boundary Conditions. Clausthal University of Technology, 2014.

    '''
    #========================= KGX MATRIX CONSTRUCTION ============================
    # Symbolic values and variables
    i, j, k, l = sp.symbols(['i', 'j', 'k', 'l'], positive=True)
    x, y = sp.symbols(['x', 'y'])

    # Shape function components
    component1 = sp.sin((i)*sp.pi*x/a) * sp.sin((j)*sp.pi*y/b)
    component2 = sp.sin((k)*sp.pi*x/a) * sp.sin((l)*sp.pi*y/b)
    dx_component1 = sp.diff(component1, x)
    dx_component2 = sp.diff(component2, x)
    integrated_component = eval(integration_KGx_component)

    # Construction of KGx
    component_list = []
    KGx = sp.zeros(n**2)

    for i_index in range(n):
        for j_index in range(n):
            for k_index in range(n):
                for l_index in range(n):
                     cur_comp = integrated_component.subs([
                         (i, 1+i_index),
                         (j, 1+j_index),
                         (k, 1+k_index),
                         (l, 1+l_index)
                     ])
                     #subs(i,1+i_index).subs(j,1+j_index).subs(k,1+k_index).subs(l,1+l_index)
                     component_list.append(cur_comp)       

    KGx = np.array(component_list).reshape(n**2, n**2).astype(float)

    #ToDo
    #Built KGy and KGxy functions

    #========================= K MATRIX CONSTRUCTION ============================
    # Unpack D
    [D11, D12, D13],\
    [D21, D22, D23],\
    [D31, D32, D33] = D

    # Shape function components (used in Nww constructions)
    dxx_component1 = sp.diff(component1, x, x)
    dxx_component2 = sp.diff(component2, x, x)
    dyy_component1 = sp.diff(component1, y, y)
    dyy_component2 = sp.diff(component2, y, y)
    _2dxy_component1 = 2*sp.diff(component1, x, y)
    _2dxy_component2 = 2*sp.diff(component2, x, y)
    component_K = _2dxy_component2*(D33*_2dxy_component1 + D13*dxx_component1 + D23*dyy_component1)\
                  + dxx_component2*(D31*_2dxy_component1 + D11*dxx_component1 + D21*dyy_component1) \
                  + dyy_component2*(D32*_2dxy_component1 + D12*dxx_component1 + D22*dyy_component1)

    integrated_component = eval(integration_K_component)

    #Construction of K
    component_list_K = []
    K = sp.zeros(n**2)

    for i_index in range(n):
        for j_index in range(n):
            for k_index in range(n):
                for l_index in range(n):
                    cur_comp = integrated_component.subs([
                        (i, 1+i_index),
                        (j, 1+j_index),
                        (k, 1+k_index),
                        (l, 1+l_index)
                    ])
                    #(i,1+i_index).subs(j,1+j_index).subs(k,1+k_index).subs(l,1+l_index)
                    component_list_K.append(cur_comp)

    K = np.array(component_list_K).reshape(n**2, n**2).astype(float)

    #==================SOLVING THE EIGENPROBLEM AND FINDING THE RESULTS ============
    #eigenproblem
    eig_values, eig_vectors = linalg.eig(K,KGx)
    
    #finding the min eigenvalue
    min_eigen_value = min((eig_values.real**2)**(1/2))
    Nx = min_eigen_value * b
        
    #Finding the corresponding C vector
    i = np.where(eig_values == min_eigen_value)
    C_vector = eig_vectors[::, i[0]].real
    
    #complete shape function
    Nw = [[]]
    for i in range(n):
        for j in range(n):
            Nw[0].append(sp.sin((i+1)*sp.pi*x/a)*sp.sin((j+1)*sp.pi*y/b))
    Nw_values = sp.Matrix(Nw)

    #ploting function
    if shape_plot:
        _plot_surface(a, b, C_vector, Nw_values)
    
    #Builting the result dictionary
    result = {'Nx':Nx}
    if eig:
        result['eigen values'] = eig_values
        result['eigen vectors'] = eig_vectors


    return result