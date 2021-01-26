# function

import sympy as sp
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from ._plot_surface_function import _plot_surface


def critical_buckling(a, b, D, method = 'sin_series', n = 3, shape_plot = True):

    #Creating shape function
    Nw = [[]]
    x = sp.symbols('x')
    y = sp.symbols('y')
    

    if method == 'sin_series':
        for i in range(n):
            for j in range(n):
                Nw[0].append(sp.sin((i+1)*sp.pi*x/a)*sp.sin((j+1)*sp.pi*y/b))


    if method =='polynomial':
        for i in range(n):
            for j in range(n):
                Nw[0].append(x**(i+1)*y**(j+1)*(x-a)*(y-b))

        
    #Shape function differentiation  
    Nw_values = sp.Matrix(Nw)
    Nw_dx = sp.diff(Nw_values,x)
    Nw_dxx = sp.diff(Nw_values,x,x)
    Nw_dy = sp.diff(Nw_values,y)

    print('it may take a while...')
    
    #Matrices construction
    Ds = sp.Matrix(D)
    Nww_values = sp.Matrix([sp.diff(Nw_values,x,x),sp.diff(Nw_values,y,y),2*sp.diff(Nw_values,x,y)])
    print('Shape function ready!')
    KGx_values = np.array(sp.integrate(sp.integrate(sp.diff(Nw_values,x).T*sp.diff(Nw_values,x),(x,0,a)),(y,0,b)))
    print('[KGx] ready')
    #KGy_values = np.array(sp.integrate(sp.integrate(sp.diff(Nw_values,y).T*sp.diff(Nw_values,y),(x,0,a)),(y,0,b)))
    #KGxy_values = np.array(sp.integrate(sp.integrate(Nw_dx.T*Nw_dy + Nw_dy.T*Nw_dx,(x,0,a)),(y,0,b)))
    K_values = np.array(sp.integrate(sp.integrate(Nww_values.T*Ds*Nww_values,(x,0,a)),(y,0,b)))
    print('[K] ready')
    
    #eigenproblem
    A=np.array(K_values,dtype=float)
    B=np.array(KGx_values,dtype=float)
    eig_values, eig_vectors = linalg.eig(A,B)
    
    #finding the min eigenvalue
    min_eigen_value = min((eig_values.real**2)**(1/2))
    Nx = min_eigen_value*a
    print(f'buckling load is {Nx}')
    
    #Finding the corresponding C vector
    i = np.where(eig_values == min_eigen_value)
    C_vector = eig_vectors[::,i[0]].real
    
    #ploting function
    if shape_plot:
        _plot_surface(a, b, C_vector, Nw_values)
    
    return Nx