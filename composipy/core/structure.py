import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
from time import time
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from itertools import product

from composipy.core.property import Property
from composipy.utils import ComposipyValidator
from composipy.pre_integrated_component.build_k import *
from composipy.pre_integrated_component.functions import sxieta
from composipy.core.structures import _generate_index

class Structure(ComposipyValidator):
      pass


class PlateStructure(Structure):
    '''
    This class defines a PlateStructure.

    Parameters
    ----------
    dproperty : composipy.Property
        Property of the plate.
    a : float
        Size of the plate parallel to the x axis.
    b : float
        Size of the plate parallel to the y axis.
    constraints : str or dict. Default: \"PINNED\"
        Plate boundary conditions.
    Nxx : float, default 0
        Linear force parallel x axis 
    Nyy : float, default 0
        Linear force parallel y axis
    Nxy : float, default 0
        Linear shear force
    m : int, default 10
        Size of shape function along x axis
    n : int, default 10
        Size of shape function along y axis


    Example
    -------
    >>> from composipy import OrthotropicMaterial, LaminateProperty, PlateStructure
    >>> ply_1 = OrthotropicMaterial(129500, 9370, 0.38, 5240, 0.2)
    >>> stacking = [90, 0, 90]
    >>> laminate = Laminate(stacking, ply_1)
    >>> constraints

    Note
    -----
    The ``constraint`` argument can be a str type \'PINNED\' or \'CLAMPED\'.
    Or a dictionary like described below:

    >>> constraints = {    
    ---     x0 = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
    ---     xa = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
    ---     y0 = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
    ---     yb = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
    --- }

    Attention, the rotations aren't considered around the axis. They are related to the shape function.
    That means, for example, the ``RZ`` turns around the x axis and moves the plate in z direction. 
    '''



    def __init__(
            self,
            dproperty,
            a,
            b,
            constraints={
                'x0': ['TX', 'TY', 'TZ'],
                'xa': ['TX', 'TY', 'TZ'],
                'y0': ['TX', 'TY', 'TZ'],
                'yb': ['TX', 'TY', 'TZ'],
            },
            Nxx=0,
            Nyy=0,
            Nxy=0,
            m=10,
            n=10):

        self.__dproperty = self._is_instance(dproperty, Property, 'dproperty')
        self.__a = self._float_number(a, n_min=0, name='a')
        self.__b = self._float_number(b, n_min=0, name='b')
        self.__constraints = constraints
        self.__Nxx = self._float_number(Nxx, name='Nxx')
        self.__Nyy = self._float_number(Nyy, name='Nyy')
        self.__Nxy = self._float_number(Nxy, name='Nxy')
        self.__m = self._int_number(m, n_min=1, name='m')
        self.__n = self._int_number(n, n_min=1, name='n')
        self.su_idx = None
        self.sv_idx = None
        self.sw_idx = None
        self.eigenvalue = None
        self.eigenvector = None

    @property
    def dproperty(self):
            return self.__dproperty
        
    @property
    def a(self):
        return self.__a
    
    @property
    def b(self):
            return self.__b

    @property
    def constraints(self):
            return self.__constraints
    
    @property
    def Nxx(self):
            return self.__Nxx
    @property
    def Nyy(self):
            return self.__Nyy
    
    @property
    def Nxy(self):
            return self.__Nxy
    
    @property
    def m(self):
            return self.__m
    
    @property
    def n(self):
            return self.__n

    
    #TODO: parallel process, sparse matrix, cython
    def calc_K_KG(self):
        '''
        Calculates the stiffness and geometrical stifness matrices.

        Returns
        -------
        K_KG : tuple
            A tuple of array containing (K, KG)
        '''

        A11, A12, A16, B11, B12, B16 = self.dproperty.ABD[0, ::]
        A12, A22, A26, B12, B22, B26 = self.dproperty.ABD[1, ::]
        A16, A26, A66, B16, B26, B66 = self.dproperty.ABD[2, ::]
        B11, B12, B16, D11, D12, D16 = self.dproperty.ABD[3, ::]
        B12, B22, B26, D12, D22, D26 = self.dproperty.ABD[4, ::]
        B16, B26, B66, D16, D26, D66 = self.dproperty.ABD[5, ::]

        k11, k12, k13, k21, k22, k23, k31, k32, k33 = [], [], [], [], [], [], [], [], []
        k33g = []

        uidx, vidx, widx = _generate_index(self.constraints, self.m, self.n)

        for i in range(self.m**2*self.n**2):
                ui, uj, uk, ul = uidx[i]
                vi, vj, vk, vl = vidx[i]
                wi, wj, wk, wl = widx[i]

                k11.append(calc_K11_ijkl(self.a, self.b, ui, uj, uk, ul, A11, A16, A66))
                k12.append(calc_k12_ijkl(self.a, self.b, ui, uj, vk, vl, A12, A16, A26, A66))
                k13.append(calc_k13_ijkl(self.a, self.b, ui, uj, wk, wl, B11, B12, B16, B26, B66))
                k21.append(calc_k21_ijkl(self.a, self.b, vi, vj, uk, ul, A12, A16, A26, A66))
                k22.append(calc_k22_ijkl(self.a, self.b, vi, vj, vk, vl, A22, A26, A66))
                k23.append(calc_k23_ijkl(self.a, self.b, vi, vj, wk, wl, B12, B16, B22, B26, B66))
                k31.append(calc_k31_ijkl(self.a, self.b, wi, wj, uk, ul, B11, B12, B16, B26, B66))
                k32.append(calc_k32_ijkl(self.a, self.b, wi, wj, vk, vl, B11, B12, B16, B22, B26, B66))
                k33.append(calc_k33_ijkl(self.a, self.b, wi, wj, wk, wl, D11, D12, D22, D16, D26, D66))
                k33g.append(calc_kG33_ijkl(self.a, self.b, wi, wj, wk, wl, self.Nxx, self.Nyy, self.Nxy))
        
        k11 = np.array(k11).reshape(self.m*self.n, self.m*self.n)
        k12 = np.array(k12).reshape(self.m*self.n, self.m*self.n)
        k13 = np.array(k13).reshape(self.m*self.n, self.m*self.n)
        k21 = np.array(k21).reshape(self.m*self.n, self.m*self.n)
        k22 = np.array(k22).reshape(self.m*self.n, self.m*self.n)
        k23 = np.array(k23).reshape(self.m*self.n, self.m*self.n)
        k31 = np.array(k31).reshape(self.m*self.n, self.m*self.n)
        k32 = np.array(k32).reshape(self.m*self.n, self.m*self.n)
        k33 = np.array(k33).reshape(self.m*self.n, self.m*self.n)
        k00 = np.zeros(self.m**2*self.n**2).reshape(self.m*self.n, self.m*self.n)
        k33g = np.array(k33g).reshape(self.m*self.n, self.m*self.n)

        K = np.vstack([
        np.hstack([k11, k12, k13]),
        np.hstack([k21, k22, k31]),
        np.hstack([k31, k32, k33])
        ])

        KG = np.vstack([
        np.hstack([k00, k00, k00]),
        np.hstack([k00, k00, k00]),
        np.hstack([k00, k00, k33g])
        ])

        return K, KG  

    def buckling_analysis(self, silent=True, num_eigvalues = 5):
        """
        Linear Buckling Analysis

        Parameters
        ----------
        silent : bool, optional
            Prints time to compute calculation
        num_eigvalues : int, optional
            Number of calculated eigenvalues.

        Returns
        -------
        eigvals, eigvecs
        """
        if not silent:
            ti = time()       
            print('calculating K and KG')
        K, KG = self.calc_K_KG()
        K, KG = csr_matrix(K), csr_matrix(KG)
        if not silent:
            print(f'K and KG calculated in {time()-ti} seconds')
        tol = 0
    
        if not silent:
            ti = time()
            print('Calculating eigenvalues')
        k = min(num_eigvalues, KG.shape[0]-2)
        mode = 'cayley'

        eigvals, eigvecs = eigsh(A=KG, k=k,
              which='SM', M=K, tol=tol, sigma=1., mode=mode)
        
        eigvals = -1./eigvals
    
        eigvals = eigvals
        eigvecs = eigvecs 
        if not silent:
             print(f'eigenvalues calculated in {time()-ti} seconds')

        self.eigenvalue, self.eigenvector = eigvals, eigvecs
        return eigvals, eigvecs
    
    def plot_eigenvalue(self, nth=0, ngridx=20, ngridy=20):
        if (not isinstance(self.eigenvalue, np.ndarray) 
                and not isinstance(self.eigenvector, np.ndarray)):
            self.buckling_analysis()
        
        c_values = self.eigenvector[:, nth] # ritz coefficients
        len_c_values = len(c_values)
        len_w = len(self.sw_idx)
        cw_values = c_values[len_c_values-len_w:len_c_values]
        xi_arr = np.linspace(-1, 1, ngridx)
        eta_arr = np.linspace(-1, 1, ngridy)
        xi_mesh, eta_mesh = np.meshgrid(xi_arr, eta_arr)
        z = np.zeros(ngridx*ngridy).reshape(ngridx, ngridy)

        for i in range(ngridx):
            for j in range(ngridy):
                sw = sxieta(self.sw_idx, xi_mesh[i, j], eta_mesh[i, j])
                wij = float((sw @ cw_values))
                z[i, j] = wij
        
        # coordinate transformation
        x_mesh = (self.a/2) * (xi_mesh+1)
        y_mesh = (self.b/2) * (eta_mesh+1)

        ax = plt.figure().add_subplot(projection='3d')
        surf = ax.plot_surface(x_mesh, y_mesh, z, cmap=cm.coolwarm)
        ax.set_xticks(np.linspace(0, max(self.a, self.b), 6))
        ax.set_yticks(np.linspace(0, max(self.a, self.b), 6))
        plt.show()

        return None