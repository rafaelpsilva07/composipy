import numpy as np

from time import time
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from itertools import product

from composipy.core.property import Property
from composipy.utils import ComposipyValidator
from composipy.pre_integrated_component.build_k import *


class Structure(ComposipyValidator):
      pass


class PlateStructure(Structure):
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


    def _compute_constraints(self):
        if self.constraints == 'PINNED':
            x0 = ['TX', 'TY', 'TZ']
            xa = ['TX', 'TY', 'TZ']
            y0 = ['TX', 'TY', 'TZ']
            yb = ['TX', 'TY', 'TZ']
        elif self.constraints == 'CLAMPED':
            x0 = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
            xa = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
            y0 = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
            yb = ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
        else:
            x0 = self.constraints['x0']
            xa = self.constraints['xa']
            y0 = self.constraints['y0']
            yb = self.constraints['yb']

        sm = [i for i in range(self.m+4)]
        sn = [i for i in range(self.n+4)]

        um, un = sm.copy(), sn.copy()
        vm, vn = sm.copy(), sn.copy()
        wm, wn = sm.copy(), sn.copy()

        # x0
        if 'TX' in x0:
              um.remove(0)
        if 'TY' in x0:
              vm.remove(0)
        if 'TZ' in x0:
              wm.remove(0)
        if 'RX' in x0:
              um.remove(1)
        if 'RY' in x0:
              vm.remove(1)
        if 'RZ' in x0:
              wm.remove(1)
        #xa
        if 'TX' in xa:
              um.remove(2)
        if 'TY' in xa:
              vm.remove(2)
        if 'TZ' in xa:
              wm.remove(2)
        if 'RX' in xa:
              um.remove(3)
        if 'RY' in xa:
              vm.remove(3)
        if 'RZ' in xa:
              wm.remove(3)
        
        #y0
        if 'TX' in y0:
              un.remove(0)
        if 'TY' in y0:
              vn.remove(0)
        if 'TZ' in y0:
              wn.remove(0)
        if 'RX' in y0:
              un.remove(1)
        if 'RY' in y0:
              vn.remove(1)
        if 'RZ' in y0:
              wn.remove(1)
        
        #yb
        if 'TX' in yb:
              un.remove(2)
        if 'TY' in yb:
              vn.remove(2)
        if 'TZ' in yb:
              wn.remove(2)
        if 'RX' in yb:
              un.remove(3)
        if 'RY' in yb:
              vn.remove(3)
        if 'RZ' in yb:
              wn.remove(3)
        
        um, un = um[0:self.m], un[0:self.n]
        vm, vn = vm[0:self.m], un[0:self.n]
        wm, wn = wm[0:self.m], un[0:self.n]

        uidx = list(product(um, un, um, un))
        vidx = list(product(vm, vn, vm, vn))
        widx = list(product(wm, wn, wm, wn))

        return (uidx, vidx, widx)
    
    #TODO: parallel process, sparse matrix, cython
    def calc_K_KG(self):
        A11, A12, A16, B11, B12, B16 = self.dproperty.ABD[0, ::]
        A12, A22, A26, B12, B22, B26 = self.dproperty.ABD[1, ::]
        A16, A26, A66, B16, B26, B66 = self.dproperty.ABD[2, ::]
        B11, B12, B16, D11, D12, D16 = self.dproperty.ABD[3, ::]
        B12, B22, B26, D12, D22, D26 = self.dproperty.ABD[4, ::]
        B16, B26, B66, D16, D26, D66 = self.dproperty.ABD[5, ::]

        k11, k12, k13, k21, k22, k23, k31, k32, k33 = [], [], [], [], [], [], [], [], []
        k33g = []

        uidx, vidx, widx = self._compute_constraints()

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
        r"""Linear Buckling Analysis

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

        return eigvals, eigvecs