import numpy as np
from scipy.linalg import eig
import sys
from itertools import product

sys.path.append('D:/repositories/composipy/composipy')

from _descriptors import _NumberDescriptor

sys.path.append('D:/repositories/composipy/composipy/pre_integrated_component')
from build_k import *


class Panel:
    def __init__(
            self,
            property,
            a,
            b,
            boundary_conditions=None,
            Nxx=0,
            Nyy=0,
            Nxy=0,
            m=10,
            n=10,):
        #self.property = property
        # self.a = _NumberDescriptor(a, n_min=0, name='a')
        # self.b = _NumberDescriptor(b, n_min=0, name='a')
        # self.boundary_conditions = boundary_conditions
        # self.Nxx = _NumberDescriptor(Nxx, name='Nxx')
        # self.Nyy = _NumberDescriptor(Nyy, name='Nyy')
        # self.Nxy = _NumberDescriptor(Nxy, name='Nxy')
        # self.m = _NumberDescriptor(m, n_min=1, name='m')
        # self.n = _NumberDescriptor(n, n_min=1, name='n')

        self.property = property
        self.a = a
        self.b = b
        self.boundary_conditions = boundary_conditions
        self.Nxx = Nxx
        self.Nyy = Nyy
        self.Nxy = Nxy
        self.m = m
        self.n = n


    def _compute_boundary_conditions(self):
        #TODO: create this function
        # bardell_functions = {
        #     'TX': ,
        #     'TY': ,
        #     'TZ' ,
        #     'RX' ,
        #     'RX' ,
        #     'RZ' ,

        # }

        sm = [i for i in range(self.m+4)]
        sn = [i for i in range(self.n+4)]

        um, un = sm.copy(), sn.copy()
        vm, vn = sm.copy(), sn.copy()
        wm, wn = sm.copy(), sn.copy()

        um.remove(0); um.remove(1); um.remove(2); um.remove(3); un.remove(0); un.remove(1); un.remove(2); un.remove(3)
        vm.remove(0); vm.remove(1); vm.remove(2); vm.remove(3); vn.remove(0); vn.remove(1); vn.remove(2); vn.remove(3)
        wm.remove(0); wm.remove(1); wm.remove(2); wm.remove(3); wn.remove(0); wn.remove(1); wn.remove(2); wn.remove(3)

        um, un = um[0:self.m], un[0:self.n]
        vm, vn = vm[0:self.m], un[0:self.n]
        wm, wn = wm[0:self.m], un[0:self.n]
        #return {'um': um, 'un': un, 'vm': vm, 'vn': vn, 'wm': wm, 'wn': wn}

        uidx = list(product(um, um, um, um))
        vidx = list(product(vm, vm, vm, vm))
        widx = list(product(wm, wm, wm, wm))

        return (uidx, vidx, widx)
    
    
    def calc_K_KG(self):
        A11, A12, A16, B11, B12, B16 = self.property.ABD[0, ::]
        A12, A22, A26, B12, B22, B26 = self.property.ABD[1, ::]
        A16, A26, A66, B16, B26, B66 = self.property.ABD[2, ::]
        B11, B12, B16, D11, D12, D16 = self.property.ABD[3, ::]
        B12, B22, B26, D12, D22, D26 = self.property.ABD[4, ::]
        B16, B26, B66, D16, D26, D66 = self.property.ABD[5, ::]

        k11, k12, k13, k21, k22, k23, k31, k32, k33 = [], [], [], [], [], [], [], [], []
        k33g = []

        uidx, vidx, widx = self._compute_boundary_conditions()

        for i in range(self.m**2*self.n**2):
                ui, uj, uk, ul = uidx[i]
                vi, vj, vk, vl = vidx[i]
                wi, wj, wk, wl = widx[i]

                k11.append(calc_K11_ijkl(self.a, self.b, ui, uj, uk, ul, A11, A16, A66))
                k12.append(calc_k12_ijkl(self.a, self.b, ui, uj, vk, vl, A12, A16, A26, A66))
                k13.append(calc_k13_ijkl())
                k21.append(calc_k21_ijkl(self.a, self.b, vi, vj, uk, ul, A12, A16, A26, A66))
                k22.append(calc_k22_ijkl(self.a, self.b, vi, vj, vk, vl, A22, A26, A66))
                k23.append(calc_k23_ijkl())
                k31.append(calc_k31_ijkl())
                k32.append(calc_k32_ijkl())
                k33.append(calc_k33_ijkl(self.a, self.b, wi, wj, wk, wl, D11, D12, D22, D16, D26, D66))
                k33g.append(calc_kG33_ijkl(self.a, self.b, wi, wj, wk, wl, self.Nxx, self.Nyy, self.Nxy))
        
        k11 = np.array(k11).reshape(self.m**2, self.m**2)
        k12 = np.array(k12).reshape(self.m**2, self.m**2)
        k13 = np.array(k13).reshape(self.m**2, self.m**2)
        k21 = np.array(k21).reshape(self.m**2, self.m**2)
        k22 = np.array(k22).reshape(self.m**2, self.m**2)
        k23 = np.array(k23).reshape(self.m**2, self.m**2)
        k31 = np.array(k31).reshape(self.m**2, self.m**2)
        k32 = np.array(k32).reshape(self.m**2, self.m**2)
        k33 = np.array(k33).reshape(self.m**2, self.m**2)
        k00 = np.zeros(n**4).reshape(self.m**2, self.m**2)
        k33g = np.array(k33g).reshape(self.m**2, self.m**2)

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

        return (K, KG)
    

    def buckling_analysis(self):
        A, B = self.calc_K_KG()
        eig_values, eig_vectors = eig(A, B)
        eig_values = abs(eig_values)
        eig_values[np.isnan(eig_values)] = 99999999999

        return eig_values.min()
    
if __name__ == '__main__':

     sys.path.append('D:/repositories/composipy')

     from composipy import Ply, Laminate

     ply1 = Ply(129500, 9370, 0.38, 5240, 0.2)
     stacking = [-45, 45, 90, 0, 0, 0, 0, 90, 45, -45]
     laminate1 = Laminate(stacking, ply1)
     print(laminate1.ABD)

     panel = Panel(laminate1, 360, 360, Nxx=1)
     print(panel.buckling_analysis())