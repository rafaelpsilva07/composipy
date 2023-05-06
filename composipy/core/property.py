import numbers

import numpy as np
import scipy as sp

from composipy.utils import ComposipyValidator
from composipy.core.material import Material, OrthotropicMaterial, IsotropicMaterial


class Property(ComposipyValidator):
    pass


class LaminateProperty(Property):
    '''
    This class creates laminate object. It needs Material objects (to define plies) and the angle information.
    Some formulation characteristics are:
    Laminate formulations ares used (see References)
    Main reference is the chapter 4 of reference 2.

    Parameters
    ----------
    stacking : list or dict
        To define a angle stacking sequence
            An iterable containing the angles (in degrees) of layup.
        To define a stack based on lamination parameters.
            {xiA: [xiA1, xiA2, xiA3, xiA4],
            xiB: [xiB1, xiB2, xiB3, xiB4],
            xiD: [xiD1, xiD2, xiD3, xiD4],
            T: thickness}
    plies : Material or list
        A single Material or a list of OrthotropicMaterial object


    Example
    -------
    >>> from composipy import OrthotropicMaterial, LaminateProperty
    >>> ply_1 = OrthotropicMaterial(129500, 9370, 0.38, 5240, 0.2)
    >>> stacking = [90, 0, 90]
    >>> laminate_1 = Laminate(stacking, ply_1)
    >>> laminate_1.ABD # retunrs an array containing bending stiffness matrix [D] of the laminate
    >>> laminate_1.xiA # lamination parameters of extension
    >>> laminate_1.xiA # lamination parameters of bending

    Example
    -------
    >>> stacking = {'xiD': [0., 0., -1., 0.], 'T': 1.0} #using lamination parameters to define D
    >>> laminate_2 = Laminate(stacking, ply_1)
    >>> laminate_2.ABD # retunrs an array containing bending stiffness matrix [D] of the laminate
    >>> laminate_2.xiA # lamination parameters of extension
    >>> laminate_2.xiA # lamination parameters of bending
    '''

    def __init__(self, stacking, plies):
        # Checking layup
        if not isinstance(stacking, dict): # implements angle stacking sequence
            if isinstance(plies, Material):
                n_plies = len(stacking)
                plies = [plies for s in range(n_plies)]
            elif len(plies) != len(stacking):
                raise ValueError('Number of plies and number of stacking must match')
            xiA = None
            xiB = None
            xiD = None
            total_thickness = None
            layup = list(zip(stacking, plies))
        else:
            try:
                xiA = stacking['xiA']
            except KeyError:
                xiA = None
            try:             
                xiB = stacking['xiB']
            except KeyError:
                xiB = None
            try: 
                xiD = stacking['xiD']
            except KeyError:
                KeyError('xiD must be a key')
            try:
                total_thickness = stacking['T']
            except KeyError:
                KeyError('T must be a key')
            layup = []

        self.stacking = stacking
        self.plies = plies
        self.layup = layup
        self._z_position = None
        self._Q_layup = None
        self._T_layup = None
        self._A = None
        self._B = None
        self._D = None
        self._ABD = None
        self._ABD_p = None
        self._xiA = xiA
        self._xiB = xiB
        self._xiD = xiD
        self._total_thickness = total_thickness

#Properties
    @property
    def z_position(self):
    
        total_thickness = 0
        for t in self.layup:
            total_thickness += t[1].thickness
        
        current_z = -total_thickness/2
        ply_position = [current_z]
        for t in self.layup:
            current_z += t[1].thickness
            ply_position.append(current_z)
        
        return ply_position
    
    @property
    def Q_layup(self):    
        if self._Q_layup is None:
            
            self._Q_layup = []
            for theta, ply in self.layup:
                c = np.cos(theta*np.pi/180)
                s = np.sin(theta*np.pi/180)

                T_real = np.array([
                    [c**2, s**2, 2*c*s],
                    [s**2, c**2, -2*c*s],
                    [-c*s, c*s, c**2-s**2]
                    ])
                T_engineering =  np.array([
                    [c**2, s**2, c*s],
                    [s**2, c**2, -c*s],
                    [-2*c*s, 2*c*s, c**2-s**2]
                    ])

                self._Q_layup.append(
                    (np.linalg.inv(T_real))
                    @ ply.Q_0 
                    @ T_engineering
                    )
        return self._Q_layup
     
    @property
    def T_layup(self):
        if self._T_layup is None:
            
            self._T_layup = []
            for theta in self.layup:
                c = np.cos(theta[0]*np.pi/180)
                s = np.sin(theta[0]*np.pi/180)

                T_real = np.array([
                    [c**2, s**2, 2*c*s],
                    [s**2, c**2, -2*c*s],
                    [-c*s, c*s, c**2-s**2]
                    ])

                T_engineering =  np.array([
                    [c**2, s**2, c*s],
                    [s**2, c**2, -c*s],
                    [-2*c*s, 2*c*s, c**2-s**2]
                    ])

                self._T_layup.append([T_real,T_engineering])
        return self._T_layup

    @property
    def A(self):
        '''[A] Matrix as numpy.ndarray '''

        if not self._xiA is None:
            raise NotImplementedError
        if self._A is None:
            self._A = np.zeros(9).reshape(3,3)

            for i in enumerate(self.Q_layup):
                zk1 = self.z_position[i[0]+1]
                zk0 = self.z_position[i[0]]
                self._A += (zk1-zk0) * i[1]
        return self._A
    
    @property
    def B(self):
        '''[B] Matrix as numpy.ndarray '''
        if not self._xiB is None:
            raise NotImplementedError
        if self._B is None:
            self._B = np.zeros(9).reshape(3,3)

            for i in enumerate(self.Q_layup):
                zk1 = self.z_position[i[0]+1]
                zk0 = self.z_position[i[0]]
                self._B += (1/2) * (zk1**2-zk0**2) * i[1]
        return self._B
    
    @property
    def D(self):
        '''[D] Matrix as numpy.ndarray '''

        if not self._xiD is None:

            U1, U2, U3, U4, U5 = self.plies.Invariants()
            xi1, xi2, xi3, xi4 = self._xiD
            T = self._total_thickness

            D11 = T**3*(U1 + U2*xi1 + U3*xi3)/12
            D12 = T**3*(-U3*xi3 + U4)/12
            D13 = T**3*(U2*xi2/2 + U3*xi4)/12
            D21 = T**3*(-U3*xi3 + U4)/12
            D22 = T**3*(U1 - U2*xi1 + U3*xi3)/12
            D23 = T**3*(U2*xi2/2 - U3*xi4)/12
            D31 = T**3*(U2*xi2/2 + U3*xi4)/12
            D32 = T**3*(U2*xi2/2 - U3*xi4)/12
            D33 = T**3*(-U3*xi3 + U5)/12

            self._D = np.array([[D11, D12, D13],
                                [D21, D22, D23],
                                [D31, D32, D33]])

        if self._D is None:
            self._D = np.zeros(9).reshape(3,3)

            for i in enumerate(self.Q_layup):
                zk1 = self.z_position[i[0]+1]
                zk0 = self.z_position[i[0]]
                self._D += (1/3) * (zk1**3-zk0**3) * i[1]
        return self._D

    @property
    def ABD(self):
        '''[ABD] Matrices as numpy.ndarray'''
        self._ABD = np.vstack([
              np.hstack([self.A, self.B]),
              np.hstack([self.B, self.D])
            ])
        return self._ABD
   
    @property
    def xiA(self):
        '''Lamination parameter xiA'''
        xiA = np.zeros(4)
        T = sum([ply.thickness for ply in self.plies])        
        for i, angle in enumerate(self.stacking):
            angle *= np.pi / 180
            zk1 = self.z_position[i+1]
            zk0 = self.z_position[i]

            xiA[0] += (zk1-zk0) * np.cos(2*angle)
            xiA[1] += (zk1-zk0) * np.sin(2*angle)
            xiA[2] += (zk1-zk0) * np.cos(4*angle)
            xiA[3] += (zk1-zk0) * np.sin(4*angle)                        
        
        self._xiA = xiA / T
        return self._xiA

    @property
    def xiD(self):
        '''Lamination parameter xiD'''
        xiD = np.zeros(4)
        T = sum([ply.thickness for ply in self.plies])        
        for i, angle in enumerate(self.stacking):
            angle *= np.pi / 180
            zk1 = self.z_position[i+1]
            zk0 = self.z_position[i]

            xiD[0] += (zk1**3-zk0**3) * np.cos(2*angle)
            xiD[1] += (zk1**3-zk0**3) * np.sin(2*angle)
            xiD[2] += (zk1**3-zk0**3) * np.cos(4*angle)
            xiD[3] += (zk1**3-zk0**3) * np.sin(4*angle)                        
        return 4 * xiD / T**3

       #Representation
    def __repr__(self):
                    
        representation = f'composipy.LaminateProperty\n'
        representation += f'stacking = {self.stacking}'
        return representation

#Comparisons
    def __eq__(self, other):
        if isinstance(other, LaminateProperty):
            return (self.layup == other.layup)
        return NotImplemented

