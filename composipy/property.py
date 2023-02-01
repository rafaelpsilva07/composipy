import numbers

import numpy as np
import scipy as sp

from composipy.validators import ComposipyValidator
from composipy.material import Material, OrthotropicMaterial, IsotropicMaterial


class Property(ComposipyValidator):
    pass


class LaminateProperty(Property):
    '''
    This class creates laminate object. It needs ply objects and the angle information.
    Some formulation characteristics are:
    Laminate formulations ares used (see References)
    Main reference is the chapter 4 of reference 2.

    Parameters
    ----------
    stacking : list
        An iterable containing the angles (in degrees) of layup.
    plies : Material or list
        A single Material or a list of OrthotropicMaterial object


    Example
    -------
    >>> from composipy import Ply, Laminate
    >>> ply_1 = OrthotropicMaterial(129500, 9370, 0.38, 5240, 0.2)
    >>> stacking = [90, 0, 90]
    >>> laminate = Laminate(stacking, ply_1)
    >>> laminate_1.D # retunrs an array containing bending stiffness matrix [D] of the laminate
    >>> laminate_1.A # retunrs an array containing stiffness matrix [A] of the laminate
    >>> laminate_1.B # retunrs an array containing coupled stiffness matrix [B] of the laminate
    >>> laminate_1.ABD() # returns an array containing ABD matrices.
    >>> laminate_1.ABD_p # returns an array containing the ABD prime matricies of the laminate

    
    References
    ----------
        1 - JONES, M. Robert. Mechanics of Composite Materials. Taylor & Francis: 2nd ed 1999.
        2 - Analysis and Design of composite structures. Class notes. ITA 2020.

    '''

    def __init__(self, stacking, plies):
        # Checking layup
        if isinstance(plies, Material):
            n_plies = len(stacking)
            plies = [plies for s in range(n_plies)]
        elif len(plies) != len(stacking):
            raise ValueError('Number of plies and number of stacking must match')

        self.stacking = stacking
        self.plies = plies
        self.layup = list(zip(stacking, plies))
        self._z_position = None
        self._Q_layup = None
        self._T_layup = None
        self._A = None
        self._B = None
        self._D = None
        self._ABD = None
        self._ABD_p = None
        self._xiA = None
        self._xiB = None

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
        return xiA / T

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

