import numbers

import numpy as np
import scipy as sp


import sys
sys.path.append('D:/repositories/composipy')

from composipy.ply_class import Ply

class Laminate:
    '''
    This class creates laminate object. It needs ply objects and the angle information.
    Some formulation characteristics are:
    Laminate formulations ares used (see References)
    Main reference is the chapter 4 of reference 2.

    Parameters
    ----------
    layup : list
        The layup instance is composed of a list containing a tuple to each ply.
        The tuple must contain the ply angle (float in degrees) with relation to the 1 direciton
        and a ply object (of Ply class).
       
        [(angle_of_ply_1, ply_1), (angle_of_ply_2, ply_2), ... (angle_of_ply_n, ply_n)]

    Returns
    -------
    None

    Example
    -------
    >>> from composipy import Ply, Laminate
    >>> ply_1 = Ply(129500, 9370, 0.38, 5240, 0.2)
    >>> layup_1 = [(90, ply_1), (0, ply_1), (90, ply_1)]
    >>> laminate = Laminate(layup_1)
    >>> laminate_1.D # retunrs a array containing bending stiffness matrix [D] of the laminate
    >>> laminate_1.A # retunrs a array containing stiffness matrix [A] of the laminate
    >>> laminate_1.B # retunrs a array containing coupled stiffness matrix [B] of the laminate
    >>> laminate_1.print_ABD() # method that prints ABD matrices of the laminate
    >>> laminate_1.ABD_p # returns an array containing the ABD prime matricies of the laminate

    
    References
    ----------
        1 - JONES, M. Robert. Mechanics of Composite Materials. Taylor & Francis: 2nd ed 1999.
        2 - Analysis and Design of composite structures. Class notes. ITA 2020.

    '''

    def __init__(self, layup):
  
        # Checking layup
        if not isinstance(layup, list):
            raise ValueError(
                'layup must be a list of tuples.\n\
                 Each tuple must contain a angle value and a Ply object'
                )       
        for ply in layup:
            if not isinstance(ply[0], numbers.Real):
                raise ValueError(f'the angle must be a real number. Check {ply}')
            if not isinstance(ply[1], Ply):
                raise ValueError(f'the ply mus be a Ply object. Check {ply}')

        self.layup = layup
        self._z_position = None
        self._Q_layup = None
        self._T_layup = None
        self._A = None
        self._B = None
        self._D = None
        self._ABD_p = None

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

                self._Q_layup.append(
                    (np.linalg.inv(T_real))
                    @ theta[1].Q_0 
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
    def ABD_p(self):
        ''' [A',B',D'], which is inverse of ABD Matrix, as numpy.ndarray '''
        if self._ABD_p is None:
            A_p = np.linalg.inv(self.A)
                  + (-np.linalg.inv(self.A) * self.B)
                  * (np.linalg.inv(self.D-self.B * np.linalg.inv(self.A) * self.B))
                  * (self.B * np.linalg.inv(self.A))
            B_p = (-np.linalg.inv(self.A) * self.B)
                  * np.linalg.inv(
                      self.D - self.B * np.linalg.inv(self.A) * self.B
                      )
            D_p = np.linalg.inv(
                self.D - self.B
                * np.linalg.inv(self.A)
                * self.B
                )
            ABD_p = sp.matrix(
                np.vstack(
                    (np.hstack((A_p,B_p)), np.hstack((B_p,D_p))))
                )
            ABD_p[np.isnan(ABD_p)] = 0
            self._ABD_p = ABD_p

        return self._ABD_p
    
    @staticmethod
    def _pprint(*args):
        nStrings = 20
        representation = ''

        for arg in args:
            nSpaces = nStrings - len(str(arg))
            if nSpaces < 0:
                nSpaces = 1
            representation += str(arg) + nSpaces * " "
        return representation

    #Representation (str not implemented)
    def __repr__(self):
        representation = ''
        for angle, ply in self.layup:
            if angle == abs(90):
                angle_repr = '|||||'
            elif angle == 0:
                angle_repr = '====='
            elif angle == 45 or angle == -45:
                angle_repr = '/////'
            else:
                angle_repr = '***'
            
            representation += self._pprint(ply.name, angle, angle_repr) + '\n'

        return representation

#Comparisons
    def __eq__(self, other):
        if isinstance(other, Laminate):
            return (self.layup == other.layup)
        return NotImplemented

#Methods
    def show_ABD(self):
        ''' This method prints ABD Matrix'''

        result = "[A]:\n" + str(self.A) + "\n"
        result += "[B]:\n" + str(self.B) + "\n"
        result += "[D]:\n" + str(self.D)
        print(result)
        return None
