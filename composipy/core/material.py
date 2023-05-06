import itertools
import numpy as np

from composipy.utils import ComposipyValidator


__all__ = ['OrthotropicMaterial', 'IsotropicMaterial']


class Material(ComposipyValidator):
    pass


class OrthotropicMaterial(Material):
    """
    Creates an OrthotropicMaterial object.


    Parameters
    ----------
    e1 : float, int
        Young modulus in the 1st direction
    e2 : float, int
        Young modulus in the 2nd direction
    v12 : float
        poisson modulus 12
    g12 : float, int
        Shear modulus
    thickness : float, int
        thickness of the ply
    t1 : float, optional, default None
        tension allowable in direction 1.
    c1 : float, optional, default None
        compression allowable in direction 1.
    t2 : float, optional, default None
        tension allowable in direction 2.
    c2 : float, optional, default None
        compression allowable in direction 2.
    s : float, optional, default None
        shear allowable.
    name : str, optional, default None
        Name of the ply
 
    Examples
    --------
    >>> from composipy import OrthotropicMaterial
    >>> ply_1 = OrthotropicMaterial(129500, 9370, 0.38, 5240, 0.2)
    >>> ply_1.Q_0 # get the compliance matrix of the lamina
    Out:
    array([[130867.31382151,   3598.19426713,      0.     ],
        [  3598.19426713,   9468.93228191,      0.        ],
        [     0.        ,      0.        ,   5240.        ]])
 

    """

    _ids = itertools.count(1) # counts number of instances in order to generate names

    def __init__(
            self, 
            e1, 
            e2, 
            v12, 
            g12, 
            thickness, 
            t1=0, 
            c1=0, 
            t2=0, 
            c2=0, 
            s=0, 
            name=None):
        self.__e1 = self._float_number(e1, n_min=0, name='e1')
        self.__e2 = self._float_number(e2, n_min=0, name='')
        self.__v12 = self._float_number(v12, n_min=-0.5, n_max=0.5, name='v12')
        self.__g12 = self._float_number(g12, n_min=0, name='g12')
        self.__thickness = self._float_number(thickness, n_min=0, name='thickness')
        self.__t1 = self._float_number(t1, n_min=0, name='t1')
        self.__c1 = self._float_number(c1, n_min=0, name='c1')
        self.__t2 = self._float_number(t2, n_min=0, name='t2')
        self.__c2 = self._float_number(c2, n_min=0, name='c2')
        self.__s = self._float_number(s, n_min=0, name='s')
        self.__v21 = self.v12 * (self.e2/self.e1)
        self.__Q_0 = None #Calculated property
        self._id = next(self._ids)
        self.__name = name

    #Propertie getters    
    @property
    def e1(self):
        return self.__e1
    @property
    def e2(self):
        return self.__e2
    @property
    def v12(self):
        return self.__v12
    @property
    def g12(self):
        return self.__g12
    @property
    def thickness(self):
        return self.__thickness
    @property 
    def t1(self):
        return self.__t1
    @property 
    def c1(self):
        return self.__c1
    @property 
    def t2(self):
        return self.__t2
    @property 
    def c2(self):
        return self.__c2
    @property
    def s(self):
        return self.__s
    @property
    def v21(self):
        return self.__v21
    @property
    def name(self):
        if self.__name is None:
            self.__name = f'ply_{self._id}'
            return self.__name
        else:
            return self.__name

    @property
    def Q_0(self):
        """
        Get the compliance matrix of the material.       
        
        Returns
        -------
        self.Q_0 : numpy.ndarray
            Compliance matrix (3x3) of the Material 
        """

        if self.__Q_0 is None:
            Q11 = self.e1 / (1-self.v12*self.v21) 
            Q12 = (self.v12*self.e2) / (1-self.v12*self.v21)
            Q22 = self.e2 / (1-self.v12*self.v21)
            Q66 = self.g12
            self.__Q_0 = np.array([[Q11, Q12, 0],
                                   [Q12, Q22, 0],
                                   [0, 0, Q66]])
        return self.__Q_0

    def Invariants(self):
        '''
        Get the material invariants

        Returns
        -------
        Invariants : numpy.ndarray
           array([U1, U2, U3, U4, U5])
        
        '''

        Q11 = self.Q_0[0, 0]
        Q12 = self.Q_0[0, 1]
        Q22 = self.Q_0[1, 1]
        Q66 = self.Q_0[2, 2]

        U1 = 1/8 * (3*Q11 + 3*Q22 + 2*Q12 + 4*Q66)
        U2 = 1/2 * (Q11 - Q22)
        U3 = 1/8 * (Q11 + Q22 - 2*Q12 - 4*Q66)
        U4 = 1/8 * (Q11 + Q22 + 6*Q12 - 4*Q66)
        U5 = 1/8 * (Q11 + Q22 - 2*Q12 + 4*Q66)
        
        return np.array([U1, U2, U3, U4, U5])


    def __repr__(self):
        return (
            f'Ply({self.name}, E1 = {self.e1}, E2 = {self.e2}, \n\
            v12 = {self.v12}, G12 = {self.g12}, thickness = {self.thickness}, \n\
            T1 = {self.t1}, C1 = {self.c1}, T2 = {self.t2}, C2 = {self.c2}, S = {self.s})')
    
#Comparisons
    def __eq__(self, other):
        if isinstance(other, OrthotropicMaterial):
            return (self.e1 == other.e1
                    and self.e2 == other.e2
                    and self.v12 == other.v12
                    and self.g12 == other.g12
                    and self.thickness == other.thickness
                    and self.t1 == other.t1
                    and self.c1 == other.c1
                    and self.t2 == other.t2
                    and self.c2 == other.c2
                    and self.s == other.s)
        else:
            raise ValueError(f'Not a instance of Material')


class IsotropicMaterial(Material):
    """Creates an instance of Ply.
    Parameters
    ----------
    e : float
        Young modulus
    v : float
        poisson modulus 12
    g : float
        Shear modulus
    thickness : float
        thickness of the ply
    t : float, default None
        tension allowable.
    c : float, default None
        compression allowable.
    s : float, default None
        shear allowable.
    name : str, default None
        Name of the ply
    """

    _ids = itertools.count(1) # counts number of instances in order to generate names

    def __init__(
            self, 
            e, 
            v, 
            g, 
            thickness, 
            t=None, 
            c=None, 
            s=None, 
            name=None):
        self.__e = self._float_number(e, n_min=0, name='e')
        self.__v = self._float_number(v, n_min=-0.5, n_max=0.5, name='v')
        self.__g = self._float_number(g, n_min=0, name='g')
        self.__thickness = self._float_number(thickness, n_min=0, name='thickness')
        self.__t = self._float_number(t, n_min=0, name='t')
        self.__c = self._float_number(c, n_min=0, name='c')
        self.__s = self._float_number(s, n_min=0, name='s')
        self.__Q_0 = None #Calculated property
        self._id = next(self._ids)
        self.__name = name

    #Propertie getters    
    @property
    def e(self):
        return self.__e

    @property
    def v(self):
        return self.__v
    @property
    def g(self):
        return self.__g
    @property
    def thickness(self):
        return self.__thickness
    @property 
    def t(self):
        return self.__t
    @property 
    def c(self):
        return self.__c
    @property
    def s(self):
        return self.__s
    @property
    def v(self):
        return self.__v
    @property
    def name(self):
        if self.__name is None:
            self.__name = f'ply_{self._id}'
            return self.__name
        else:
            return self.__name

    @property
    def Q_0(self):
        """
        Get the compliance matrix of the material.       
        
        Returns
        -------
        self.Q_0 : numpy.ndarray
            Compliance matrix (3x3) of the Material 
        """

        if self.__Q_0 is None:
            Q11 = self.e / (1-self.v*self.v) 
            Q12 = (self.v*self.e) / (1-self.v*self.v)
            Q22 = self.e / (1-self.v*self.v)
            Q66 = self.g
            self.__Q_0 = np.array([[Q11, Q12, 0],
                                   [Q12, Q22, 0],
                                   [0, 0, Q66]])
        return self.__Q_0

    def __repr__(self):
        return (
            f'IsotropicMaterial({self.name}, E = {self.e}, \n\
            v = {self.v}, G = {self.g}, thickness = {self.thickness}, \n\
            T = {self.t}, C = {self.c}, S = {self.s})')

#Comparisons
    def __eq__(self, other):
        if isinstance(other, IsotropicMaterial):
            return (self.e == other.e
                    and self.v == other.v
                    and self.g == other.g
                    and self.thickness == other.thickness
                    and self.t == other.t
                    and self.c == other.c
                    and self.s == other.s)
        else:
            raise ValueError(f'Not a instance of Material')