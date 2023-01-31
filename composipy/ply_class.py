import itertools
import numbers
import numpy as np
import sys

sys.path.append('D:/repositories/composipy')

from composipy._descriptors import _NumberDescriptor


class Ply:
    """Creates an instance of Ply.

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
    t1 : float, default None
        tension allowable in direction 1.
    c1 : float, default None
        compression allowable in direction 1.
    t2 : float, default None
        tension allowable in direction 2.
    c2 : float, default None
        compression allowable in direction 2.
    s : float, default None
        shear allowable.
    name : str, default None
        Name of the ply

    Examples
    --------
    >>> from composipy import Ply
    >>> ply_1 = Ply(129500, 9370, 0.38, 5240, 0.2)
    >>> ply_1.Q_0 # get the compliance matrix of the lamina
    Out:
    array([[130867.31382151,   3598.19426713,      0.     ],
        [  3598.19426713,   9468.93228191,      0.        ],
        [     0.        ,      0.        ,   5240.        ]])

    References
    ----------
    1. JONES, M. Robert. Mechanics of Composite Materials. Taylor & Francis: 2nd ed 1999.
    2. Analysis and Design of composite structures. Class notes. ITA 2020.
        
    """

    _ids = itertools.count(1) # counts number of instances in order to generate names

    def __init__(
            self, 
            e1, 
            e2, 
            v12, 
            g12, 
            thickness, 
            t1=None, 
            c1=None, 
            t2=None, 
            c2=None, 
            s=None, 
            name=None):
        # self.e1 = _NumberDescriptor(e1, n_min=0, name='e1')
        # self.e2 = _NumberDescriptor(e2, n_min=0, name='e2')
        # self.v12 = _NumberDescriptor(v12, n_min=-0.5, n_max=0.5, name='v12')
        # self.g12 = _NumberDescriptor(g12, n_min=0, name='g12')
        # self.thickness = _NumberDescriptor(thickness, n_min=0, name='thickness')
        # self.t1 = _NumberDescriptor(t1, n_min=0, name='t1')
        # self.c1 = _NumberDescriptor(c1, n_min=0, name='c1')
        # self.t2 = _NumberDescriptor(t2, n_min=0, name='t2')
        # self.c2 = _NumberDescriptor(c2, n_min=0, name='c2')
        # self.s = _NumberDescriptor(s, n_min=0, name='s')
        self.e1 = e1
        self.e2 = e2
        self.v12 = v12
        self.g12 = g12
        self.thickness = thickness
        self.t1 = t1
        self.c1 = c1
        self.t2 = t2
        self.c2 = c2
        self.s = s
        self.v21 = self.v12 * (self.e2/self.e1)
        self._Q_0 = None #Calculated property
        self._id = next(self._ids)
        self.__name = name

    @property
    def name(self):
        if self.__name is None:
            self.__name = f'ply_{self._id}'
            return self.__name
        else:
            return self.__name
 
    #Calculated properties (have only getter)   
    @property
    def Q_0(self):
        """Get the compliance matrix of the instance.
        This is for a lamina under plane stress
        It uses the engineering constants of the instance.

        Parameters
        ----------
        self : self
            Instance of Ply
        
        Returns
        -------
        self.Q_0 : numpy.ndarray
            Compliance matrix of the ply
        """

        Q11 = self.e1 / (1-self.v12*self.v21) 
        Q12 = (self.v12*self.e2) / (1-self.v12*self.v21)
        Q22 = self.e2 / (1-self.v12*self.v21)
        Q66 = self.g12
        self._Q_0 = np.array([[Q11, Q12, 0],
                                [Q12, Q22, 0],
                                [0, 0, Q66]])
        return self._Q_0

    def __repr__(self):
        return (
            f'Ply({self.name}, E1 = {self.e1}, E2 = {self.e2}, \n\
            v12 = {self.v12}, G12 = {self.g12}, thickness = {self.thickness}, \n\
            T1 = {self.t1}, C1 = {self.c1}, T2 = {self.t2}, C2 = {self.c2}, S = {self.s})')
    
#Comparisons
    def __eq__(self, other):
        if isinstance(other, Ply):
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
            raise ValueError(f'Not a instance of Ply')
