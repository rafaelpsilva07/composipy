import itertools
import numbers
import numpy as np


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

    _ids = itertools.count(1) # counts number of instances

    def __init__(self, e1, e2, v12, g12, thickness,  t1, c1, t2, c2, s, name=None):
        if not isinstance(e1,numbers.Real) or e1<0:
            raise ValueError('e1 must be a positive number')    
        if not isinstance(e2,numbers.Real) or e2<0:
            raise ValueError('e2 must be a positive number')
        if not isinstance(v12,numbers.Real) or v12 <-0.5 or v12>0.5:
            raise ValueError('v12 must be a number between -0.5 and 0.5')
        if not isinstance(g12,numbers.Real) or g12<0:
            raise ValueError('g12 must be a positive number')
        if not isinstance(thickness,numbers.Real) or thickness<0:
            raise ValueError('thickness must be a positive number')
        if not isinstance(t1,numbers.Real) or t1<0:
            raise ValueError('t1 must be a positive number')
        if not isinstance(c1,numbers.Real) or c1>0:
            raise ValueError('c1 must be a negative number')
        if not isinstance(t2,numbers.Real) or t2<0:
            raise ValueError('t2 must be a positive number')
        if not isinstance(c2,numbers.Real) or c2>0:
            raise ValueError('c2 must be a negative number')
        if not isinstance(s,numbers.Real) or s<0:
            raise ValueError('c2 must be a positive number')

        self.__e1 = e1
        self.__e2 = e2
        self.__v12 = v12
        self.__g12 = g12
        self.__thickness = thickness
        self.__t1 = t1
        self.__c1 = c1
        self.__t2 = t2
        self.__c2 = c2
        self.__s = s
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

        if self.__Q_0 is None:
            Q11 = self.e1 / (1-self.v12*self.v21) 
            Q12 = (self.v12*self.e2) / (1-self.v12*self.v21)
            Q22 = self.e2 / (1-self.v12*self.v21)
            Q66 = self.g12
            self.__Q_0 = np.array([[Q11, Q12, 0],
                                   [Q12, Q22, 0],
                                   [0, 0, Q66]])
        return self.__Q_0

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

if __name__ == '__main__':
    
    E1 = 129500
    E2 = 9370
    v12 = 0.38
    G12 = 5240
    t = 0.2
    t1 = 500
    c1 = -200
    t2 = 50
    c2 = -250
    s = 25

    p1 = Ply(E1, E2, v12, G12, t, t1, c1, t2, c2, s)
    print(p1)
    print(p1.Q_0)