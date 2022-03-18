import itertools
import numbers
import numpy as np


class Load:
    """Creates an instace of in-plane loads per unit of lenght.

    Parameters
    ----------
    Nx : float, int
        Normal load per unit lenght in x direction
    Ny : float, int
        Normal load per unit lenght in y direction
    Nxy : float
        Normal load per unit lenght in xy direction
    Mx : float, int
        Moment per unit lenght in x direction
    My : float, int
        Moment per unit lenght in y direction
    Mxy : float, int
        Moment per unit lenght in xy direction
   

    Examples
    --------
    >>> from composipy import Load
    >>> load_1 = Load(500, 250, 20, 1000, 400, 100)
    >>> ply_1.Q_0 # get the compliance matrix of the lamina
    Out:
    array([500, 250, 20, 1000, 400, 100])

    References
    ----------
    1. JONES, M. Robert. Mechanics of Composite Materials. Taylor & Francis: 2nd ed 1999.
    2. Analysis and Design of composite structures. Class notes. ITA 2020.
        
    """

    _ids = itertools.count(1) # counts number of instances

    def __init__(self, Nx, Ny, Nxy, Mx, My, Mxy):
        if not isinstance(Nx,numbers.Real):
            raise ValueError('Nx must be a real number')    
        if not isinstance(Ny,numbers.Real):
            raise ValueError('Ny must be a real number')
        if not isinstance(Nxy,numbers.Real):
            raise ValueError('Nxy must be a real number')
        if not isinstance(Mx,numbers.Real):
            raise ValueError('Mx must be a real number')
        if not isinstance(Mx,numbers.Real):
            raise ValueError('My must be a real number')
        if not isinstance(Mxy,numbers.Real):
            raise ValueError('Mxy must be a real number')

        self.__Nx = Nx
        self.__Ny = Ny
        self.__Nxy = Nxy
        self.__Mx = Mx
        self.__My = My
        self.__Mxy = Mxy
        self.__Loads = None #Calculated property
        self._id = next(self._ids)

    #Load getters    
    @property
    def Nx(self):
        return self.__Nx
    @property
    def Ny(self):
        return self.__Ny
    @property
    def Nxy(self):
        return self.__Nxy
    @property
    def Mx(self):
        return self.__Mx
    @property
    def My(self):
        return self.__My
    @property
    def Mxy(self):
        return self.__Mxy

    #Load array   
    @property
    def Load(self):
        """Create an array for the in-plane loads.

        Parameters
        ----------
        self : self
            Instance of Load
        
        Returns
        -------
        self.Load : numpy.ndarray
            In-plane load vector
        """

        if self._Load is None:
            self.__Load = np.array([self.Nx, self.Ny, self.Nxy, self.Mx, self.My, self.Mxy])
        return self.__Load



    def __repr__(self):
        return (
            f'Load(Nx = {self.Nx}, Ny = {self.Ny},  Nxy = {self.Nxy} \n\
                     Mx = {self.Mx}, My = {self.My}, Mxy = {self.Mxy})')
    
#Comparisons
    def __eq__(self, other):
        if isinstance(other, Load):
            return (self.Nx == other.Ny
                    and self.Ny == other.Ny
                    and self.Nxy == other.Nxy
                    and self.Mx == other.My
                    and self.My == other.My
                    and self.Mxy == other.Mxy)
        else:
            raise ValueError(f'Not a instance of Load')

if __name__ == '__main__':
    
    Nx = 500
    Ny = 250
    Nxy = 20
    Mx = 1000
    My = 400
    Mxy = 100

    load1 = Load(Nx, Ny, Nxy, Mx, My, Mxy)
    print(load1)
    print(load1.Load)