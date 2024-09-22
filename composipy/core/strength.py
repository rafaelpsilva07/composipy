import numpy as np


class LaminateStrength():
    '''
    Creates a LaminateStrength object to evaluate strength.

    Parameters
    ----------
    dproperty : LaminateProperty
        A laminate property object
    Nx : float, int, optional, default 0
        Membrane load in x direction.
    Ny : float, int, optional, default 0
        Membrane load in y direction.    
    Nxy : float, int, optional, default 0
        Membrane load in xy direction.
    Mx : float, int, optional, default 0
        Moment in x direction.
    My : float, int, optional, default 0
        Moment in y direction.
    Mxy : float, int, optional, default 0
        Moment in xy direction.
    '''


    def __init__(
            self, 
            dproperty, 
            Nx=0,
            Ny=0,
            Nxy=0,
            Mx=0,
            My=0,
            Mxy=0):
        self.dproperty = dproperty
        self.Nx = Nx
        self.Ny = Ny
        self.Nxy = Nxy
        self.Mx = Mx
        self.My = My
        self.Mxy = Mxy

