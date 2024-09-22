import numpy as np


__all__ = ['OrthotropicMaterial', 'IsotropicMaterial']


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


    def epsilon0(self):
        '''
        Calculates the strains of the laminate
        Returns
        -------
        epsilon0: numpy ndarray
            [epsilonx0, epsilony0, epsilonxy0, kappax0, kappay0, kappaxy0]
        '''

        ABD = self.dproperty.ABD
        abd = np.linalg.inv(ABD) #equivalent relations on pg 153 of Daniel
        N = np.array([self.Nx, self.Ny, self.Nxy, self.Mx, self.My, self.Mxy])
        return abd @ N


    def epsilonk(self):
        '''
        Calculates the strains for each lamina
        Returns
        -------
        epsilonk: list
            [epsilon1, epsilon2, ..., epsilonk], with epsilon as np.ndarray()
            For each epsilon we have
            epsilonk = np.ndarray([epsilonx, epsilony, gammas])k, coordinate directions
            For reference refer to page 145 of Daniel equation 5.8
        '''       
        nplies = len(self.dproperty.stacking)
        z = self.dproperty.z_position
        zmid = [(z[i]+z[i+1]) / 2
               for i in range(nplies)]           
        epsilon0 = self.epsilon0()
        epsilon = np.array([epsilon0[0],
                            epsilon0[1],
                            epsilon0[2],])
        kappa = np.array([epsilon0[3],
                          epsilon0[4],
                          epsilon0[5],])        
        epsilonk = []
        for h in zmid:
            epsilonk.append(
                epsilon + h * kappa
            )
        return epsilonk


    def stressk(self):
        '''
        Calculates the strains for each lamina
        Returns
        -------
        epsilonk: list
            [epsilon1, epsilon2, ..., epsilonk], with epsilon as np.ndarray()
            For each epsilon we have
            epsilonk = np.ndarray([epsilonx, epsilony, gammas])k, coordinate directions
            For reference refer to page 145 of Daniel equation 5.8
        '''
        Qk = self.dproperty.Q_layup
        epsilonk = self.epsilonk()
        
        stressk = []
        for Q, epsilon in zip(Qk, epsilonk):
            stressk.append(
                Q @ epsilon
            )
        return stressk

    def epsilonk_123(self):
        ## TODO: compute strains in principal lamina directions
        ## See equation 3.58 to perform the transformation

        stacking = self.dproperty.stacking
        epsilonk = self.epsilonk()
        
        epsilonk_123 = []
        for theta, epsilon in zip(stacking, epsilonk):
            c = np.cos(theta*np.pi/180)
            s = np.sin(theta*np.pi/180)

            T = np.array([
                [c**2, s**2, 2*c*s],
                [s**2, c**2, -2*c*s],
                [-c*s, c*s, c**2-s**2]
                ])

            epsilonk_123.append(
                T @ epsilon
            )
        return epsilonk_123


    def strength(self):
        ## TODO: chose a criteria
        ## Study and apply criteria

        ## page 34 of Daniel provides material properties
        pass