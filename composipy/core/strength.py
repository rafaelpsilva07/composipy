import numpy as np
import pandas as pd

#TODO:
# tests
# improve documentation
# strength criteria

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


    def epsilon0(self): #ok
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


    def _epsilonk(self): #ok
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
        z = self.dproperty.z_position[::-1] #so bot is negative and top is positive
        zmid = [(z[i], z[i+1]) #tuple with top and bot
               for i in range(nplies)]
        epsilon0 = self.epsilon0()
        epsilon = np.array([epsilon0[0],
                            epsilon0[1],
                            epsilon0[2],])
        kappa = np.array([epsilon0[3],
                          epsilon0[4],
                          epsilon0[5],])        
        epsilonk = []
        for htop, hbot in zmid:
            epsilonk.append(
                (epsilon + htop * kappa,
                 epsilon + hbot * kappa)
            )
        return epsilonk


    def _stressk(self):
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
        epsilonk = self._epsilonk()
        
        stressk = []
        for Q, epsilon in zip(Qk, epsilonk):
            epsilontop, epsilonbot = epsilon
            stressk.append(
                (Q @ epsilontop,
                Q @ epsilonbot)
            )
        return stressk

    def epsilonk_123(self):

        stacking = self.dproperty.stacking
        epsilonk = self._epsilonk()
        
        epsilonk_123 = []
        for theta, epsilon in zip(stacking, epsilonk):
            epsilontop, epsilonbot = epsilon           
            c = np.cos(theta*np.pi/180)
            s = np.sin(theta*np.pi/180)
            epsilontop[2] /= 2 # tensorial shear strain (see nasa pg 50)
            epsilonbot[2] /= 2

            T = np.array([
                [c**2, s**2, 2*c*s],
                [s**2, c**2, -2*c*s],
                [-c*s, c*s, c**2-s**2]
                ])

            cur_epsilontop = T @ epsilontop
            cur_epsilonbot = T @ epsilonbot
            cur_epsilontop[2] *= 2
            cur_epsilonbot[2] *= 2
            epsilonk_123.append((cur_epsilontop, cur_epsilonbot)) # engineering shear strain (see nasa pg 50)
        
        return epsilonk_123


    def _stressk_123(self):

        stacking = self.dproperty.stacking
        stressk = self._stressk()
        
        stressk_123 = []
        for theta, stress in zip(stacking, stressk):
            stresstop, stressbot = stress   
            c = np.cos(theta*np.pi/180)
            s = np.sin(theta*np.pi/180)
            #stresstop[2] /= 2 # tensorial shear strain (see nasa pg 50)
            #stressbot[2] /=

            T = np.array([
                [c**2, s**2, 2*c*s],
                [s**2, c**2, -2*c*s],
                [-c*s, c*s, c**2-s**2]
                ])

            cur_stresstop = T @ stresstop
            cur_stressbot = T @ stressbot
            #cur_stresstop[2] *= 2
            #cur_stressbot[2] *= 2
            stressk_123.append((cur_stresstop, cur_stressbot)) 
        
        return stressk_123


    def calculate_strain(self):
        '''
        Returns
        -------
        strains : pd.Dataframe
            ply by ply strains in plate direction and material direction       
        '''
        epsilonk = self._epsilonk()
        epsilonk_123 = self.epsilonk_123()
        stacking = self.dproperty.stacking

        cur_ply = 1
        data = {}
        data['ply'] = []
        data['position'] = []
        data['angle'] = []       
        data['epsilonx'] = []
        data['epsilonx'] = []
        data['epsilony'] = []
        data['epsilony'] = []
        data['gammaxy'] = []
        data['gammaxy'] = []
        data['epsilon1'] = []
        data['epsilon1'] = []
        data['epsilon2'] = []
        data['epsilon2'] = []
        data['gamma12'] = []
        data['gamma12'] = []
        for epsilon, epsilon123, theta in zip(epsilonk, epsilonk_123, stacking):
            epsilontop, epsilonbot = epsilon
            epsilon123top, epsilon123bot = epsilon123

            data['ply'].append(cur_ply)
            data['ply'].append(cur_ply)
            data['position'].append('top')
            data['position'].append('bot')
            data['angle'] = theta
            data['angle'] = theta            
            data['epsilonx'].append(epsilontop[0]) #plate direction
            data['epsilonx'].append(epsilonbot[0])
            data['epsilony'].append(epsilontop[1])
            data['epsilony'].append(epsilonbot[1])
            data['gammaxy'].append(epsilontop[2])
            data['gammaxy'].append(epsilonbot[2])
            data['epsilon1'].append(epsilon123top[0]) #material direction
            data['epsilon1'].append(epsilon123bot[0])
            data['epsilon2'].append(epsilon123top[1])
            data['epsilon2'].append(epsilon123bot[1])
            data['gamma12'].append(epsilon123top[2])
            data['gamma12'].append(epsilon123bot[2])
            cur_ply += 1
        pd.set_option('display.precision', 2)
        return pd.DataFrame(data)


    def calculate_stress(self):
        '''
        Returns
        -------
        stress : pd.Dataframe
            ply by ply stress in plate direction and material direction       
        '''

        stressk = self._stressk()
        stressk_123 = self._stressk_123()
        stacking = self.dproperty.stacking

        cur_ply = 1
        data = {}
        data['ply'] = []
        data['position'] = []
        data['angle'] = []
        data['sigmax'] = []
        data['sigmax'] = []
        data['sigmay'] = []
        data['sigmay'] = []
        data['tauxy'] = []
        data['tauxy'] = []
        data['sigma1'] = []
        data['sigma1'] = []
        data['sigma2'] = []
        data['sigma2'] = []
        data['tau12'] = []
        data['tau12'] = []
        for sigma, sigma123, theta in zip(stressk, stressk_123, stacking):
            sigmatop, sigmabot = sigma
            sigma123top, sigma123bot = sigma123

            data['ply'].append(cur_ply)
            data['ply'].append(cur_ply)
            data['position'].append('top')
            data['position'].append('bot')
            data['angle'] = theta
            data['angle'] = theta
            data['sigmax'].append(sigmatop[0]) #plate direction
            data['sigmax'].append(sigmabot[0])
            data['sigmay'].append(sigmatop[1])
            data['sigmay'].append(sigmabot[1])
            data['tauxy'].append(sigmatop[2])
            data['tauxy'].append(sigmabot[2])
            data['sigma1'].append(sigma123top[0]) #material direction
            data['sigma1'].append(sigma123bot[0])
            data['sigma2'].append(sigma123top[1])
            data['sigma2'].append(sigma123bot[1])
            data['tau12'].append(sigma123top[2])
            data['tau12'].append(sigma123bot[2])
            cur_ply += 1
        pd.set_option('display.precision', 2)
        return pd.DataFrame(data)
