import numpy as np
import pandas as pd

#TODO:
# add more advanced strength criteria

__all__ = ['LaminateStrength']


class LaminateStrength():
    '''
    Creates a LaminateStrength object to evaluate strength.
    The class is capable of calculating strain and stress at the midplane, as well
    as ply-by-ply at the top and bottom, both in the laminate and material directions.

    Parameters
    ----------
    dproperty : LaminateProperty
        A laminate property object
    Nxx : float, int, optional, default 0
        Membrane load in x direction.
    Nyy : float, int, optional, default 0
        Membrane load in y direction.    
    Nxy : float, int, optional, default 0
        Membrane load in xy direction.
    Mxx : float, int, optional, default 0
        Moment in x direction.
    Myy : float, int, optional, default 0
        Moment in y direction.
    Mxy : float, int, optional, default 0
        Moment in xy direction.
    '''

    def __init__(
            self, 
            dproperty, 
            Nxx=0,
            Nyy=0,
            Nxy=0,
            Mxx=0,
            Myy=0,
            Mxy=0):
        self.dproperty = dproperty
        self.Nxx = Nxx
        self.Nyy = Nyy
        self.Nxy = Nxy
        self.Mxx = Mxx
        self.Myy = Myy
        self.Mxy = Mxy
        self._positions = None #used for plotting.


    def epsilon0(self): #ok
        '''
        Calculates the strains of the laminate
        
        Returns
        -------
        epsilon0 : numpy ndarray
            [epsilonx0, epsilony0, epsilonxy0, kappax0, kappay0, kappaxy0]
        '''

        ABD = self.dproperty.ABD
        abd = np.linalg.inv(ABD) #equivalent relations on pg 153 of Daniel
        N = np.array([self.Nxx, self.Nyy, self.Nxy, self.Mxx, self.Myy, self.Mxy])
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
        z = self.dproperty.z_position
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

    def _epsilonk_123(self):

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

    def _allowablek_123(self):
        '''
        Finds the allowable stresses for each lamina
        Returns
        -------
        allowablek_123: list
            [allowable1, allowable2, ..., allowablek], with allowable as np.ndarray()
            For each allowable we have
            allowablek = np.ndarray([T1, C1, T2, C2, S])k
        '''
        laminate = self.dproperty.layup

        allowablek_123 = []
        for layer in laminate:
            angle,ply = layer
            maxstresses = np.array([ply.t1,
                                    ply.c1,
                                    ply.t2,
                                    ply.c2,
                                    ply.s,])
            allowablek_123.append(maxstresses)

        return allowablek_123

    def calculate_strain(self):
        '''
        Calculates strain ply by at laminate direction and material direction.

        Returns
        -------
        strains : pd.Dataframe
            ply by ply strains in plate direction and material direction       

        Note
        ----
        The sequence of the DataFrame starts from the BOTTOM OF THE LAYUP to the TOP  OF THE LAYUP.
        When defining the laminate, the first element of the list corresponds to the bottom-most layer. This is especially important for non-symmetric laminates.
        '''
        epsilonk = self._epsilonk()
        epsilonk_123 = self._epsilonk_123()
        stacking = self.dproperty.stacking
        z = self.dproperty.z_position

        cur_ply = 1
        data = {}
        data['ply'] = []
        data['z'] = [] # so the user can check and plot ply position
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
        for i, epsk_eps123_stck in  enumerate(zip(epsilonk, epsilonk_123, stacking)):
            epsilon, epsilon123, theta = epsk_eps123_stck
            epsilontop, epsilonbot = epsilon
            epsilon123top, epsilon123bot = epsilon123

            data['ply'].append(cur_ply)
            data['ply'].append(cur_ply)
            data['z'].append(z[i])
            data['z'].append(z[i+1])
            data['position'].append('bot')
            data['position'].append('top')
            data['angle'].append(theta)
            data['angle'].append(theta)            
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
        Calculates stress ply by at laminate direction and material direction.

        Returns
        -------
        stress : pd.Dataframe
            ply by ply stress in plate direction and material direction       

        Note
        ----
        The sequence of the DataFrame starts from the BOTTOM OF THE LAYUP to the TOP  OF THE LAYUP.
        When defining the laminate, the first element of the list corresponds to the bottom-most layer. This is especially important for non-symmetric laminates.
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
            data['angle'].append(theta)
            data['angle'].append(theta)
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

    def _calculate_maxstressallowable(self):
        '''
        Calculates max stress allowable for each ply.

        Returns
        -------
        allowable : pd.Dataframe
            ply by ply max stress allowables     

        Note
        ----
        The sequence of the DataFrame starts from the BOTTOM OF THE LAYUP to the TOP  OF THE LAYUP.
        When defining the laminate, the first element of the list corresponds to the bottom-most layer. This is especially important for non-symmetric laminates.
        '''

        allowablek_123 = self._allowablek_123()
        stacking = self.dproperty.stacking

        cur_ply = 1
        data = {}
        data['ply'] = []
        data['position'] = []
        data['angle'] = []
        data['angle'] = []
        data['t1'] = []
        data['t1'] = []
        data['c1'] = []
        data['c1'] = []
        data['t2'] = []
        data['t2'] = []
        data['c2'] = []
        data['c2'] = []
        data['s'] = []
        data['s'] = []
        for allowable123, theta in zip(allowablek_123, stacking):
            data['ply'].append(cur_ply)
            data['ply'].append(cur_ply)
            data['position'].append('top')
            data['position'].append('bot')
            data['angle'].append(theta)
            data['angle'].append(theta)
            data['t1'].append(allowable123[0])
            data['t1'].append(allowable123[0])
            data['c1'].append(allowable123[1])
            data['c1'].append(allowable123[1])
            data['t2'].append(allowable123[2])
            data['t2'].append(allowable123[2])
            data['c2'].append(allowable123[3])
            data['c2'].append(allowable123[3])
            data['s'].append(allowable123[4])
            data['s'].append(allowable123[4])
            cur_ply += 1
        pd.set_option('display.precision', 2)
        return pd.DataFrame(data)

    def calculate_maxstressmargin(self):
        '''
        Calculates margin of safety for each ply according to the max stress criteria.

        Returns
        -------
        margin : pd.Dataframe
            ply by ply max stress margin of safety     

        Note
        ----
        The sequence of the DataFrame starts from the BOTTOM OF THE LAYUP to the TOP  OF THE LAYUP.
        When defining the laminate, the first element of the list corresponds to the bottom-most layer. This is especially important for non-symmetric laminates.
        '''
        df_allowable = self._calculate_maxstressallowable()
        df_stress = self.calculate_stress()

        margin_t1 = df_allowable['t1'] / df_stress['sigma1'].mask(df_stress['sigma1'] < 0, 0) - 1 # replace negative stress with zero to yield infinite margins in tension
        margin_c1 = df_allowable['c1'] / df_stress['sigma1'].mask(df_stress['sigma1'] > 0, 0).abs() - 1 # replace positive stress with zero to yield infinite margins in compression. Also flip stress sign since allowables are all positive
        margin_t2 = df_allowable['t2'] / df_stress['sigma2'].mask(df_stress['sigma2'] < 0, 0) - 1
        margin_c2 = df_allowable['c2'] / df_stress['sigma2'].mask(df_stress['sigma2'] > 0, 0).abs() - 1
        margin_s = df_allowable['s'] / df_stress['tau12'].abs() - 1

        df_margin = pd.DataFrame({})
        df_margin['ply'] = df_allowable['ply']
        df_margin['position'] = df_allowable['position']
        df_margin['angle'] = df_allowable['angle']
        df_margin['margin_t1'] = margin_t1
        df_margin['margin_c1'] = margin_c1
        df_margin['margin_t2'] = margin_t2
        df_margin['margin_c2'] = margin_c2
        df_margin['margin_s'] = margin_s
        df_margin['critical'] = df_margin[['margin_t1','margin_c1','margin_t2','margin_c2','margin_s']].idxmin(axis=1) # find the critical failure mode of each ply
        pd.set_option('display.precision', 2)

        return df_margin