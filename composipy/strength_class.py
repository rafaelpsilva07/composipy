
import numpy as np
import scipy as sp

from composipy.laminate_class import Laminate
from composipy.load_class import Load 

class Strength:
    '''
    This class creates the Strength object. It's inputs are Laminate and Load objects.
    The strength Object returns laminate strenght parameteres, including strain and stress in laminate and load coordinates. 
    It also returns the Tsai-Wu failure index for all plies in the current load conditions and the First Ply Failure
    analysis (based on Tsai-Wu failure criterion).

    Parameters
    ----------
    laminate : Laminate
        The laminate is an object that contains the layup scheme with its mechanical proprieties.      
        [(angle_of_ply_1, ply_1), (angle_of_ply_2, ply_2), ... (angle_of_ply_n, ply_n)]

    load: Load
        Load is and object that contains the in-plane forces and moments in x and y directions.
        Load(Nx, Ny, Nxy, Mx, My, Mxy)

    Example
    -------
    >>> from composipy import Ply, Laminate, Load, Strength
    >>> ply_1 = Ply(129500, 9370, 0.38, 5240, 0.2)
    >>> layup_1 = [(90, ply_1), (0, ply_1), (90, ply_1)]
    >>> laminate_1 = Laminate(layup_1)
    >>> load_1 = Load(100, 200, 50, 1000, 250, 100)
    >>> str_analysis = Strength(laminate_1, load_1)
    >>> str_analysis.mid_strain_xy # Returns the mid-plane strain and slope of the laminate in xy coordinates
    >>> str_analysis.mid_strain_12 # Returns the mid-plane strain and slope of the laminate in 12 coordinates
    >>> str_analysis.stress_xy # Returns the maximum stress for all plies in the xy coordinates
    >>> str_analysis.stress_12 # Returns the maximum stress for all plies in the 12 coordinates
    >>> str_analysis.TW_i # Returns the Tsai-Wu failure index for all plies 
    >>> str_anaysis.FPF # Prints the First Ply Failure load, failed plies, TW_i for the laminate and Strength Ratio
    
    References
    ----------
        1 - JONES, M. Robert. Mechanics of Composite Materials. Taylor & Francis: 2nd ed 1999.
        2 - Analysis and Design of composite structures. Class notes. ITA 2020.

    '''

    def __init__(self, laminate, load):
  
        # Checking layup
        if not isinstance(laminate, Laminate):
            raise ValueError(
                'laminate must be a Laminate object. Check {laminate}'
                )       
        if not isinstance(load, Load):
            raise ValueError('load must be a Load object. Check {load}')

        self.laminate = laminate
        self.load = load
        self._mid_strain_xy = None
        self._mid_strain_12 = None
        self._stress_12 = None
        self._stress_xy = None
        self._TW_i = None
        self._FPF = None

#Properties
# TODO we can create methods instead properties. Like get_strain_xy and so on.
# These methods might be documented since they are accessible by the user.
    @property 
    def mid_strain_xy(self):
        ''' [strain_x, strain_y, strain_xy, slope_x, slope_y, slope_xy] '''
        if self._mid_strain_xy is None:
            load_array = sp.matrix(
                np.vstack((
                    self.load.Nx,
                    self.load.Ny,
                    self.load.Nxy,
                    self.load.Mx,
                    self.load.My,
                    self.load.Mxy)))
            self._mid_strain_xy = self.laminate.ABD_p * load_array
        return self._mid_strain_xy

    @property 
    def mid_strain_12(self):
        ''' [strain_1, strain_2, strain_12, slope_1, slope_2, slope_12]'''
        if self._mid_strain_12 is None:
            mid_strain_12_lam = []
            for T_ply in self.laminate.T_layup:
                mid_strain_12_ply = T_ply[1] * self.mid_strain_xy[0:3]
                mid_strain_12_lam.append(mid_strain_12_ply)
            self._mid_strain_12 = mid_strain_12_lam
        return self._mid_strain_12

    @property
    def stress_xy(self):
        ''' [sigma_x, sigma_y, tau_xy] for each ply'''
        if self._stress_xy is None:    
            ply_num = 0
            stress_xy_lam = []
            for q_bar in self.laminate.Q_layup:
                stress_xy_ply = q_bar
                                * self.mid_strain_xy[0:3] 
                                + self.laminate.z_position[ply_num]
                                * self.mid_strain_xy[3:7]
                stress_xy_lam.append(stress_xy_ply)
                ply_num += 1
            self._stress_xy = stress_xy_lam
        return self._stress_xy

    @property 
    def stress_12(self):
        ''' [sigma_1, sigma_2, tau_12] for each ply'''
        if self._stress_12 is None:
            stress_12_lam = []
            ply_num = 0
            for T_ply in self.laminate.T_layup:
                stress_12_ply = T_ply[0] * self.stress_xy[ply_num]
                stress_12_lam.append(stress_12_ply)
                ply_num += 1
            self._stress_12 = stress_12_lam
        return self._stress_12    
    
    @property
    def TW_i(self):
        ''' Tsai-Wu index (F12 was aproximated using Hoffman criterion)'''
        if self._TW_i is None:
            TW_i_lam = []
            ply_num = 0
            for X in self.laminate.layup:
                F1 = 1/X[1].t1 + 1/X[1].c1
                F2 = 1/X[1].t2 + 1/X[1].c2
                F11 = -1/(X[1].t1*X[1].c1)
                F22 = -1/(X[1].t2*X[1].c2)
                F66 = 1/X[1].s**2
                F12 = 1/(2*X[1].t1*X[1].c1) 
                TW_i = float(F1 * self.stress_12[ply_num][0]
                             + F2 * self.stress_12[ply_num][1] 
                             + F11 * self.stress_12[ply_num][0]**2
                             + F22 * self.stress_12[ply_num][1]**2 
                             + F66 * self.stress_12[ply_num][2]**2 
                             + 2 * F12 * self.stress_12[ply_num][0] * self.stress_12[ply_num][1])
                TW_i_lam.append(TW_i)
                ply_num += 1
            self._TW_i = TW_i_lam
        return self._TW_i

    @property 
    def FPF(self):
        ''' First Ply Failure analysis '''
        if self._FPF is None:
            i = 1
            failed_plies = []
            while len(failed_plies) == 0:
                ABD_p = self.laminate.ABD_p
                load_array = i * sp.matrix(np.vstack((self.load.Nx,
                                                      self.load.Ny,
                                                      self.load.Nxy,
                                                      self.load.Mx,
                                                      self.load.My,
                                                      self.load.Mxy)))
                mid_strain_xy = ABD_p * load_array
                stress_xy_lam = []
                stress_12_lam = []
                TW_i_lam = []
                for x in range(len(self.laminate.layup)):
                    stress_xy_ply = self.laminate.Q_layup[x]
                                    * mid_strain_xy[0:3] 
                                    + self.laminate.z_position[x]
                                    * mid_strain_xy[3:7]
                    stress_xy_lam.append(stress_xy_ply)
                    stress_12_ply = self.laminate.T_layup[x][0] * stress_xy_ply
                    stress_12_lam.append(stress_12_ply)

                    F1 = 1 / self.laminate.layup[x][1].t1 
                        + 1 / self.laminate.layup[x][1].c1
                    F2 = 1 / self.laminate.layup[x][1].t2 
                        + 1 / self.laminate.layup[x][1].c2
                    F11 = -1 / (
                        self.laminate.layup[x][1].t1
                        * self.laminate.layup[x][1].c1
                        )
                    F22 = -1 / (
                        self.laminate.layup[x][1].t2
                        * self.laminate.layup[x][1].c2
                        )
                    F66 = 1 / self.laminate.layup[x][1].s**2
                    F12 = 1 / (
                        2 * self.laminate.layup[x][1].t1
                        * self.laminate.layup[x][1].c1
                        ) 
                    
                    TW_i = float(F1*stress_12_ply[0] 
                                 + F2*stress_12_ply[1] 
                                 + F11*stress_12_ply[0]**2 
                                 + F22*stress_12_ply[1]**2 
                                 + F66*stress_12_ply[2]**2
                                 + 2*F12*stress_12_ply[0] * stress_12_ply[1])
                    TW_i_lam.append(TW_i)

                for j in range(len(TW_i_lam)):
                    if TW_i_lam[j] > 1:
                        failed_plies.append(j+1)

                i = i * 1.01
                if i > 100: break

            self._SR = i
            print("FPF Load: \n", load_array)
            print("Failed Ply(ies): \n", failed_plies)
            print("Laminate Tsai-Wu Index in FPF: \n", TW_i_lam)  
            print("Strength Ratio for FPF: \n", self._SR)
        return None

    #Representation
    def __repr__(self):
        representation = ''
        for  ply_num in range(len(self.laminate.layup)):
            representation += f'Ply {ply_num+1} \nStress_12 \n{self.stress_12[ply_num]} \nTsai-Wu Index: {self.TW_i[ply_num]} \n ============== \n'
        return representation

#Comparisons
    def __eq__(self, other):
        if isinstance(other, Strength):
            return (self.laminate == other.laminate
                    and self.load == other.load)
        return NotImplemented
