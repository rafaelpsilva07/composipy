import numpy as np
import numbers

class Ply:
    '''This function builts a ply object'''
    
    def __init__(self, e1, e2, v12, g12, thickness):

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

        self._e1 = e1
        self._e2 = e2
        self._v12 = v12
        self._g12 = g12
        self._thickness = thickness
        self._Q_0 = None #Calculated property
        self._v21 = None #Calculated property

#Properties getters and setters
    @property
    def e1(self):
        return self._e1
    
    @e1.setter
    def e1(self,value):
        self._e1 = value
        self._Q_0 = None
        self._v21 = None

    @property
    def e2(self):
        return self._e2
    
    @e2.setter
    def e2(self,value):
        self._e2 = value
        self._Q_0 = None
        self._v21 = None
        
    @property
    def v12(self):
        return self._v12
    
    @v12.setter
    def v12(self,value):
        self._v12 = value
        self._Q_0 = None
        self._v21 = None

    @property
    def g12(self):
        return self._g12
    
    @g12.setter
    def g12(self,value):
        self._g12 = value
        self._Q_0 = None
        self._v21 = None

        
    @property
    def thickness(self):
        return self._thickness
    
    @thickness.setter
    def thickness(self,value):
        self._thickness = value
        self._Q_0 = None
        self._v21 = None

    #Calculated properties (have only getter)
    @property
    def v21(self):
        if self._v21 is None:
            #value not cached - calculate it
            self._v21 = self.v12*self.e2/self.e1
        return self._v21
    
    @property
    def Q_0(self):
        if self._Q_0 is None:
            self._Q_0 = np.array([[self.e1/(1-self.v12*self.v21),self.v12*self.e2/(1-self.v12*self.v21),0],
                                        [self.v12*self.e2/(1-self.v12*self.v21),self.e2/(1-self.v12*self.v21),0],
                                        [0,0,self.g12]])
        return self._Q_0
    
    def __repr__(self):
        return f'Ply(e1={self.e1}, e2={self.e2}, e3={self.v12},\
 g12={self.g12}, thickness={self.thickness}) object'
    
#Comparisons
    def __eq__(self, other):
        if isinstance(other, Ply):
            return (self.e1 == other.e1
                    and self.e2 == other.e2
                    and self.v12 == other.v12
                    and self.g12 == other.g12
                    and self.thickness == other.thickness)
        return NotImplemented

