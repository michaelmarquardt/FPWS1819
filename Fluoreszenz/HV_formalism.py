import numpy as np

# Basis
HH  = HV_State([1.,0.],[1.,0.])
HV  = HV_State([1.,0.],[0.,1.])
VH  = HV_State([0.,1.],[1.,0.])
VV  = HV_State([0.,1.],[0.,1.])
basis   = [HH,HV,VH,VV]

class HV_State:
    """
    A 2 photon state in (H,V) (horizontal, vertical) expression:
        (Vectors )
        H = (1,0)
        V = (0,1)
    example: |H>x|V>:
        bra1    = (1,0)
        bra2    = (0,1)
        ket1    = (1,0)T (transposed)
        ket2    = (0,1)T
    """
    
    def __init__(self,s1,s2):
        """
        Takes two list like variables for photon 1 and 2
        """
        self.bra1 = np.matrix([s1])
        self.bra2 = np.matrix([s2])
        self.ket1 = self.bra1.getH()
        self.ket2 = self.bra2.getH()
    
    def braket(self,HV_state):
        """
        <self|HV_state>
        """
        return self.bra1*HV_state.ket1 + self.bra2*HV_state.ket2
    
    def ketbra(self,HV_state):
        """
        |self><HV_state|
        """
        O1  = self.ket1*HV_state.bra1
        O2  = self.ket2*HV_state.bra2
        return HV_Operator(O1,O2)
    
    def toket(self,HV_operator):
        """
        HV_operator|self> and <self|HV_operator*
        """
        self.bra1   = HV_operator.O1*self.bra1
        self.bra2   = HV_operator.O2*self.bra2
        self.ket1 = self.bra1.getH()
        self.ket2 = self.bra2.getH()
    
    def tobra(self,HV_operator):
        """
        <self|HV_operator and HV_operator*|self> 
        """
        self.ket1   = self.ket1*HV_operator.O1
        self.ket2   = self.ket2*HV_operator.O2
        self.bra1 = self.ket1.getH()
        self.bra2 = self.ket2.getH()

class State:
    """
    A complex state expressed in a superposition of HV_States
    """
    hv  = []
    
    def __init__(self,hv_list=[]):
        """
        Takes list of HV_states (list is important)
        """
        for HV_state in hv_list:
            self.add(HV_state)
    
    def add(self,HV_state):
        """
        Adds a purestate to the superposition hv
        """
        self.hv.add(HV_state)
    
    def braket(self,state):
        """
        Calculates the <bra|ket> with <self,state>
        """
        result  = 0.
        for hv1 in self.hv:
            for hv2 in state.hv:
                result += hv1.braket(hv2)
        return result
    
    def ketbra(self,state):
        """
        Calculates the |ket><bra| with |self><state|
        It returns a Ketbra
        """
        result  = Operator()
        for hv1 in self.hv:
            for hv2 in state.hv:
                result.add(hv1.ketbra(hv2))
        return result
    
    def toket(self,operator):
        """
        operator|self> and <self|operator*
        """
        result  = State()
        for hv1 in self.hv:
            for hv2 in operator.hv:
                result.add(hv1.toket(hv2))
        return result
    
    def tobra(self,operator):
        """
        <self|operator and operator*|self> 
        """
        result  = State()
        for hv1 in self.hv:
            for hv2 in operator.hv:
                result.add(hv1.tobra(hv2))
        return result
   
    def contract(self):
        """
        Contracts the list hv to a minimal set
        """
        #TODO
        """
        new_hv  = []
        old_hv  = self.hv
        k       = 0
        while k < len(self.hv)-1:
            for l in range(k+1,len(self.hv)):
        """     
    

class HV_Operator:
    """
    A operator which can be expressed in the HV_State formalism
    """
    def __init__(self,O1,O2):
        """
        Operators for photon 1 and 2 as (2x2)-matrix (ketbra)
        """
        self.O1 = O1
        self.O2 = O2
    
    def matmul(self,HV_operator):
        """
        self*HV_operator
        """
        O1  = self.O1*HV_operator.O1
        O2  = self.O2*HV_operator.O2
        return HV_Operator(O1,O2)

class Operator:
    """
    A complex operator expressed in a superposition of HV_Operators
    """
    hv  = []
    
    def __init__(self,hv_list=[]):
        """
        Takes list of HV_states (list is important)
        """
        for HV_operator in hv_list:
            self.add(HV_operator)
    
    def add(self,HV_operator):
        """
        Adds a purestate to the superposition hv
        """
        self.hv.add(HV_operator)
    
    def matmul(self,operator):
        """
        self*operator
        """
        result  = Operator()
        for hv1 in self.hv:
            for hv2 in operator.hv:
                result.add(hv1.matmul(hv2))
        return result
    
    def trace(self):
        result  = 0.
        for hv1 in self.hv:
            for hv2 in basis:
                result += hv2.braket(hv2.toket(hv1))
        return result
