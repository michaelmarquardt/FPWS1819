import numpy as np
from HV_formalism import *

# Definitions
# Basis
HH  = HV_State([1.,0.],[1.,0.])
HV  = HV_State([1.,0.],[0.,1.])
VH  = HV_State([0.,1.],[1.,0.])
VV  = HV_State([0.,1.],[0.,1.])
basis   = [HH,HV,VH,VV]


def R(a):
    return np.matrix([[ np.cos(a), np.sin(a)],
                     [-np.sin(a), np.cos(a)]])

def J(a):
    e   = np.matrix([[1.,0.],[0.,0.]])
    return R*J*R.transpose()

def C(a,b,rho):
    
    
    
    
