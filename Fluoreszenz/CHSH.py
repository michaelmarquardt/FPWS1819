import numpy as np
from HV_formalism import *

# Definitions
def R(a):
    return np.matrix([[ np.cos(a), np.sin(a)],
                      [-np.sin(a), np.cos(a)]])

def J(a):
    e   = np.matrix([[1.,0.],[0.,0.]])
    return R(a)*e*R(a).transpose()

def A(a,b):
    return Operator(J(a),J(b))

def C(a,b,rho):
    print("calculate C")
    return rho.matmul(A(a,b)).trace()

def E(a,b,rho):
    C1  = C(a,b,rho)
    C2  = C(a,b+np.pi/2.,rho)
    C3  = C(a+np.pi/2.,b,rho)
    C4  = C(a+np.pi/2.,b+np.pi/2.,rho)
    head    =  (C1-C2-C3+C4)
    foot    =  (C1+C2+C3+C4)
    return head/foot

def S(a1,b1,a2,b2,rho):
    return E(a1,a2,rho)-E(a1,b2,rho)+E(a2,b1,rho)+E(a2,b2,rho)

rad = np.pi/180.
a1  = 22.5*rad
b1  = 0.*rad
a2  = -22.5*rad
b2  = -45.*rad

phip    = State()
phim    = State()
psip    = State()
psim    = State()

r2  = 1./np.sqrt(2.)

phip.add(HV_State([r2,0.],[r2,0.]))
phip.add(HV_State([0.,r2],[0.,r2]))
phim.add(HV_State([r2,0.],[r2,0.]))
phim.add(HV_State([0.,-r2],[0.,-r2]))
psip.add(HV_State([r2,0.],[0.,r2]))
psip.add(HV_State([0.,r2],[r2,0.]))
psim.add(HV_State([r2,0.],[0.,r2]))
psim.add(HV_State([0.,-r2],[-r2,0.]))

rhopsip = psip.ketbra(psip)

print(S(a1,b1,a2,b2,rhopsip))
    
    
