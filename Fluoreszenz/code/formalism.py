import numpy as np

"""
States are represented in the HV basis of states:
    ( |H>|H> )
    ( |V>|V> )
    ( |H>|V> )
    ( |V>|H> )
"""

# Definitions
def R(a):
    return np.matrix([[ np.cos(a), np.sin(a)],
                      [-np.sin(a), np.cos(a)]])

def J(a):
    e   = np.matrix([[1.,0.],[0.,0.]])
    return R(a)*e*R(a).transpose()

def A(a,b):
    Ja  = J(a)
    Jb  = J(b)
    return np.matrix([[Ja[0,0]*Jb[0,0],Ja[0,1]*Jb[0,1],Ja[0,0]*Jb[0,1],Ja[0,1]*Jb[0,0]],
                      [Ja[1,0]*Jb[1,0],Ja[1,1]*Jb[1,1],Ja[1,0]*Jb[1,1],Ja[1,1]*Jb[1,0]],
                      [Ja[0,0]*Jb[1,0],Ja[0,1]*Jb[1,1],Ja[0,0]*Jb[1,1],Ja[0,1]*Jb[1,0]],
                      [Ja[1,0]*Jb[0,0],Ja[1,1]*Jb[0,1],Ja[1,0]*Jb[0,1],Ja[1,1]*Jb[0,0]]])

def C(a,b,rho):
    # Tr(rho*A)
    op  = rho*A(a,b)
    return op.trace()[0,0]

def E(a,b,rho):
    C1  = C(a,b,rho)
    C2  = C(a,b+np.pi/2.,rho)
    C3  = C(a+np.pi/2.,b,rho)
    C4  = C(a+np.pi/2.,b+np.pi/2.,rho)
    head    =  (C1-C2-C3+C4)
    foot    =  (C1+C2+C3+C4)
    return head/foot

def S(a1,b1,a2,b2,rho):
    return E(a1,b1,rho)-E(a1,b2,rho)+E(a2,b1,rho)+E(a2,b2,rho)

rad = np.pi/180.
a1  = 22.5*rad
b1  = 0.*rad
a2  = -22.5*rad
b2  = -45.*rad

phip    = np.matrix([[1.],
                     [1.],
                     [0.],
                     [0.]])/np.sqrt(2.)
phim    = np.matrix([[1.],
                     [-1.],
                     [0.],
                     [0.]])/np.sqrt(2.)
psip    = np.matrix([[0.],
                     [0.],
                     [1.],
                     [1.]])/np.sqrt(2.)
psim    = np.matrix([[0.],
                     [0.],
                     [1.],
                     [-1.]])/np.sqrt(2.)

bell    = [phip,phim,psip,psim]

bellmat = np.matrix([[1.,1.,0.,0.],
                     [1.,-1.,0.,0.],
                     [0.,0.,1.,1.],
                     [0.,0.,1.,-1.]])/np.sqrt(2.)


print("2**1.5   = {}".format(2**1.5))
print("")
for k in range(4):
    print("bell = {}".format(bell[k].H))
    rho = bell[k]*bell[k].H
    #print("rho")
    #print(rho)
    print("S    = {}".format(S(a1,b1,a2,b2,rho)))
    print()

"""
rho = np.zeros((4,4))
ps  = np.array([54,10,4,300])
ps  = ps/np.linalg.norm(ps)
print("ps   = {}".format(ps))
for k in (0,1,2,3):
    rho += bell[k]*bell[k].H*ps[k]
print("S    = {}".format(S(a1,b1,a2,b2,rho)))
print()
"""

# Basic states:
pol_dict    = {
    "H" : np.array((1.,0.)),
    "V" : np.array((0.,1.)),
    "D" : np.array((1.,1.))*2**-0.5,
    "R" : np.array((1.,-1.j))*2**-0.5,
    "L" : np.array((1.,1.j))*2**-0.5
    }

def basis(P1,P2):
    """
    Give Photon states P1 and P2 as H, V, D, R or L
    Returns expression as matrix
    """
    P1s = pol_dict[P1]
    P2s = pol_dict[P2]
    return np.matrix([[P1s[0]*P2s[0]],
                      [P1s[1]*P2s[1]],
                      [P1s[0]*P2s[1]],
                      [P1s[1]*P2s[0]]])
