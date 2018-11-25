import sys
sys.path.append('../../pylibs')
from export_data import *
texfile("../report/values.tex")
set_dig(3)
set_pre("")
usepgf()
import matplotlib.pyplot as plt

import numpy as np
from scipy.optimize import curve_fit
import scipy.constants as con
import subprocess

DATA    = "../Messdaten/"
PLOT    = "../figures/"
subprocess.call(["mkdir","-p",PLOT])

#   Define errors
#######################################################################################
#U, I    = getdata("error1.dat")
dU  = 0.01      # V
dI  = 1.       # muA
dp  = 5.        # mbar

si("dU",dU,None,"V",2)
si("dI",dI,None,"muA",0)
si("dp",dp,None,"mbar",0)

#   Functions
#######################################################################################
def getdata(name):
    U, I    = np.genfromtxt(open(DATA+name), usecols=(1,2), skip_header=1, unpack=True)
    return U, I*10**6   # U in V and I in muA

def I_single(U, Te, Iisat, A):
    return A*np.exp(con.e*U/Te) + Iisat

def fit_single(U, I, dI, U_min=-200., U_max=0.):
    r_use   = np.where(U>=U_min)
    r_use   = np.where(U[r_use]<=U_max)
    U_use   = U[r_use]
    I_use   = I[r_use]
    dI_use  = np.ones_like(I_use)*dI
    # Initial guess
    U_e     = U_use[np.where(I_use<=I_use[0]/np.e)][0]
    Te0     = con.e*(U_use[0]-U_e)
    Iisat0  = I_use[-1]
    A0      = (I_use[0]-Iisat0)*np.exp(-con.e*U_use[0]/Te0)
    p0      = [Te0, Iisat0, A0]
    # Fit
    popt, pcov  = curve_fit(I_single, U_use, I_use, sigma=dI_use, p0=p0)
    return popt, np.diag(popt)

def changecolor(color, rgb=(1,1,1)):
    R   = int(color[1:3],16)
    G   = int(color[3:5],16)
    B   = int(color[5:7],16)
    R   = min(255,R+rgb[0])
    R   = max(R,0)
    G   = min(255,G+rgb[1])
    G   = max(G,0)
    B   = min(255,B+rgb[2])
    B   = max(B,0)
    return "#{:02x}{:02x}{:02x}".format(R,G,B)

#   Single
#######################################################################################
ps      = [0.3, 1.,5.,10.,30.]      # mbar
datas   = ["single_Ar_7.dat",
           "single_Ar_5.dat",
           "single_Ar_4.dat",
           "single_Ar_8.dat",
           "single_Ar_6.dat"]
U_maxs  = [-39., -55., -42., -70., -65.]

# Start main loop
for p, data, U_max in zip(ps, datas, U_maxs):
    U, I    = getdata(data)
    # Exponential fit to the curves
    popt, perr  = fit_single(U, I, dI, U_max=U_max)
    x           = np.linspace(-150,-20,1000)
    
    # Plot of all characteristic curves
    p1  = plt.plot(x, I_single(x, *popt), zorder=3)
    plt.errorbar(U, I, xerr=dU, yerr=dI, 
        fmt='.', 
        label=r"$p=\SI{{{:.1f}}}{{\milli\bar}}$".format(p), 
        color=changecolor(p1[0].get_color(),rgb=(40,40,40)),
        zorder=2
        )
    plt.axvline(U_max,linestyle='--',color=changecolor(p1[0].get_color(),rgb=(40,40,40)))
    
plt.xlabel(r"$U$ [\si{\volt}]")
plt.ylabel(r"$I$ [\si{\ampere}]")
plt.xlim((-130,-20))
plt.ylim((-50,1050))
plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_KL.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_KL.pdf", bbox_inches="tight")
plt.close()

