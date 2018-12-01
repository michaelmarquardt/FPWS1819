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

#   Define constants
#######################################################################################
mAr = 39.948    # u (mass of argon)
atu = con.physical_constants["atomic mass constant"][0] #kg/u
si("mAr",mAr,None,"u",3)
mAr = atu*mAr   # kg
S   = 1.        # mm2
si("S",S,None,"mm2",0)
S   = S*10**-6  # m2
Troom   = 300.  # K
si("Troom",Troom,None,"K",0)

#   Define errors
#######################################################################################
#U, I    = getdata("error1.dat")
dU  = 0.01      # V
dI  = 1.        # muA
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
    return A*np.exp(U/Te) + Iisat   # Te in eV

def fit_single(U, I, dI, U_min=-200., U_max=0.):
    r_use   = np.where(U>=U_min)
    r_use   = np.where(U[r_use]<=U_max)
    U_use   = U[r_use]
    I_use   = I[r_use]
    dI_use  = np.ones_like(I_use)*dI
    # Initial guess
    U_e     = U_use[np.where(I_use<=I_use[0]/np.e)][0]
    Te0     = (U_use[0]-U_e)
    Iisat0  = I_use[-1]
    A0      = (I_use[0]-Iisat0)*np.exp(-U_use[0]/Te0)
    p0      = [Te0, Iisat0, A0]
    # Fit
    popt, pcov  = curve_fit(I_single, U_use, I_use, sigma=dI_use, p0=p0)
    return popt, np.sqrt(np.diag(pcov))

#   Single
#######################################################################################
set_pre("A")
ps      = np.array([0.3, 1.,5.,10.,30.])      # mbar
datas   = ["single_Ar_7.dat",
           "single_Ar_5.dat",
           "single_Ar_4.dat",
           "single_Ar_8.dat",
           "single_Ar_6.dat"]
U_maxs  = [-39., -55., -42., -70., -65.]
Tes     = []
dTes    = []
Iisats  = []
dIisats = []
nes     = []
dnes    = []
alphas  = []
dalphas = []
debs    = []
ddebs   = []
wpes    = []
dwpes   = []
wpis    = []
dwpis   = []

# Start main loop
f   = plt.figure()
for p, data, U_max, a in zip(ps, datas, U_maxs, alphabet):
    U, I    = getdata(data)
    # Exponential fit to the curves
    popt, perr  = fit_single(U, I, dI, U_max=U_max)
    si("p"+a,p,dp,"mbar",1)
    si("Te"+a,popt[0],perr[0],"eV",2)
    si("Iisat"+a,popt[1],perr[1],"muA",2)
    si("A"+a,popt[2]/1000.,perr[2]/1000.,"mA",2)
    Tes.append(popt[0])
    dTes.append(perr[0])
    Iisats.append(popt[1])
    dIisats.append(perr[1])
    
    # Calculate electron density
    ne  = -popt[1]*np.sqrt(mAr/con.e/popt[0])/0.61/con.e/S   #1_m3
    dne = np.abs(perr[1]/popt[1]*ne) + np.abs(0.5*perr[0]/popt[0]*ne)
    si("ne"+a,ne,dne,"1_m3",2,e=0)
    nes.append(ne)
    dnes.append(dne)
    
    # Calculate degree of ionisation
    nG  = 10**8*p/con.k/Troom
    dnG = 10**8*dp/con.k/Troom
    si("nG"+a,nG,dnG,"1_m3",2,e=0)
    alpha   = ne/nG
    dalpha  = ne*dnG/nG**2
    si("alpha"+a,alpha,dalpha,"",2,e=0)
    alphas.append(alpha)
    dalphas.append(dalpha)
    
    # Calculate debye length and plasma frequency
    deb     = np.sqrt(con.epsilon_0*popt[0]/ne/con.e)*10**9    # nm
    ddeb    = np.abs(0.5*deb*perr[0]/popt[0]) + np.abs(0.5*deb*dne/ne)
    si("deb"+a,deb,ddeb,"nm",0)
    debs.append(deb)
    ddebs.append(ddeb)
    wpe     = np.sqrt(4*np.pi*ne*con.e**2/con.m_e)*10**-6   # MHz
    dwpe    = np.abs(0.5*wpe*dne/ne)
    si("wpe"+a,wpe,dwpe,"MHz",0)
    wpes.append(wpe)
    dwpes.append(dwpe)
    wpi     = np.sqrt(4*np.pi*ne*con.e**2/mAr)*10**-3   # kHz
    dwpi    = np.abs(0.5*wpi*dne/ne)
    si("wpi"+a,wpi,dwpi,"MHz",0)
    wpis.append(wpi)
    dwpis.append(dwpi)
    
    # Plot of all characteristic curves
    x           = np.linspace(-150,-20,1000)
    p1  = plt.plot(x, I_single(x, *popt), zorder=3)
    plt.errorbar(U, I, xerr=dU, yerr=dI, 
        fmt='.', 
        label=r"$p=\SI{{{:.1f}}}{{\milli\bar}}$".format(p), 
        color=changecolor(p1[0].get_color(),rgb=(40,40,40)),
        zorder=2
        )
    plt.axvline(U_max,linestyle='--',color=changecolor(p1[0].get_color(),rgb=(40,40,40)))

plt.xlabel(r"$U$ [\si{\volt}]")
plt.ylabel(r"$I$ [\si{\micro\ampere}]")
plt.xlim((-130,-20))
plt.ylim((-50,1050))
plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_KL.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_KL.pdf", bbox_inches="tight")
plt.close()

alphas  = np.array(alphas)
dalphas = np.array(dalphas)

f   = plt.figure()
perr    = np.ones((2,5))*dp
perr[0,perr[0,:]>=ps] = ps[perr[0,:]>=ps]*0.99999999
plt.errorbar(ps, Tes, yerr=dTes, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$T_e$ [\si{\electronvolt}]")
plt.xscale("log")
plt.xlim((0.2,40))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_Te.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_Te.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, nes, yerr=dnes, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$n_e$ [\si{\per\metre\tothe{3}}]")
plt.xscale("log")
plt.xlim((0.2,40))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_ne.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_ne.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
aerr   = np.array([dalphas,dalphas])
aerr[0,aerr[0,:]>=alphas] = alphas[aerr[0,:]>=alphas]*0.99999999
plt.errorbar(ps, alphas, yerr=aerr, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\alpha")
plt.xscale("log")
plt.xlim((0.2,40))
plt.ylim((10**-8,10**-4))
plt.yscale("log")
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_alpha.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_alpha.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, debs, yerr=ddebs, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\lambda_\text{D}$ [\si{\nano\metre}]")
plt.xscale("log")
plt.xlim((0.2,40))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_deb.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_deb.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, wpes, yerr=dwpes, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\omega_{\text{p},e}$ [\si{\mega\hertz}]")
plt.xscale("log")
plt.xlim((0.2,40))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_wpe.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_wpe.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, wpis, yerr=dwpis, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\omega_{\text{p},i}$ [\si{\mega\hertz}]")
plt.xscale("log")
plt.xlim((0.2,40))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_wpi.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_wpi.pdf", bbox_inches="tight")
plt.close()
