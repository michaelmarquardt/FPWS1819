import sys
sys.path.append('../../pylibs')
from export_data import *
texfile("../report/values.tex")
set_dig(3)
set_pre("")
usepgf(size=(1.3,3.))
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

# Constants
si("k",con.k/con.e,None,"eV_K",3,e=0)
si("e",con.e,None,"C",3,e=0)
si("me",con.m_e,None,"kg",3,e=0)
si("eO",con.epsilon_0,None,"As_Vm",3,e=0)
si("u",atu,None,"kg",3,e=0)

#   Functions
#######################################################################################
def getdata(name):
    U, I    = np.genfromtxt(open(DATA+name), usecols=(1,2), skip_header=1, unpack=True)
    return U, I*10**6   # U in V and I in muA

def getdata_cut(name):
    U, I    = np.genfromtxt(open(DATA+name), usecols=(1,2), skip_header=1, unpack=True)
    I      *= 10**6
    r_use   = np.where(np.logical_or((U>-60.),(I<400)))
    return U[r_use], I[r_use]   # U in V and I in muA

def I_single(U, Te, Iisat, A):
    return A*np.exp(U/Te) + Iisat   # Te in eV

def fit_single(U, I, dI, U_min=-200., U_max=0.):
    r_use   = np.where(U>=U_min)
    U_use   = U[r_use]
    I_use   = I[r_use]
    r_use   = np.where(U_use<=U_max)
    U_use   = U_use[r_use]
    I_use   = I_use[r_use]
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

def linear(x, c0, c1):
    return c1*x + c0

def fit_double(U, I, dI, U_min, U_max):
    r_use   = np.where(U>=U_min)
    U_use   = U[r_use]
    I_use   = I[r_use]
    r_use   = np.where(U_use<=U_max)
    U_use   = U_use[r_use]
    I_use   = I_use[r_use]
    dI_use  = np.ones_like(I_use)*dI
    # Initial guess
    c1  = (I_use[-1]-I_use[0])/(U_use[-1]-U_use[0])
    c0  = (I_use[-1]-c1*U_use[-1])
    p0  = [c0, c1]
    # Fit
    popt, pcov  = curve_fit(linear, U_use, I_use, sigma=dI_use, p0=p0)
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
U_maxs  = [-39., -55., -42., -62., -65.]
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
    si("Iisat"+a,np.abs(popt[1]),perr[1],"muA",2)
    si("A"+a,popt[2]/1000.,perr[2]/1000.,"mA",0)
    Tes.append(popt[0])
    dTes.append(perr[0])
    Iisats.append(popt[1])
    dIisats.append(perr[1])
    siline("popt"+a,[p, dp, popt[0], perr[0], popt[1], perr[1], popt[2]/1000., perr[2]/1000.],[1,2,2,0])
    
    # Calculate electron density
    ne  = -popt[1]*np.sqrt(mAr/con.e/popt[0])/0.61/con.e/S*10**-6       #1_m3
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
    deb     = np.sqrt(con.epsilon_0*popt[0]/ne/con.e)*10**6    # mum
    ddeb    = np.abs(0.5*deb*perr[0]/popt[0]) + np.abs(0.5*deb*dne/ne)
    si("deb"+a,deb,ddeb,"mum",0)
    debs.append(deb)
    ddebs.append(ddeb)
    wpe     = np.sqrt(ne*con.e**2/con.epsilon_0/con.m_e)*10**-9   # GHz
    dwpe    = np.abs(0.5*wpe*dne/ne)
    si("wpe"+a,wpe,dwpe,"GHz",2)
    wpes.append(wpe)
    dwpes.append(dwpe)
    wpi     = np.sqrt(ne*con.e**2/con.epsilon_0/mAr)*10**-6   # MHz
    dwpi    = np.abs(0.5*wpi*dne/ne)
    si("wpi"+a,wpi,dwpi,"MHz",2)
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
#print(alphas)
aerr[0,aerr[0,:]>=alphas] = alphas[aerr[0,:]>=alphas]*0.99999999
plt.errorbar(ps, alphas, yerr=aerr, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\alpha")
plt.xscale("log")
plt.xlim((0.2,40))
plt.ylim((10**-14,10**-10))
plt.yscale("log")
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_alpha.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_alpha.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, debs, yerr=ddebs, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\lambda_\text{D}$ [\si{\micro\metre}]")
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
plt.ylabel(r"$\omega_{e}$ [\si{\giga\hertz}]")
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
plt.ylabel(r"$\omega_{i}$ [\si{\mega\hertz}]")
plt.xscale("log")
plt.xlim((0.2,40))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"single_wpi.pgf", bbox_inches="tight")
plt.savefig(PLOT+"single_wpi.pdf", bbox_inches="tight")
plt.close()

#   Double
#######################################################################################
# Special plot
f   = plt.figure()
plt.axvline(0,linestyle='--',color='k')
plt.axhline(0,linestyle='--',color='k')
U, I    = getdata("double_Ar_1.dat")
plt.errorbar(U, I, xerr=dU, yerr=dI, 
    fmt='.',
    )
plt.xlabel(r"$U$ [\si{\volt}]")
plt.ylabel(r"$I$ [\si{\micro\ampere}]")
#plt.xlim((-130,-20))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"double_KL_1.pgf", bbox_inches="tight")
plt.savefig(PLOT+"double_KL_1.pdf", bbox_inches="tight")
plt.close()

set_pre("B")
ps      = np.array([0.7, 1.,3.,5.,10.])      # mbar
datas   = ["double_Ar_6.dat",
           "double_Ar_3.dat",
           "double_Ar_2.dat",
           "double_Ar_4.dat",
           "double_Ar_5.dat"]
U_mins  = np.array([[-80., -80., -80., -80., -80.],
                    [ -5.,  -5.,  -5.,  -5.,  -5.],
                    [ 20.,  20.,  20.,  20.,  20.]]).transpose()
U_maxs  = np.array([[-25., -25., -25., -25., -25.],
                    [  5.,   5.,   5.,   5.,   5.],
                    [100., 100., 100., 100., 100.]]).transpose()
#U_min   = [-80., -5., 20.]
#U_max   = [-25., 5., 100.]
Tes         = []
dTes        = []
Ii1sats     = []
dIi1sats    = []
Ii2sats     = []
dI21sats    = []
dIdUs       = []
ddIdUs      = []
nes         = []
dnes        = []
alphas      = []
dalphas     = []
debs        = []
ddebs       = []
wpes        = []
dwpes       = []
wpis        = []
dwpis       = []

# Start main loop
f   = plt.figure()
plt.axvline(0,linestyle='--',color='k')
plt.axhline(0,linestyle='--',color='k')

for p, data, U_min, U_max, a in zip(ps, datas, U_mins, U_maxs, alphabet):
    U, I    = getdata(data)
    
    # Linear fits to the curves
    popt0, perr0    = fit_double(U, I, dI, U_min=U_min[0], U_max=U_max[0])
    popt1, perr1    = fit_double(U, I, dI, U_min=U_min[1], U_max=U_max[1])
    popt2, perr2    = fit_double(U, I, dI, U_min=U_min[2], U_max=U_max[2])
    
    I1  = np.abs(popt0[0])
    I2  = np.abs(popt2[0])
    dI1 = perr0[0]
    dI2 = perr2[0]
    dIdU    = popt1[1]
    ddIdU   = perr1[1]
    si("p"+a,p,dp,"mbar",1)
    si("Iiasat"+a,I1,dI1,"muA",2)
    si("Iibsat"+a,I2,dI2,"muA",2)
    si("Sa"+a,popt0[1],perr0[1],"muA_V",3)
    si("Sb"+a,popt2[1],perr2[1],"muA_V",3)
    si("dIdU"+a,dIdU,ddIdU,"muA_V",2)
    si("A"+a,popt1[0],perr1[0],"muA",2)
    Ii1sats.append(I1)
    dIi1sats.append(dI1)
    Ii2sats.append(I2)
    dI21sats.append(dI2)
    dIdUs.append(dIdU)
    ddIdUs.append(ddIdU)
    siline("popt"+a,[p, dp, 
                     I1, dI1, 
                     popt0[1], perr0[1],
                     I2, dI2, 
                     popt2[1], perr2[1],
                     popt1[0], perr1[0], 
                     popt1[1], perr1[1],
                     ],[1,2,3,2,3,2,2])
    
    # Calculate electron temperature
    Te  = I1*I2/(I1+I2)/dIdU        # eV
    dTe = np.abs(I2**2/(I1+I2)**2/dIdU*dI1) + np.abs(I1**2/(I1+I2)**2/dIdU*dI2) + np.abs(ddIdU*Te/dIdU)
    si("Te"+a,Te,dTe,"eV",2)
    Tes.append(Te)
    dTes.append(dTe)
    
    # Calculate electron density
    ne  = I2*np.sqrt(mAr/con.e/Te)/0.61/con.e/S*10**-6   #1_m3
    dne = np.abs(dI2/I2*ne) + np.abs(0.5*dTe/Te*ne)
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
    deb     = np.sqrt(con.epsilon_0*Te/ne/con.e)*10**6    # mum
    ddeb    = np.abs(0.5*deb*dTe/Te) + np.abs(0.5*deb*dne/ne)
    si("deb"+a,deb,ddeb,"mum",0)
    debs.append(deb)
    ddebs.append(ddeb)
    wpe     = np.sqrt(ne*con.e**2/con.epsilon_0/con.m_e)*10**-9   # GHz
    dwpe    = np.abs(0.5*wpe*dne/ne)
    si("wpe"+a,wpe,dwpe,"GHz",2)
    wpes.append(wpe)
    dwpes.append(dwpe)
    wpi     = np.sqrt(ne*con.e**2/con.epsilon_0/mAr)*10**-6   # MHz
    dwpi    = np.abs(0.5*wpi*dne/ne)
    si("wpi"+a,wpi,dwpi,"MHz",2)
    wpis.append(wpi)
    dwpis.append(dwpi)
    
    # Plot of all characteristic curves
    if a == "a" or a == "c" or a == "e":
        x0  = np.array([-130.,  10.])
        x1  = np.array([-130., 130.])
        x2  = np.array([ -10., 130.])
        p1  = plt.plot(x0, linear(x0, *popt0), zorder=3)
        plt.plot(x1, linear(x1, *popt1), zorder=2, color=p1[0].get_color())
        plt.plot(x2, linear(x2, *popt2), zorder=1, color=p1[0].get_color())
        plt.errorbar(U, I, xerr=dU, yerr=dI, 
            fmt='.', 
            label=r"$p=\SI{{{:.1f}}}{{\milli\bar}}$".format(p), 
            color=changecolor(p1[0].get_color(),rgb=(40,40,40)),
            zorder=0
            )

plt.annotate(r"$S_1 U -I_{i,1,\text{sat}}$", xy=(-100, -35))
plt.annotate(r"$\frac{\dd I}{\dd U}\Big\lvert_\text{fl} U -I_\text{off}$", xy=(20, 5))
plt.annotate(r"$S_2 U +I_{i,2,\text{sat}}$", xy=(55, 35))
plt.xlabel(r"$U$ [\si{\volt}]")
plt.ylabel(r"$I$ [\si{\micro\ampere}]")
plt.xlim((-130,130))
plt.ylim((-40,40))
plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"double_KL.pgf", bbox_inches="tight")
plt.savefig(PLOT+"double_KL.pdf", bbox_inches="tight")
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
plt.xlim((0.4,20))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"double_Te.pgf", bbox_inches="tight")
plt.savefig(PLOT+"double_Te.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, nes, yerr=dnes, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$n_e$ [\si{\per\metre\tothe{3}}]")
plt.xscale("log")
plt.xlim((0.4,20))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"double_ne.pgf", bbox_inches="tight")
plt.savefig(PLOT+"double_ne.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
aerr   = np.array([dalphas,dalphas])
aerr[0,aerr[0,:]>=alphas] = alphas[aerr[0,:]>=alphas]*0.99999999
plt.errorbar(ps, alphas, yerr=aerr, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\alpha")
plt.xscale("log")
plt.xlim((0.4,20))
plt.ylim((10**-13,5*10**-11))
plt.yscale("log")
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"double_alpha.pgf", bbox_inches="tight")
plt.savefig(PLOT+"double_alpha.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, debs, yerr=ddebs, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\lambda_\text{D}$ [\si{\micro\metre}]")
plt.xscale("log")
plt.xlim((0.4,20))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"double_deb.pgf", bbox_inches="tight")
plt.savefig(PLOT+"double_deb.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, wpes, yerr=dwpes, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\omega_{e}$ [\si{\giga\hertz}]")
plt.xscale("log")
plt.xlim((0.4,20))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"double_wpe.pgf", bbox_inches="tight")
plt.savefig(PLOT+"double_wpe.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
plt.errorbar(ps, wpis, yerr=dwpis, xerr=perr, fmt='x')
plt.xlabel(r"$p$ [\si{\milli\bar}]")
plt.ylabel(r"$\omega_{i}$ [\si{\mega\hertz}]")
plt.xscale("log")
plt.xlim((0.4,20))
#plt.ylim((-50,1050))
#plt.legend()
plt.tight_layout()
plt.savefig(PLOT+"double_wpi.pgf", bbox_inches="tight")
plt.savefig(PLOT+"double_wpi.pdf", bbox_inches="tight")
plt.close()

#   Radial distance
#######################################################################################
IEs      = [5, 10, 20]  # mA
p        = 3.           # mbar
drs      = [0., 1., 2., 3., 4., 5., 6., 7.]   # mm
datass   = [["radial_Ar_17.dat",
             "radial_Ar_18.dat",
             "radial_Ar_19.dat",
             "radial_Ar_20.dat",
             "radial_Ar_21.dat",
             "radial_Ar_22.dat",
             "radial_Ar_23.dat",
             "radial_Ar_24.dat"],
            ["radial_Ar_16.dat",
             "radial_Ar_15.dat",
             "radial_Ar_14.dat",
             "radial_Ar_13.dat",
             "radial_Ar_12.dat",
             "radial_Ar_11.dat",
             "radial_Ar_10.dat",
             "radial_Ar_9.dat"],
            ["radial_Ar_1.dat",
             "radial_Ar_2.dat",
             "radial_Ar_3.dat",
             "radial_Ar_4.dat",
             "radial_Ar_5.dat",
             "radial_Ar_6.dat",
             "radial_Ar_7.dat",
             "radial_Ar_8.dat"]]
U_maxs   = [-60., -45., -45.]
aTes     = []
adTes    = []
aIisats  = []
adIisats = []
anes     = []
adnes    = []
aalphas  = []
adalphas = []
adebs    = []
addebs   = []
awpes    = []
adwpes   = []
awpis    = []
adwpis   = []
for IE, datas, U_max, A in zip(IEs, datass, U_maxs, Alphabet[2:]): 
    set_pre(A)
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
    for data, a, dr in zip(datas, alphabet, drs):
        U, I    = getdata_cut(data)
        # Exponential fit to the curves
        popt, perr  = fit_single(U, I, dI, U_max=U_max)
        si("p"+a,p,dp,"mbar",1)
        si("Te"+a,popt[0],perr[0],"eV",2)
        si("Iisat"+a,np.abs(popt[1]),perr[1],"muA",2)
        si("A"+a,popt[2]/1000.,perr[2]/1000.,"mA",2)
        Tes.append(popt[0])
        dTes.append(perr[0])
        Iisats.append(popt[1])
        dIisats.append(perr[1])
        siline("popt"+a,[dr, None, popt[0], perr[0], popt[1], perr[1], popt[2]/1000., perr[2]/1000.],[0,2,2,0])
        
        # Calculate electron density
        ne  = -popt[1]*np.sqrt(mAr/con.e/popt[0])/0.61/con.e/S*10**-6   #1_m3
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
        deb     = np.sqrt(con.epsilon_0*popt[0]/ne/con.e)*10**6    # mum
        ddeb    = np.abs(0.5*deb*perr[0]/popt[0]) + np.abs(0.5*deb*dne/ne)
        si("deb"+a,deb,ddeb,"mum",0)
        debs.append(deb)
        ddebs.append(ddeb)
        wpe     = np.sqrt(ne*con.e**2/con.epsilon_0/con.m_e)*10**-9   # GHz
        dwpe    = np.abs(0.5*wpe*dne/ne)
        si("wpe"+a,wpe,dwpe,"GHz",2)
        wpes.append(wpe)
        dwpes.append(dwpe)
        wpi     = np.sqrt(ne*con.e**2/con.epsilon_0/mAr)*10**-6   # MHz
        dwpi    = np.abs(0.5*wpi*dne/ne)
        si("wpi"+a,wpi,dwpi,"MHz",2)
        wpis.append(wpi)
        dwpis.append(dwpi)
        
        # Plot of all characteristic curves
        x           = np.linspace(-150,-20,1000)
        p1  = plt.plot(x, I_single(x, *popt), zorder=3)
        plt.errorbar(U, I, xerr=dU, yerr=dI, 
            fmt='.', 
            label=r"$\Delta r=\SI{{{:.0f}}}{{\milli\metre}}$".format(dr), 
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
    plt.savefig(PLOT+"radial_KL_I{:.0f}.pgf".format(IE), bbox_inches="tight")
    plt.savefig(PLOT+"radial_KL_I{:.0f}.pdf".format(IE), bbox_inches="tight")
    plt.close()

    alphas  = np.array(alphas)
    dalphas = np.array(dalphas)
    
    aTes.append(Tes)
    adTes.append(dTes)
    aIisats.append(Iisats)
    adIisats.append(dIisats)
    anes.append(nes)
    adnes.append(dnes)
    aalphas.append(alphas)
    adalphas.append(dalphas)
    adebs.append(debs)
    addebs.append(ddebs)
    awpes.append(wpes)
    adwpes.append(dwpes)
    awpis.append(wpis)
    adwpis.append(dwpis)

f   = plt.figure()
for k in [0,1,2]:
    plt.errorbar(drs, aTes[k], yerr=adTes[k], fmt='x', label=r"$I_\text{{E}}=\SI{{{:.0f}}}{{\milli\ampere}}$".format(IEs[k]))
plt.xlabel(r"$\Delta r$ [\si{\milli\metre}]")
plt.ylabel(r"$T_e$ [\si{\electronvolt}]")
#plt.xscale("log")
plt.xlim((-1,10))
#plt.ylim((-50,1050))
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(PLOT+"radial_Te.pgf", bbox_inches="tight")
plt.savefig(PLOT+"radial_Te.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
for k in [0,1,2]:
    plt.errorbar(drs, anes[k], yerr=adnes[k], fmt='x', label=r"$I_\text{{E}}=\SI{{{:.0f}}}{{\milli\ampere}}$".format(IEs[k]))
plt.xlabel(r"$\Delta r$ [\si{\milli\metre}]")
plt.ylabel(r"$n_e$ [\si{\per\metre\tothe{3}}]")
#plt.xscale("log")
plt.xlim((-1,10))
#plt.ylim((-50,1050))
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(PLOT+"radial_ne.pgf", bbox_inches="tight")
plt.savefig(PLOT+"radial_ne.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
for k in [0,1,2]:
    aerr   = np.array([adalphas[k],adalphas[k]])
    aerr[0,aerr[0,:]>=aalphas[k]] = aalphas[k][aerr[0,:]>=aalphas[k]]*0.99999999
    plt.errorbar(drs, aalphas[k], yerr=aerr, fmt='x', label=r"$I_\text{{E}}=\SI{{{:.0f}}}{{\milli\ampere}}$".format(IEs[k]))
plt.xlabel(r"$\Delta r$ [\si{\milli\metre}]")
plt.ylabel(r"$\alpha")
#plt.xscale("log")
plt.xlim((-1,10))
#plt.ylim((10**-7,10**-5))
#plt.ylim((0,225*10**-8))
#plt.yscale("log")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(PLOT+"radial_alpha.pgf", bbox_inches="tight")
plt.savefig(PLOT+"radial_alpha.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
for k in [0,1,2]:
    plt.errorbar(drs, adebs[k], yerr=addebs[k], fmt='x', label=r"$I_\text{{E}}=\SI{{{:.0f}}}{{\milli\ampere}}$".format(IEs[k]))
plt.xlabel(r"$\Delta r$ [\si{\milli\metre}]")
plt.ylabel(r"$\lambda_\text{D}$ [\si{\micro\metre}]")
#plt.xscale("log")
plt.xlim((-1,10))
#plt.ylim((-50,1050))
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(PLOT+"radial_deb.pgf", bbox_inches="tight")
plt.savefig(PLOT+"radial_deb.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
for k in [0,1,2]:
    plt.errorbar(drs, awpes[k], yerr=adwpes[k], fmt='x', label=r"$I_\text{{E}}=\SI{{{:.0f}}}{{\milli\ampere}}$".format(IEs[k]))
plt.xlabel(r"$\Delta r$ [\si{\milli\metre}]")
plt.ylabel(r"$\omega_{e}$ [\si{\giga\hertz}]")
#plt.xscale("log")
plt.xlim((-1,10))
#plt.ylim((-50,1050))
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(PLOT+"radial_wpe.pgf", bbox_inches="tight")
plt.savefig(PLOT+"radial_wpe.pdf", bbox_inches="tight")
plt.close()

f   = plt.figure()
for k in [0,1,2]:
    plt.errorbar(drs, awpis[k], yerr=adwpis[k], fmt='x', label=r"$I_\text{{E}}=\SI{{{:.0f}}}{{\milli\ampere}}$".format(IEs[k]))
plt.xlabel(r"$\Delta r$ [\si{\milli\metre}]")
plt.ylabel(r"$\omega_{i}$ [\si{\mega\hertz}]")
#plt.xscale("log")
plt.xlim((-1,10))
#plt.ylim((-50,1050))
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(PLOT+"radial_wpi.pgf", bbox_inches="tight")
plt.savefig(PLOT+"radial_wpi.pdf", bbox_inches="tight")
plt.close()
