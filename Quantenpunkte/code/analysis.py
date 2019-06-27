# This file calculates the correlation function
import numpy as np
import subprocess
from scipy.optimize import curve_fit
from export_data import *
texfile("../report/values.tex")
set_dig(3)
set_pre("")
usepgf(size=(1.1,2.3))
import matplotlib.pyplot as plt

# Create folders and files
##########################
subprocess.call(["mkdir","-p","../plots"])
texfile("../report/values.tex")
PLOT    = "../plots/"
DAT     = "../messdaten/"

# Parameter
###########
actime  = {
    "deckenlampe"   : 0.1,
    "neon"          : 0.1,
    "na"            : 0.1,
    "InAs"          : 1.,
    "InP"           : 1.
    }
naum    = [0, 20, 40, 60, 80, 100, 150, 200, 250, 300, 400, 500, 550, 600]
angles  = [0, 10, 20, 25, 30, 35, 40, 45, 50, 60, 70]
temps   = [10, 30, 50, 70, 90, 110, 130, 150, 170]

# Umrechnung der Winkel in die powersitaeten
############################################
powers  = angles.copy()
ang, I  = np.loadtxt(DAT+"umrechnung.txt",unpack=True)
l       = 0
for k in range(len(angles)):
    if angles[k]%10 != 0:
        l   += 1
        powers[k]   = (I[k-l]+I[k-l+1])*0.5
    else:
        powers[k]   = I[k-l]

plt.figure("umrechnung")
plt.plot(angles, powers, "x-")
plt.xlabel("Polarisator Winkel [\si{\degree}]")
plt.ylabel("$I$ [mW]")
plt.savefig(PLOT+"umrechnung.pgf")
plt.savefig(PLOT+"umrechnung.pdf")


# Deckenlampe
#############
lam, I  = np.loadtxt(DAT+"deckenlampe.asc",unpack=True)
I       /= actime["deckenlampe"]
plt.figure("deckenlampe")
plt.plot(lam,I,"-")
plt.xlabel("$\lambda$ [nm]")
plt.ylabel("$I$ [counts/s]")
plt.savefig(PLOT+"deckenlampe.pgf")
plt.savefig(PLOT+"deckenlampe.pdf")
plt.close()

# Neon Dampflampe
#################
lam, I  = np.loadtxt(DAT+"neon.asc",unpack=True)
I       /= actime["neon"]
plt.figure("neon")
plt.plot(lam,I,"-")
plt.xlabel("$\lambda$ [nm]")
plt.ylabel("$I$ [counts/s]")
plt.savefig(PLOT+"neon.pgf")
plt.savefig(PLOT+"neon.pdf")
plt.close()

# Natrium Dampflampe
####################
f1  = plt.figure("na")
for um in naum:
    lam, I  = np.loadtxt(DAT+"na{:d}um.asc".format(um),unpack=True)
    I       /= actime["na"]
    plt.plot(lam, I, "-", label=r"\SI{{{:d}}}{{\micro\metre}}".format(um))
plt.xlim((585,600))
#plt.ylim((0,600000))
plt.xlabel("$\lambda$ [nm]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"na.pgf")
plt.savefig(PLOT+"na.pdf")
plt.close()

# InAs Quantenpunkte
####################
# Bereite fit vor
def multigauss(x, a0, a1, a2, a3, mu0, mu1, mu2, mu3, s0, s1, s2, s3):
    """
    Up to 4 gaussian like
    f(x)    = a*exp(-(x-mu)**2/(2*sigma**2))
    """
    f   = np.zeros_like(x)
    for ai, mui, si in zip([a0, a1, a2, a3], [mu0, mu1, mu2, mu3], [s0, s1, s2, s3]):
        f   += ai*np.exp(-(x-mui)**2/(2.*si**2))
    return f

def fitto(lam,I):
    # Fit
    lb  = np.where(lam>850)[0][0]
    mu0 = [921.,902.,884.,866.]
    s0  = [5.,5.,5.,5.]
    a0  = []
    for mu in mu0:
        a0.append(I[np.where(lam>=mu)[0][0]])
    p0  = a0+mu0+s0

    # Disable some gaussian using bounds
    bounds  = (np.full(12,-np.inf),np.full(12,np.inf))
    # Initial bounds
    bounds[0][0:4]  = [0.,0.,0.,0.]
    bounds[0][4:8]  = [910.,890.,875.,860.]
    bounds[1][4:8]  = [940.,920.,895.,880.]
    sucess  = False
    num     = 3
    while sucess == False:
        try:
            popt, pcov  = curve_fit(multigauss, lam[lb:], I[lb:], p0=p0, bounds=bounds)
            sucess      = True
        except:
            if num < 0:
                raise IndexError("Fit function does not converge for 0 gaussian.")
            bounds[0][num]      = 0.
            bounds[1][num]      = 0.1
            p0[num]             = 0.0
            bounds[0][num+4]    = mu0[num]*0.99
            bounds[1][num+4]    = mu0[num]*1.01
            bounds[0][num+8]    = s0[num]*0.99
            bounds[1][num+8]    = s0[num]*1.01
            print("Fit function does not converge: Set Peak {:d} to 0.".format(num))
            num -= 1
    perr        = np.sqrt(np.diag(pcov))
    for k in range(4):
        if perr[k]  > popt[k]:
            perr[k] = popt[k]
    return popt, perr

# Prepare parameters for collection
shell   = np.zeros((4,len(angles)))
dshell  = np.zeros((4,len(angles)))

# Execute fit
plt.figure("InAs")
for angle, power, k in zip(angles[::-1], powers[::-1], range(len(angles))):
    print("angle = {:d}".format(angle))
    lam, I  = np.loadtxt(DAT+"InAs{:d}deg.asc".format(angle),unpack=True)
    I       /= actime["InAs"]
    popt, perr  = fitto(lam,I)
    #print(popt)
    #print(perr)
    shell[:,k]  = popt[:4]
    dshell[:,k] = perr[:4]
    p1  = plt.plot(
        lam, 
        I, 
        "x", 
        zorder  = 0,
        alpha  = 0.2
        )
    plt.plot(
        lam, 
        multigauss(lam, *popt), 
        "-", 
        label   = r"$P=\SI{{{:.2f}}}{{\milli\watt}}$".format(power),
        color   = changecolor(p1[0].get_color(),rgb=(0,0,0)),
        zorder  = 1
        )
plt.xlabel("$\lambda$ [nm]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"InAs_power.pgf")
plt.savefig(PLOT+"InAs_power.pdf")
plt.close()

# Comparison plot
plt.figure("InAs_overpower")
plt.errorbar(powers, shell[0,:], dshell[0,:], fmt="x-", label="s-Schale")
plt.errorbar(powers, shell[1,:], dshell[1,:], fmt="x-", label="p-Schale")
plt.errorbar(powers, shell[2,:], dshell[2,:], fmt="x-", label="d-Schale")
plt.errorbar(powers, shell[3,:], dshell[3,:], fmt="x-", label="f-Schale")
plt.xlabel("$P$ [mW]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"InAs_overpower.pgf")
plt.savefig(PLOT+"InAs_overpower.pdf")
plt.close()

# Prepare parameters for collection
shellt  = np.zeros((4,len(temps)))
dshellt = np.zeros((4,len(temps)))

# Execute fit
plt.figure("InAs_temperature")
for temp, k in zip(temps, range(len(temps))):
    print("temp = {:d}".format(temp))
    lam, I  = np.loadtxt(DAT+"InAs{:d}K.asc".format(temp),unpack=True)
    I       /= actime["InAs"]
    popt, perr  = fitto(lam,I)
    #print(popt)
    #print(perr)
    shellt[:,k]     = popt[:4]
    dshellt[:,k]    = perr[:4]
    p1  = plt.plot(
        lam, 
        I, 
        "x", 
        zorder  = 0,
        alpha  = 0.2
        )
    plt.plot(
        lam, 
        multigauss(lam, *popt), 
        "-", 
        label   = r"$T=\SI{{{:d}}}{{\kelvin}}$".format(temp),
        color   = changecolor(p1[0].get_color(),rgb=(0,0,0)),
        zorder  = 1
        )
plt.xlabel("$\lambda$ [nm]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"InAs_temperature.pgf")
plt.savefig(PLOT+"InAs_temperature.pdf")
plt.close()

# Comparison plot
plt.figure("InAs_overpower")
plt.errorbar(temps, shellt[0,:], dshellt[0,:], fmt="x-", label="s-Schale")
plt.errorbar(temps, shellt[1,:], dshellt[1,:], fmt="x-", label="p-Schale")
plt.errorbar(temps, shellt[2,:], dshellt[2,:], fmt="x-", label="d-Schale")
plt.errorbar(temps, shellt[3,:], dshellt[3,:], fmt="x-", label="f-Schale")
plt.xlabel("$P$ [mW]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"InAs_overtemp.pgf")
plt.savefig(PLOT+"InAs_overtemp.pdf")
plt.close()

# InP Quantenpunkte
###################
# Bereite fit vor
def gauss(x, a, mu, s):
    return a*np.exp(-(x-mu)**2/(2.*s**2))

def fitto(lam,I):
    # Fit
    lb  = np.where(lam>780)[0][0]
    a0  = np.max(I[:lb])
    mu0 = lam[np.argmax(I[:lb])]
    s0  = 6.
    p0  = [a0,mu0,s0]
    popt, pcov  = curve_fit(gauss, lam[:lb], I[:lb], p0=p0)
    perr        = np.sqrt(np.diag(pcov))
    return popt, perr, lb

# Prepare parameters for collection
sInP   = np.zeros(len(angles))
dsInP  = np.zeros(len(angles))

# Execute fit
plt.figure("InP")
for angle, power, k in zip(angles[::-1], powers[::-1], range(len(angles))):
    print("angle = {:d}".format(angle))
    lam, I  = np.loadtxt(DAT+"InP{:d}deg.asc".format(angle),unpack=True)
    I       /= actime["InP"]
    popt, perr, lb  = fitto(lam,I)
    #print(popt)
    #print(perr)
    sInP[k]   = popt[0]
    dsInP[k]  = perr[0]
    p1  = plt.plot(
        lam, 
        I, 
        "x", 
        zorder  = 0,
        alpha  = 0.2
        )
    plt.plot(
        lam[:lb], 
        gauss(lam[:lb], *popt), 
        "-", 
        label   = r"$P=\SI{{{:.2f}}}{{\milli\watt}}$".format(power),
        color   = changecolor(p1[0].get_color(),rgb=(0,0,0)),
        zorder  = 1
        )
plt.xlabel("$\lambda$ [nm]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"InP_power.pgf")
plt.savefig(PLOT+"InP_power.pdf")
plt.close()

# Prepare parameters for collection
sInPt   = np.zeros(len(temps))
dsInPt  = np.zeros(len(temps))

# Execute fit
plt.figure("InP")
for temp, k in zip(temps, range(len(temps))):
    print("temp = {:d}".format(temp))
    lam, I  = np.loadtxt(DAT+"InP{:d}K.asc".format(temp),unpack=True)
    I       /= actime["InP"]
    popt, perr, lb  = fitto(lam,I)
    #print(popt)
    #print(perr)
    sInPt[k]    = popt[0]
    dsInPt[k]   = perr[0]
    p1  = plt.plot(
        lam, 
        I, 
        "x", 
        zorder  = 0,
        alpha  = 0.2
        )
    plt.plot(
        lam[:lb], 
        gauss(lam[:lb], *popt), 
        "-", 
        label   = r"$T=\SI{{{:d}}}{{\kelvin}}$".format(temp),
        color   = changecolor(p1[0].get_color(),rgb=(0,0,0)),
        zorder  = 1
        )
plt.xlabel("$\lambda$ [nm]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"InP_temperature.pgf")
plt.savefig(PLOT+"InP_temperature.pdf")
plt.close()

# Comparison plots
##################
plt.figure("overpower")
plt.errorbar(powers, shell[0,:], dshell[0,:], fmt="x-", label="s (InAs)")
plt.errorbar(powers, shell[1,:], dshell[1,:], fmt="x-", label="p (InAs)")
plt.errorbar(powers, shell[2,:], dshell[2,:], fmt="x-", label="d (InAs)")
plt.errorbar(powers, shell[3,:], dshell[3,:], fmt="x-", label="f (InAs)")
plt.errorbar(powers, sInP, dsInP, fmt="x-", label="s (InP)")
plt.xlabel("$P$ [mW]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"overpower.pgf")
plt.savefig(PLOT+"overpower.pdf")
plt.close()

plt.figure("overtemp")
plt.errorbar(temps, shellt[0,:], dshellt[0,:], fmt="x-", label="s (InAs)")
plt.errorbar(temps, shellt[1,:], dshellt[1,:], fmt="x-", label="p (InAs)")
plt.errorbar(temps, shellt[2,:], dshellt[2,:], fmt="x-", label="d (InAs)")
plt.errorbar(temps, shellt[3,:], dshellt[3,:], fmt="x-", label="f (InAs)")
plt.errorbar(temps, sInPt, dsInPt, fmt="x-", label="s (InP)")
plt.xlabel("$P$ [mW]")
plt.ylabel("$I$ [counts/s]")
plt.legend()
plt.savefig(PLOT+"overtemp.pgf")
plt.savefig(PLOT+"overtemp.pdf")
plt.close()
