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

# Umrechnung der Winkel in die Intensitaeten
############################################
intens  = angles.copy()
ang, I  = np.loadtxt(DAT+"umrechnung.txt",unpack=True)
l       = 0
for k in range(len(angles)):
    if angles[k]%10 != 0:
        l   += 1
        intens[k]   = (I[k-l]+I[k-l+1])*0.5
    else:
        intens[k]   = I[k-l]

plt.figure("umrechnung")
plt.plot(angles, intens, "x-")
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



"""   
    # plot
    ts  = np.linspace(-500,500,10000)
    plt.plot(t,g,"x",label="Messdaten")
    plt.plot(ts,g_exp(ts,popt[0],popt[1],popt[2]),label="Fit")
    plt.ylim((0,2))
    plt.xlabel(r"$\tau$ [\si{\nano\second}]")
    plt.ylabel(r"$g^{(2)}$")
    if not k == 1:
        plt.legend(loc="lower left")
    plt.savefig("../plots/g_{}.pgf".format(k+1))
    plt.savefig("../plots/g_{}.pdf".format(k+1))
    plt.close()

   """
