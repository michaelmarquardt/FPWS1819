# This file calculates the correlation function
import numpy as np
import subprocess
from scipy.optimize import curve_fit
from export_data import *
texfile("../report/values.tex")
set_dig(3)
set_pre("")
usepgf()
import matplotlib.pyplot as plt

#   Calculation
#######################################################################################
Theory  = \
"""
g(t)    = (CN(t) - (1-rho**2))/rho**2
CN(t)   = c(t)/(N1*N2*omega*T)
rho     = S/(S+B)
S       = N1+N2-B
c(t)    = data (total counts signal 1 - signal 2)
N1      = count rate signal 1
N2      = count rate signal 2
omega   = bin width (0.5 ns)
T       = measurement time
S       = count rate signal
B       = count rate background
""" 
#print(Theory)

# overall constants
omega   = 0.5*10**-9    # s
B       = 1100.         # counts/s

# NV center dependend properties
N1      = np.array([11000.,12000.,120000.])   # counts/s
N2      = np.array([20000.,21000.,200000.])   # counts/s
T       = np.array([6128.20,6525.34,460.99])  # s
f       = N2/N1

# preparing fit functions
def g_exp(t,tau,a,h):
    return h - (h-a)*np.exp(-np.abs(t)/tau)

# calculations
subprocess.call(["mkdir","-p","../plots"])
for k in (0,1):
    print("measurement {}".format(k+1))
    t, c    = np.loadtxt("../Photonenstatistik_Messdaten/messung_{}.dat".format(k+1),unpack=True)
    
    # optimization of the count rates N1 and N2:
    popt    = [0.,0.,2.]
    while abs(popt[2]-1.) > 1e-7:
        # calculations
        S       = N1[k] + N2[k] - B
        rho     = S/(S+B)
        CN      = c/(N1[k]*N2[k]*omega*T[k])
        g       = (CN-(1-rho**2))/rho**2
        
        # fit (separate for + and -)
        p0  = [10.,g.min(),g.mean()]
        popt,  pcov = curve_fit(g_exp,t,g,p0=p0)

        # rescale N1 and N2
        N1[k]  *= np.sqrt(popt[2])
        N2[k]  *= np.sqrt(popt[2])
    
    print("   T         = {}".format(T[k]))
    print("   N1        = {}".format(N1[k]))
    print("   N2        = {}".format(N2[k]))
    print("   f         = {}".format(f[k]))
    
    # calculations
    S       = N1[k] + N2[k] - B
    print("    S        = {:n} counts/s".format(S))
    rho     = S/(S+B)
    print("    rho      = {:.3f}".format(rho))
    CN      = c/(N1[k]*N2[k]*omega*T[k])
    g       = (CN-(1-rho**2))/rho**2
    
    # fit (separate for + and -)
    p0  = [10.,g.min(),g.mean()]
    popt,  pcov = curve_fit(g_exp,t,g,p0=p0)
    perr        = np.sqrt(np.diag(pcov))
    print("    popt     = {}".format(popt))
    print("    perr     = {}".format(perr))    
    
    # plot
    ts  = np.linspace(-500,500,10000)
    plt.plot(t,g,"x",label="Messdaten")
    plt.plot(ts,g_exp(ts,popt[0],popt[1],popt[2]),label="Fit")
    plt.ylim((0,2))
    plt.xlabel(r"Zeitdifferenz $t$ [\si{\nano\second}]")
    plt.ylabel(r"Korrelationsfunktion $g$")
    plt.legend(loc="upper right")
    plt.savefig("../plots/g_{}.pgf".format(k+1))
    plt.savefig("../plots/g_{}.pdf".format(k+1))
    plt.close()

k   = 2
print("measurement {}".format(k+1))
t, c    = np.loadtxt("../Photonenstatistik_Messdaten/messung_{}.dat".format(k+1),unpack=True)

# optimization of the count rates N1 and N2:
g   = c
while g.mean()-1. > 1e-7:
    # calculations
    S       = N1[k] + N2[k] - B
    rho     = S/(S+B)
    CN      = c/(N1[k]*N2[k]*omega*T[k])
    g       = (CN-(1-rho**2))/rho**2

    # rescale N1 and N2
    N1[k]  *= np.sqrt(g.mean())
    N2[k]  *= np.sqrt(g.mean())

print("   T         = {}".format(T[k]))
print("   N1        = {}".format(N1[k]))
print("   N2        = {}".format(N2[k]))
print("   f         = {}".format(f[k]))

# calculations
S       = N1[k] + N2[k] - B
print("    S        = {:n} counts/s".format(S))
rho     = S/(S+B)
print("    rho      = {:.3f}".format(rho))
CN      = c/(N1[k]*N2[k]*omega*T[k])
#CN     /= CN.mean()
g       = (CN-(1-rho**2))/rho**2
#g      /= g.mean()

# plot
tsp  = np.linspace(0,500,10000)
tsm  = np.linspace(-500,0,10000)
plt.plot(t,g,"x")
plt.ylim((0,2))
plt.xlabel(r"Zeitdifferenz $t$ [\si{\nano\second}]")
plt.ylabel(r"Korrelationsfunktion $g$")
plt.savefig("../plots/g_{}.pgf".format(k+1))
plt.savefig("../plots/g_{}.pdf".format(k+1))
plt.close()

   
