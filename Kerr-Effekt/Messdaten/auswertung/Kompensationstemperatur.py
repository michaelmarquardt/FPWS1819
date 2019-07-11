import matplotlib.pyplot as plt
import scipy as sc
from scipy.optimize import curve_fit
import numpy as np

rawdata=np.genfromtxt("Kompensationstemperatur", delimiter=",")
rawdata[:,0]+=273.15
#Fitten: 
def fit(T,H_O,T_comp):
	return (H_O*T)/(T-T_comp)

popt,pcov = curve_fit(fit,rawdata[:,0],rawdata[:,3], p0 = [344,300])

print(popt)
print(pcov)



plt.figure()
#plt.plot(rawdata[:,0], rawdata[:,1], '.', label="links")
#plt.plot(rawdata[:,0], rawdata[:,2], ',', label="rechts")
plt.plot(rawdata[:,0], rawdata[:,3], 'x', label="Messwerte")
plt.plot(rawdata[:,0], fit(rawdata[:,0],*popt), label = "Fitfunktion ")
plt.ylabel(r'$\mathrm{Koerzitivfeldstärke}\   [\mathrm{mT}] $')
plt.xlabel(r'$ \mathrm{Temperatur}\   [\mathrm{K}]$' )
#plt.errorbar(rawdata_1[:,0], rawdata_1[:,1], xerr=0.2, yerr=0.02, fmt='.', capsize=3, label="2. Filament ist aus")
#plt.errorbar(rawdata_2[:,0], rawdata_2[:,1], xerr=0.2, yerr=0.02, fmt='.', capsize=3, label="2. Filament ist an")
plt.legend()
plt.savefig("Kompensationstemperatur.pdf")
#plt.show()
plt.close()



