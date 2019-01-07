import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
from scipy.optimize import curve_fit
# in .dat: frequenz in [kHz] und Wellenlaenge in [cm], in diesen Einheiten auch Rechnen
# das umrechen in w und k
# BEACHTE: Fitparameter muss mit 10 multipliziert werden, damit Geschwindigkeit rauskommt, wegen Wahl der Einheiten
# fehler fuer Wellenlaengenmessung
Delta_lambda = 0.1   # Einheit: cm
def fit(x,a,b):
    return a*x+b
rawdata_1 = np.genfromtxt("3_3_Dispersion_1.dat", delimiter=',')
rawdata_2 = np.genfromtxt("3_3_Dispersion_2.dat", delimiter=',')
rawdata_3 = np.genfromtxt("3_3_Dispersion_3.dat", delimiter=',')
w_1 = 2*np.pi* rawdata_1[:,1]      
k_1 = 2*np.pi/(rawdata_1[:,0])
w_2 = 2*np.pi* rawdata_2[:,1]      
k_2 = 2*np.pi/(rawdata_2[:,0])
w_3 = 2*np.pi* rawdata_3[:,1]      
k_3 = 2*np.pi/(rawdata_3[:,0])
Delta_K_1=2*np.pi/(rawdata_1[:,0]**2) * Delta_lambda
Delta_K_2=2*np.pi/(rawdata_2[:,0]**2) * Delta_lambda
Delta_K_3=2*np.pi/(rawdata_3[:,0]**2) * Delta_lambda
#Delta_K_1=Delta_K_2=Delta_K_3=0.1
popt_1, pcov_1 = curve_fit(fit, k_1,w_1)
popt_2, pcov_2 = curve_fit(fit, k_2,w_2)
popt_3, pcov_3 = curve_fit(fit, k_3,w_3)
print(popt_1)
print(pcov_1)
print(popt_2)
print(pcov_2)
print(popt_3)
print(pcov_3)
plt.figure()
plt.errorbar(k_1, w_1, xerr=Delta_K_1,  fmt='r.', capsize=5, label="Druck: 5.5 mm")
plt.errorbar(k_2, w_2, xerr=Delta_K_2,  fmt='b.', capsize=5, label="Druck: 7.3 mm")
plt.errorbar(k_3, w_3, xerr=Delta_K_3,  fmt='g.', capsize=5, label="Druck: 6.5 mm")
plt.xlabel(r'$\mathrm{Wellenvektor}\ k\   [\mathrm{cm}^{-1}] $')
plt.ylabel(r'$\mathrm{Frequenz}\ \omega\ [\mathrm{kHz}] $')
plt.plot(k_1,fit(k_1,*popt_1), 'r--', label="Fit: Druck 5.5 mm")
plt.plot(k_2,fit(k_2,*popt_2), 'b--', label="Fit: Druck 7.3 mm")
plt.plot(k_3,fit(k_3,*popt_3), 'g--', label="Fit: Druck 6.5 mm")
plt.legend()
plt.savefig("3_3_Dispersion.pdf")
plt.show()
plt.close()

