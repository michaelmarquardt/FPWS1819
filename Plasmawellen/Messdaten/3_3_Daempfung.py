import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
from scipy.optimize import curve_fit
# in .dat Position in [cm] und paeak to peak in [mV]
# peak tp peak in amplitude umrechnen
# peak to peak fehler: 10 mV
Delta_Amplitude = 5
Delta_ort = 0.1 # ort fehler 0.1 cm
rawdata_1 = np.genfromtxt("3_3_Daempfung_1.dat", delimiter=',')
rawdata_2 = np.genfromtxt("3_3_Daempfung_2.dat", delimiter=',')
rawdata_3 = np.genfromtxt("3_3_Daempfung_3.dat", delimiter=',')
Amplitude_1=0.5*rawdata_1[:,0]
Amplitude_2=0.5*rawdata_2[:,0]
Amplitude_3=0.5*rawdata_3[:,0]
ort_1 = rawdata_1[:,1]
ort_2 = rawdata_2[:,1]
ort_3 = rawdata_3[:,1]
plt.yscale("log")
plt.errorbar(ort_1, Amplitude_1, yerr=Delta_Amplitude,  fmt='r.', capsize=5, label="Druck: 5.5 mm")
plt.errorbar(ort_2, Amplitude_2, yerr=Delta_Amplitude,  fmt='b.', capsize=5, label="Druck: 7.3 mm")
plt.errorbar(ort_3, Amplitude_3, yerr=Delta_Amplitude,  fmt='y.', capsize=5, label="Druck: 6.5 mm")
def fit(x,a,b):
    return a*np.exp(-x/b)
init_vals = [130, 10]
popt_1, pcov_1 = curve_fit(fit, ort_1, Amplitude_1, p0=init_vals)
popt_2, pcov_2 = curve_fit(fit, ort_2, Amplitude_2, p0=init_vals)
popt_3, pcov_3 = curve_fit(fit, ort_3, Amplitude_3, p0=init_vals)
print(popt_1)
print(pcov_1)
print(popt_2)
print(pcov_2)
print(popt_3)
print(pcov_3)
plt.plot(ort_1,fit(ort_1,*popt_1), 'r--', label="Fit: Druck 5.5 mm")
plt.plot(ort_2,fit(ort_2,*popt_2), 'b--', label="Fit: Druck 7.3 mm")
plt.plot(ort_3,fit(ort_3,*popt_3), 'y--', label="Fit: Druck 6.5 mm")
plt.ylabel("Amplitude [mV]")
plt.xlabel("Position [cm]")
plt.legend()
plt.savefig("3_3_Daempfung.pdf")
plt.show()
plt.close()

