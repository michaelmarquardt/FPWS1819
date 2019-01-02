import matplotlib.pyplot as plt
import numpy as np
# plasmafrequanz: w_pe= sqrt(e/(epsilon_0 * me) sqrt(n_e)
# Naeherung w_pe= 56.5 sqrt(n_e)
# beachte w=2 pi f
#Ungenauigkeit der Messung:
df=5*10**6


# Spannung variieren
rawdata_1=np.genfromtxt("3_2_Spannung.dat", delimiter=",")
n_1=((2*np.pi*rawdata_1[:,1]*10**6)/(56.5))**2
Delta_n_1=((2*np.pi)/(56.5))**2*2*rawdata_1[:,1]*10**6*df
plt.figure()
#plt.plot(rawdata_1[:,0],n_1, 'x')
plt.errorbar(rawdata_1[:,0], n_1, xerr=0.2, yerr=Delta_n_1, fmt='.', capsize=3,)
plt.xlabel(r'$\mathrm{Entladungsspannung}\ I_{\mathrm{D}}\  [\mathrm{A}] $')
plt.ylabel(r'$ \mathrm{Plasmadichte}\  n\  [\mathrm{m}^{-3}] $' )
plt.savefig("3_2_Spannung.pdf")
#plt.show()
plt.close()
# Strom variieren
rawdata_2=np.genfromtxt("3_2_Strom.dat", delimiter=",")
n_2=((2*np.pi*rawdata_2[:,1]*10**6)/(56.5))**2
Delta_n_2=((2*np.pi)/(56.5))**2*2*rawdata_2[:,1]*10**6*df
plt.figure()
plt.errorbar(rawdata_2[:,0], n_2, xerr=0.02, yerr=Delta_n_2, fmt='.', capsize=3,)
#plt.plot(rawdata_2[:,0],n_2, 'x')
plt.xlabel(r'$\mathrm{Heizstrom}\ I_{\mathrm{H}}\  [\mathrm{A}] $')
plt.ylabel(r'$ \mathrm{Plasmadichte}\  n\  [\mathrm{m}^{-3}] $' )
plt.savefig("3_2_Strom.pdf")
#plt.show()
plt.close()
# Druck variieren
rawdata_3=np.genfromtxt("3_2_Druck.dat", delimiter=",")
n_3=((2*np.pi*rawdata_3[:,1]*10**6)/(56.5))**2
Delta_n_3=((2*np.pi)/(56.5))**2*2*rawdata_3[:,1]*10**6*df
plt.figure()
plt.errorbar(rawdata_3[:,0], n_3, xerr=0.2, yerr=Delta_n_3, fmt='.', capsize=3,)
#plt.plot(rawdata_3[:,0],n_3, 'x')
plt.xlabel(r'$\mathrm{Druck }\ p\   [\mathrm{mm}] $')
plt.ylabel(r'$ \mathrm{Plasmadichte}\  n\  [\mathrm{m}^{-3}] $' )
plt.savefig("3_2_Druck.pdf")
#plt.show()
plt.close()
