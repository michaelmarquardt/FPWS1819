import numpy as np
import matplotlib.pyplot as plt
# Kerr Winkel bei saettigung der magnetisierung 
winkel=np.genfromtxt("Winkel")
plt.figure()
plt.plot(winkel[:,0], winkel[:,1], 'x')
plt.ylabel(r'$\mathrm{Kerr-Winkel}\   [{}^\circ] $')
plt.xlabel(r'$ \mathrm{Temperatur}\   [{}^\circ \mathrm{C}]$' )
#plt.legend()
plt.savefig("Kerr_Winkel.pdf")
plt.close()
