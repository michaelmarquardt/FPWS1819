import matplotlib.pyplot as plt
import numpy as np
rawdata_1_1=np.genfromtxt("Oszilloskop_1_1.dat", delimiter=",")
rawdata_1_2=np.genfromtxt("Oszilloskop_1_2.dat", delimiter=",")
rawdata_2_1=np.genfromtxt("Oszilloskop_1_1.dat", delimiter=",")
rawdata_2_2=np.genfromtxt("Oszilloskop_1_2.dat", delimiter=",")
rawdata_3_1=np.genfromtxt("Oszilloskop_1_1.dat", delimiter=",")
rawdata_3_2=np.genfromtxt("Oszilloskop_1_2.dat", delimiter=",")


plt.figure()
plt.plot(rawdata_1_1[:,0], rawdata_1_1[:,1], label="Langmuir Sonde")
plt.plot(rawdata_1_2[:,0], rawdata_1_2[:,1], label="Gitter")
plt.xlabel(r'$\mathrm{Zeit}\ t\  [\mu \mathrm{s}] $')
plt.ylabel("Amplitude [a.u.]")
plt.legend()
#plt.show()
plt.savefig("oszi_1.pdf")
plt.close()

plt.figure()
plt.plot(rawdata_2_1[:,0], rawdata_2_1[:,1], label="Langmuir Sonde")
plt.plot(rawdata_2_2[:,0], rawdata_2_2[:,1], label="Gitter")
plt.xlabel(r'$\mathrm{Zeit}\ t\  [\mu \mathrm{s}] $')
plt.ylabel("Amplitude [a.u.]")
plt.legend()
#plt.show()
plt.savefig("oszi_2.pdf")
plt.close()


plt.figure()
plt.plot(rawdata_3_1[:,0], rawdata_3_1[:,1], label="Langmuir Sonde")
plt.plot(rawdata_3_2[:,0], rawdata_3_2[:,1], label="Gitter")
plt.xlabel(r'$\mathrm{Zeit}\ t\  [\mu \mathrm{s}] $')
plt.ylabel("Amplitude [a.u.]")
plt.legend()
#plt.show()
plt.savefig("oszi_3.pdf")
plt.close()




