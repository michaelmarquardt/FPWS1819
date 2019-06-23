import numpy as np
import formalism as f
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from numpy.linalg import eig
np.set_printoptions(linewidth=120)

"""
The values are impoerted from ../messwerte/.
    psi_minus psi_plus tomography
psi_plus:
    alpha_(prime)_beta_(prime).txt
alpha_beta.txt: (coincidence count rates at the detector)
    alpha_beta
    alpha_beta_or
    alpha_or_beta
    alpha_or_beta_or
"""

texfile = open("../report/values.tex","w")
DATA    = "../messwerte/"
angles  = ["alpha_beta.txt",
           "alpha_beta_prime.txt",
           "alpha_prime_beta.txt",
           "alpha_prime_beta_prime.txt"]
tables  = ["plus","minus","tomography"]

alphas  = [22.5,22.5,112.5,112.5,22.5,22.5,112.5,112.5,-22.5,-22.5,67.5,67.5,-22.5,-22.5,67.5,67.5]
betas   = [0,90,0,90,-45,45,-45,45,0,90,0,90,-45,45,-45,45]


for state, table in zip((DATA+"psi_plus/",DATA+"psi_minus/"),tables): 
    # Prepare arrays
    big_array   = np.zeros((16,5))
    E           = np.zeros(4)
    dE          = np.zeros(4)
    
    # Loop over array compositions
    for k in range(4):
        # Import data
        data    = np.loadtxt(state+angles[k],delimiter=",")
        
        # Save to big array
        big_array[4*k:4*k+4,:]  = data.copy()
        
        # Calculate mean coincidence count rates
        C   = np.mean(data,axis=1)
        dC  = np.std(data,axis=1)
        
        # Calculate expectation values
        E[k]    = (C[0]-C[1]-C[2]+C[3])/np.sum(C)
        dE[k]   = np.sum(dC)/np.sum(C)*(1.+np.abs(E[k]))
        
        print(C)
        print(dC)
        print(E)
        print(dE)
        
    # Calculate CHSH inequality
    S   = E[0]-E[1]+E[2]+E[3]
    dS  = np.sum(dE)
    
    print("S    = {:.3f} +- {:.3f}".format(S,dS))
    
    texfile.write("\\newcommand{{\\bigtable{:s}}}{{\n".format(table))
    texfile.write("\\begin{tabular}{S[table-format=+3.1]S[table-format=+2]|S[table-format=3]S[table-format=3]S[table-format=3]S[table-format=3]S[table-format=3]}\n")
    texfile.write("{$\\alpha$ [\si{\degree}]}&{$\\beta$ [\si{\degree}]}&\multicolumn{5}{c}{$C(\\alpha,\\beta)$ [counts/s]}\\\\\\hline\n")
    for l in range(16):
        texfile.write("{:.1f}&{:d}&{:.0f}&{:.0f}&{:.0f}&{:.0f}&{:.0f}\\\\".format(
            alphas[l],
            betas[l],
            big_array[l,0],
            big_array[l,1],
            big_array[l,2],
            big_array[l,3],
            big_array[l,4]))
        if (l+1)%4==0:
            texfile.write(r"\hline")
        texfile.write("\n")
    texfile.write("\\end{tabular}}\n\n")
    
#   Tomography
##############
"""
Calculate the density matrix rho.
Then use formalism to calculate S with this rho.
"""
# state key: Photon1 Photon2 pol1 lamb1 pol2 lamb2
state_key   = [["H", "H",   0,   0,   0,   0],
               ["H", "V",   0,   0,  90,   0],
               ["V", "V",  90,   0,  90,   0],
               ["V", "H",  90,   0,   0,   0],
               ["R", "H",  90,  45,   0,   0],
               ["R", "V",  90,  40,  90,   0],
               ["D", "V",  45,  45,  90,   0],
               ["D", "H",  45,  45,   0,   0],
               ["D", "R",  45,  45,  45,   0],
               ["D", "D",  45,  45,  45,  45],
               ["R", "D",  45,   0,  45,  45],
               ["H", "D",   0,   0,  45,  45],
               ["V", "D",  90,   0,  45,  45],
               ["V", "L",  90,   0,   0,  45],
               ["H", "L",   0,   0,   0,  45],
               ["R", "L",   0, -45,   0,  45]]
# The state is a superposition of all those states.
# Their probability P can be calculated from the count rates C in the datafile.
data    = np.loadtxt(DATA+"tomography/tomography.txt",delimiter=",")

# Big table
texfile.write("\\newcommand{\\bigtabletom}{\n")
texfile.write("\\begin{tabular}{S[table-format=2]|ll|S[table-format=+2]S[table-format=+2]|S[table-format=2]S[table-format=2]||S[table-format=3]S[table-format=3]S[table-format=3]S[table-format=3]S[table-format=3]}\n")
texfile.write("&{S1}&{S2}&{P1 [\si{\degree}]}&{L1 [\si{\degree}]}&{P2 [\si{\degree}]}&{L2 [\si{\degree}]}&\multicolumn{5}{c}{$C(\\alpha,\\beta)$ [counts/s]}\\\\\\hline\n")
for l in range(16):
    texfile.write("{:d}&$\state{{{:s}}}$&$\state{{{:s}}}$&{:d}&{:d}&{:d}&{:d}&{:.0f}&{:.0f}&{:.0f}&{:.0f}&{:.0f}\\\\".format(
        l+1,
        state_key[l][0],
        state_key[l][1],
        state_key[l][2],
        state_key[l][3],
        state_key[l][4],
        state_key[l][5],
        data[l,0],
        data[l,1],
        data[l,2],
        data[l,3],
        data[l,4]))
    texfile.write("\n")
texfile.write("\\hline\n\\end{tabular}}\n\n")

# Calculate mean coincidence count rates
C   = np.mean(data,axis=1)
dC  = np.std(data,axis=1)

# Calculate the probabilities
P       = C/np.sum(C)
dP      = 1./np.sum(C)*(dC+np.sum(dC)*np.abs(P[k]))

#print(C)
#print(dC)
#print(P)
#print(dP)

def complexceil(val,dig=3):
    ee  = 10**dig
    return np.ceil(ee*val.real)/ee+1.j*np.ceil(ee*val.imag)/ee

def printcomplexmatrix(name,mat,dmat=np.matrix([False]),dig=3):
    texfile.write("\\newcommand{\\"+name+"}{\n")
    texfile.write("\\begin{pmatrix}\n")
    if dmat.any()!=False:
        pstring = "{{0.real:.{0:d}f}}{{0.imag:+.{0:d}f}}i\\pm {{1.real:.{0:d}f}}{{1.imag:+.{0:d}f}}i".format(dig)
        for k in range(mat.shape[0]):
            texfile.write(pstring.format(mat[k,0],complexceil(dmat[k,0],dig)))
            for l in range(1,mat.shape[1]):
                texfile.write(" &"+pstring.format(mat[k,l],complexceil(dmat[k,l],dig)))
            if k+1 != mat.shape[0]:
                texfile.write("\\\\")
            texfile.write("\n")
    else:
        pstring = "{{0.real:.{0:d}f}}{{0.imag:+.{0:d}f}}i".format(dig)
        for k in range(mat.shape[0]):
            texfile.write(pstring.format(mat[k,0]))
            for l in range(1,mat.shape[1]):
                texfile.write(" &"+pstring.format(mat[k,l]))
            if k+1 != mat.shape[0]:
                texfile.write("\\\\")
            texfile.write("\n")
    texfile.write("\\end{pmatrix}}\n\n")

# Calculate density matrix
rho     = np.matrix(
    [[0.+0.j,0.,0.,0.],
     [0.,0.,0.,0.],   
     [0.,0.,0.,0.],
     [0.,0.,0.,0.]])
drho    = np.matrix(
    [[0.+0.j,0.,0.,0.],
     [0.,0.,0.,0.],   
     [0.,0.,0.,0.],
     [0.,0.,0.,0.]])
for k in range(16):
    state   = f.basis(state_key[k][0],state_key[k][1])
    ketbra  = state*state.H
    rho    += P[k]*ketbra
    drho   += dP[k]*(abs(ketbra.real)+1.j*abs(ketbra.imag))

print("rho")
print(rho)
print("drho")
print(drho)

printcomplexmatrix("rhonormal",rho,dig=3)
printcomplexmatrix("drhonormal",drho,dig=3)

w, v    = eig(rho)
print("Eigenvalues")
print(w)
print("Eigenvectors")
print(v)
print("rho (bell basis)")
print(f.bellmat*rho*f.bellmat)

printcomplexmatrix("rhobell",f.bellmat*rho*f.bellmat,dig=3)
printcomplexmatrix("drhobell",f.bellmat*drho*f.bellmat,dig=3)

#print("S    = {} +- {}".format(f.S(f.a1,f.b1,f.a2,f.b2,rho),f.S(f.a1,f.b1,f.a2,f.b2,drho)))

"""
# 2D plot of alpha and beta with a'=a-45deg b'=b-45deg
abrange = np.linspace(0.,2*np.pi,1000)
X, Y    = np.meshgrid(abrange,abrange)
Z       = np.zeros((1000,1000))
for k in range(1000):
    for l in range(1000):
        Z[k,l]    = np.real(f.S(X[k,l],Y[k,l],X[k,l]-np.pi/4.,Y[k,l]-np.pi/4.,rho))

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
"""

texfile.close()
