import numpy as np

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
        
        #print(C)
        #print(dC)
        #print(E)
        #print(dE)
        
    # Calculate CHSH inequality
    S   = E[0]-E[1]+E[2]+E[3]
    dS  = np.sum(dE)
    
    print("S    = {:.2f} +- {:.2f}".format(S,dS))
    
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
    texfile.write("\\end{tabular}}\n")
    



