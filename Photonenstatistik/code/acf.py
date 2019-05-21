# This file calculates the correlation function
import numpy as np
import matplotlib
import subprocess
from scipy.optimize import curve_fit

#   Matplotlib and pgf settings
#######################################################################################
matplotlib.rcParams["errorbar.capsize"]     = 4.
matplotlib.rcParams.update()

def figsize(hnumber,vnumber,linewidth=441.01733,textheight=660.10394):
        """
        Calculate figure size from number of horizontal and vertical plots.
        """
        fig_width = linewidth / 72.27 / hnumber
        fig_height = textheight / 72.27 / vnumber
        return [fig_width,fig_height]

#matplotlib.use("pgf")

if False:
    matplotlib.rcParams['axes.grid']        = True
    matplotlib.rcParams['grid.linestyle']   = '--'

pgf_with_pdflatex = {
    "pgf.texsystem"     : "pdflatex",  # change this if using xetex or lautex
    "text.usetex"       : True,        # use LaTeX to write all text
    "font.family"       : "serif",
    "font.serif"        : [],          # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif"   : [],
    "font.monospace"    : [],
    "figure.figsize"    : figsize(1.1,2.3),
    "savefig.pad_inches": 0,
    "savefig.bbox"      : 'tight',
    "axes.xmargin"      : .0,
    "axes.ymargin"      : .0,
    "pgf.preamble"      : [
         r"\usepackage[utf8]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         r"\usepackage[detect-all,locale=DE,separate-uncertainty=true, exponent-product=\cdot]{siunitx}",
         r"\usepackage{mathtools}",
         r"\usepackage{physics}",
         r"\usepackage{amsmath}",
         r"\usepackage{upgreek}"
         ]
    }
matplotlib.rcParams.update(pgf_with_pdflatex)
import matplotlib.pyplot as plt

#   Calculation
############################################################################################
Theory  = \
"""
g(t)    = (CN(t) - (1-rho**2))/rho**2
CN(t)   = c(t)/(N1*N2*omega*T)
rho     = S/(S+B)
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
N1      = [11000.,12000.,120000.]   # counts/s
N2      = [20000.,21000.,200000.]   # counts/s
T       = [6128.20,6525.34,460.99]  # s

# preparing fit functions
def g_exp(t,tau,a,h):
    return h - a*np.exp(-np.abs(t)/tau)

# calculations
subprocess.call(["mkdir","-p","../plots"])
for k in (0,1,2):
    print("measurement {}".format(k+1))
    t, c    = np.loadtxt("../Photonenstatistik_Messdaten/messung_{}.dat".format(k+1),unpack=True)
    
    # calculations
    S       = N1[k] + N2[k] - B
    print("    S        = {:n} counts/s".format(S))
    rho     = S/(S+B)
    print("    rho      = {:.3f}".format(rho))
    CN      = c/(N1[k]*N2[k]*omega*T[k])
    #CN     /= CN.mean()
    g       = (CN-(1-rho**2))/rho**2
    #g      /= g.mean()
    
    # fit (separate for + and -)
    gp  = g[t>=0.]
    tp  = t[t>=0.]
    tm  = t[t<=0.]
    gm  = g[t<=0.]
    p0  = [10.,g.min(),g.mean()]
    poptp,  pcovp   = curve_fit(g_exp,tp,gp,p0=p0)
    poptm,  pcovm   = curve_fit(g_exp,tm,gm,p0=p0)
    print("    poptp    = {}".format(poptp))
    print("    pcovp    = {}".format(np.diag(pcovp)))
    print("    poptm    = {}".format(poptm))
    print("    pcovm    = {}".format(np.diag(pcovm)))
    
    # plot
    tsp  = np.linspace(0,500,10000)
    tsm  = np.linspace(-500,0,10000)
    plt.plot(t,g,"x")
    plt.plot(tsp,g_exp(tsp,poptp[0],poptp[1],poptp[2]))
    plt.plot(tsm,g_exp(tsm,poptm[0],poptm[1],poptm[2]))
    plt.ylim((0,2))
    plt.savefig("../plots/g_{}.pgf".format(k+1))
    plt.savefig("../plots/g_{}.pdf".format(k+1))
    plt.close()



   
