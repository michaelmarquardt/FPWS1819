import numpy as np
import matplotlib
import subprocess
import siunits

'''
Usage of the file:
>>  from export_data import *
>>  texfile(<"path_to_file">)
>>  set_dig(<int>)
>>  set_pre(<"A">)
>>  usepgf(grid=True)               # Includes the header for using pgf in pyplot
>>  import matplotlib.pyplot as plt

Note that Latex cannot interpret numbers in the newcommand identifier which is why you 
shall only use latters.
'''

#   Define file/folder to write in
#######################################################################################
def texfile(filename="../report/values.tex"):
    global TEX
    TEX     = filename
    texfile = open(TEX, "w")
    texfile.close()
    return TEX
    
#   Defaults
#######################################################################################
def set_dig(dig):
    '''Default round to dig digits.'''
    global DIG
    DIG = dig
    return DIG

def set_pre(pre):
    '''Default name prefix (for example to number a measurement).'''
    global PRE
    PRE = pre
    return PRE

# Initialize some globals
global DIG
global PRE
global VAL
global SI
global USEPGF

DIG = 3
PRE = ""
VAL     = []
SI      = []
USEPGF  = False

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

def usepgf(grid=False):
    global USEPGF
    USEPGF = True
    
    matplotlib.use("pgf")
    
    if grid:
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

#   Write out value
#######################################################################################
def val(name, val, dig=None, e=None):
    r'''
    Write out data as normal number.
    Use e to set an exponent:
        val(a,1200,dig=1,e=2) -> \val"+PRE+"a = 12.0e2
    When e = 0:
        val(a,1200,dig=1,e=0) -> \val"+PRE+"a = 1.2e3
    '''
    global TEX
    global DIG
    global PRE
    global VAL
    
    # Test if name already used
    if any(PRE+name == Val for Val in VAL):
        raise NameError(r"Identifier \val"+PRE+name+" already used!")
    VAL.append(PRE+name)
    
    # Write to TEX
    if dig == None:
        dig = DIG
    texfile = open(TEX, "a")
    if e == None:
        prstr   = "\\newcommand{{\\val"+PRE+"{0:s}}}{{{1:."+str(dig)+"f}}}\n"
        texfile.write(prstr.format(name,round(val,dig)))
    elif e == 0:
        E   = int(np.floor(np.log10(val)))
        prstr   = "\\newcommand{{\\val"+PRE+"{0:s}}}{{{1:."+str(dig)+"f}e{2:d}}}\n"
        texfile.write(prstr.format(name,round(val*10**-E,dig),E))
    else:
        prstr   = "\\newcommand{{\\val"+PRE+"{0:s}}}{{{1:."+str(dig)+"f}e{2:d}}}\n"
        texfile.write(prstr.format(name,round(val*10**-e,dig),e))

def si(name, val, err=None, unit="", dig=None, e=None):
    r'''
    Write out data as si number with error if err != None.
    Use e to set an exponent:
        val(a,1200,dig=1,e=2) -> \val"+PRE+"a = 12.0e2
    When e = 0:
        val(a,1200,dig=1,e=0) -> \val"+PRE+"a = 1.2e3
    Shortcuts for unit can be defined in the file siunits.py.
    '''
    global TEX
    global DIG
    global PRE
    global SI
    
    # Test if name already used
    if any(PRE+name == Si for Si in SI):
        raise NameError(r"Identifier \si"+PRE+name+" already used!")
    SI.append(PRE+name)
    
    # Write to TEX
    texfile = open(TEX, "a")
    if dig == None:
        dig = DIG
    if err == None:
        if e == None:
            prstr   = "\\newcommand{{\\si{0:s}}}{{\\SI{{{1:."+str(dig)+"f}}}{{{2:s}}}}}\n"
            texfile.write(prstr.format(name,round(val,dig),siunits.unit(unit)))
        elif e == 0:
            E   = int(np.floor(np.log10(val)))
            prstr   = "\\newcommand{{\\si{0:s}}}{{\\SI{{{1:."+str(dig)+"f}e{2:d}}}{{{3:s}}}}}\n"
            texfile.write(prstr.format(name,round(val*10**-E,dig),E,siunits.unit(unit)))
        else:
            prstr   = "\\newcommand{{\\si{0:s}}}{{\\SI{{{1:."+str(dig)+"f}e{2:d}}}{{{3:s}}}}}\n"
            texfile.write(prstr.format(name,round(val*10**-e,dig),e,siunits.unit(unit)))
    else:
        if e == None:
            prstr   = "\\newcommand{{\\sierr{0:s}}}{{\\SI{{{1:."+str(dig)+"f}\\pm {2:."+str(dig)+"f}}}{{{3:s}}}}}\n"
            texfile.write(prstr.format(name,round(val,dig),round(err,dig),siunits.unit(unit)))
        elif e == 0:
            E   = int(np.floor(np.log10(val)))
            prstr   = "\\newcommand{{\\sierr{0:s}}}{{\\SI{{{1:."+str(dig)+"f}\\pm {2:."+str(dig)+"f}e{3:d}}}{{{4:s}}}}}\n"
            texfile.write(prstr.format(name,round(val*10**-E,dig),round(err*10**-E,dig),E,siunits.unit(unit)))
        else:
            prstr   = "\\newcommand{{\\sierr{0:s}}}{{\\SI{{{1:."+str(dig)+"f}\\pm {2:."+str(dig)+"f}e{3:d}}}{{{4:s}}}}}\n"
            texfile.write(prstr.format(name,round(val*10**-e,dig),round(err*10**-e,dig),e,siunits.unit(unit)))
    texfile.close()
