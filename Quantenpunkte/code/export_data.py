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
global SILINE
global USEPGF

DIG = 3
PRE = ""
VAL     = []
SI      = []
SILINE  = []
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

def usepgf(size=(1.1,2.3), grid=False):
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
        "figure.figsize"    : figsize(size[0],size[1]),
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

def ceil(val,dig):
    return np.ceil(val*10**dig)*10**-dig

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
    
    # Without error
    if e == None:
        prstr   = "\\newcommand{{\\si"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}}}{{{2:s}}}}}\n"
        texfile.write(prstr.format(name,round(val,dig),siunits.unit(unit)))
    elif e == 0:
        E   = int(np.floor(np.log10(np.abs(val))))
        prstr   = "\\newcommand{{\\si"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}e{2:d}}}{{{3:s}}}}}\n"
        texfile.write(prstr.format(name,round(val*10**-E,dig),E,siunits.unit(unit)))
    else:
        prstr   = "\\newcommand{{\\si"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}e{2:d}}}{{{3:s}}}}}\n"
        texfile.write(prstr.format(name,round(val*10**-e,dig),e,siunits.unit(unit)))
    
    # With error
    if err != None:
        
        # Value with error
        if e == None:
            prstr   = "\\newcommand{{\\sifull"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}\\pm {2:."+str(dig)+"f}}}{{{3:s}}}}}\n"
            texfile.write(prstr.format(name,round(val,dig),ceil(err,dig),siunits.unit(unit)))
        elif e == 0:
            E   = int(np.floor(np.log10(np.abs(val))))
            prstr   = "\\newcommand{{\\sifull"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}\\pm {2:."+str(dig)+"f}e{3:d}}}{{{4:s}}}}}\n"
            texfile.write(prstr.format(name,round(val*10**-E,dig),ceil(err*10**-E,dig),E,siunits.unit(unit)))
        else:
            prstr   = "\\newcommand{{\\sifull"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}\\pm {2:."+str(dig)+"f}e{3:d}}}{{{4:s}}}}}\n"
            texfile.write(prstr.format(name,round(val*10**-e,dig),ceil(err*10**-e,dig),e,siunits.unit(unit)))
        
        # Only error
        if e == None:
            prstr   = "\\newcommand{{\\sierr"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}}}{{{2:s}}}}}\n"
            texfile.write(prstr.format(name,ceil(err,dig),siunits.unit(unit)))
        elif e == 0:
            E   = int(np.floor(np.log10(np.abs(val))))
            prstr   = "\\newcommand{{\\sierr"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}e{2:d}}}{{{3:s}}}}}\n"
            texfile.write(prstr.format(name,ceil(err*10**-E,dig),E,siunits.unit(unit)))
        else:
            prstr   = "\\newcommand{{\\sierr"+PRE+"{0:s}}}{{\\SI{{{1:."+str(dig)+"f}e{2:d}}}{{{3:s}}}}}\n"
            texfile.write(prstr.format(name,ceil(err*10**-e,dig),e,siunits.unit(unit)))
        
        
    texfile.close()

#   sitabular
#######################################################################################
def siline(name, paramlist, varlist, errlist=[None]):
    """
    Uses a list varlist with variables and a errorlist with errors.
    and a paramlist with the digits to round.
    use errlist[i] = 0 for no error.
    """
    global TEX
    global PRE
    global SILINE
    
    # Test if name already used
    if any(PRE+name == Si for Si in SILINE):
        raise NameError(r"Identifier \siline"+PRE+name+" already used!")
    SILINE.append(PRE+name)
    
    if errlist == [None]:
        errlist = [None]*len(paramlist)
    
    # Write to TEX
    texfile = open(TEX, "a")
    texfile.write(r"\newcommand{\siline"+PRE+name+"}{")
    for i in range(len(paramlist)):
        if errlist[i] == None:
            prstr   = "{0:."+str(paramlist[i])+"f}"
            texfile.write(prstr.format(round(varlist[i],paramlist[i])))
        else:
            prstr   = "{0:."+str(paramlist[i])+"f}\\pm {1:."+str(paramlist[i])+"f}"
            texfile.write(prstr.format(round(varlist[i],paramlist[i]),ceil(errlist[i],paramlist[i])))
        if not i==len(paramlist)-1:
            texfile.write(" &")
    texfile.write("}\n")
    texfile.close()

#   Strings for iteration usage
#######################################################################################
alphabet    = ['a','b','c','d','e','f','g','h','i','j','k','l','m',
               'n','o','p','q','r','s','t','u','v','w','x','y','z']
Alphabet    = ['A','B','C','D','E','F','G','H','I','J','K','L','M'
               'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

#   Change a color in '#rrggbb' mode
#######################################################################################
def changecolor(color, rgb=(1,1,1)):
    R   = int(color[1:3],16)
    G   = int(color[3:5],16)
    B   = int(color[5:7],16)
    R   = min(255,R+rgb[0])
    R   = max(R,0)
    G   = min(255,G+rgb[1])
    G   = max(G,0)
    B   = min(255,B+rgb[2])
    B   = max(B,0)
    return "#{:02x}{:02x}{:02x}".format(R,G,B)
