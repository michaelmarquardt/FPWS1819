from collections import OrderedDict

def unit(ID):
    try:
        return silist[ID]
    except:
        return ID
        
#   SI-units list
#################
silist = OrderedDict()

silist[""]              = r""
silist["m"]             = r"\metre"
silist["m_s"]           = r"\metre\per\second"
silist["V"]             = r"\volt"
silist["A"]             = r"\ampere"
silist["muA"]           = r"\micro\ampere"
silist["bar"]           = r"\bar"
silist["mbar"]          = r"\milli\bar"
