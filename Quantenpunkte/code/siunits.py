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
silist["mum"]           = r"\micro\metre"
silist["nm"]            = r"\nano\metre"
silist["1_m3"]          = r"\per\metre\tothe{3}"
silist["m_s"]           = r"\metre\per\second"
silist["V"]             = r"\volt"
silist["A"]             = r"\ampere"
silist["muA"]           = r"\micro\ampere"
silist["bar"]           = r"\bar"
silist["mbar"]          = r"\milli\bar"
silist["K"]             = r"\kelvin"
silist["dC"]            = r"\degreeCelsius"
silist["eV"]            = r"\electronvolt"
silist["u"]             = r"\atomicmassunit"
silist["mm2"]           = r"\milli\metre\squared"
silist["MHz"]           = r"\mega\hertz"
silist["muA_V"]         = r"\micro\ampere\per\volt"
silist["C"]             = r"\coulomb"
silist["eV_K"]          = r"\electronvolt\per\kelvin"
silist["kg"]            = r"\kilo\gram"
silist["As_Vm"]         = r"\ampere\second\per\volt\per\metre"
silist["ns"]            = r"\nano\second"
silist["mus"]           = r"\micro\second"
silist["ps"]            = r"\pico\second"
