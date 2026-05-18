from isotope_constants import *

def m_Th230_U238_a_from_concs(Th230_pmolg, U238_ppb):
    '''From ID estimated concs calculate measured activity ratio'''

    Th230_U238_atom = (Th230_pmolg / (U238_ppb / U238_atom)) / 1000

    return Th230_U238_atom * lambda_230 / lambda_238
