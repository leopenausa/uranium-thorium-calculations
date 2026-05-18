from isotope_constants import *
import numpy as np

def calc_Th230_U238_act(d234U, T):
    ''' Function that calculates 230Th/238U activity ratios from delta234U
    measured and age (years)'''

    # Compute first initial d234U_0
    d234U_ini = d234U * np.e**(lambda_234 * T)

    # Compute different terms in age equation (closed system)
    term_1 = np.e**(-lambda_230 * T)
    term_2 = lambda_230 / (lambda_230 - lambda_234)
    term_3 = 1 - np.e**((lambda_234 - lambda_230) * T)

    # Compute final solution
    solution = (term_1 - (d234U/1000 * term_2 * term_3) - 1) * -1

    #return solution
    return solution
