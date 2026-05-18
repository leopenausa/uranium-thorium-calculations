from isotope_constants import MAX_AGE, CONVERGENCE_DIFF
from calc_Th230_U238_act import calc_Th230_U238_act

def age_calc_Th230_U238_CS(T_ini, m_Th230_U238_a, m_d234U, diff=CONVERGENCE_DIFF):
    ''' Iterate age until computed 230Th/238U activity is within diff of measured.
    Returns None if computed age exceeds MAX_AGE (model limit). '''

    T_cal = T_ini
    Th230_238U_a_estim = calc_Th230_U238_act(m_d234U, T_cal)

    while abs(Th230_238U_a_estim - m_Th230_U238_a) > diff:
        if (Th230_238U_a_estim - m_Th230_U238_a) <= 0:
            T_cal += 1
        else:
            T_cal -= 1

        if T_cal >= MAX_AGE:
            return None
        Th230_238U_a_estim = calc_Th230_U238_act(m_d234U, T_cal)

    return T_cal
