import numpy as np
from collections import namedtuple

from isotope_constants import (
    lambda_230, Th230_Th232_ini,
    DETRITAL_CORRECTION_ITERS, DETRITAL_INI_RATIO_RSD,
)
from age_calc_Th230_U238_CS import age_calc_Th230_U238_CS
from m_Th230_U238_a_from_concs import m_Th230_U238_a_from_concs

Measurement = namedtuple('Measurement', ['value', 'error', 'error_type'])


def _sample(m):
    if m.error_type == 'rsd':
        return np.random.normal(m.value, m.value * m.error)
    return np.random.normal(m.value, m.error)


def workers(u238, th230, th232, d234u, T_est=0):
    """Run one MC walk. Returns age in years, or None if age exceeds model limit."""
    U238_ppb   = _sample(u238)
    Th230_pmol = _sample(th230)
    Th232_pmol = _sample(th232)
    m_d234     = _sample(d234u)

    Th230_Th232_ini_r = np.random.normal(
        Th230_Th232_ini, Th230_Th232_ini * DETRITAL_INI_RATIO_RSD
    )

    ini_230_238_a = m_Th230_U238_a_from_concs(Th230_pmol, U238_ppb)
    age_uncorr = age_calc_Th230_U238_CS(T_est, ini_230_238_a, m_d234)

    ini_Yr = age_uncorr
    for _ in range(DETRITAL_CORRECTION_ITERS):
        if ini_Yr is None:
            break
        Th230_232_now = Th230_Th232_ini_r * np.e ** (-lambda_230 * ini_Yr)
        conc_230 = Th230_pmol - (Th230_232_now * Th232_pmol)
        ini_230_238_a = m_Th230_U238_a_from_concs(conc_230, U238_ppb)
        ini_Yr = age_calc_Th230_U238_CS(ini_Yr, ini_230_238_a, m_d234)

    return (age_uncorr, ini_Yr)
