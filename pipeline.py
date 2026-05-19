"""
Preprocessing and MC orchestration for the U-Th age dating pipeline.

Usage:
    import pipeline
    data = pipeline.load_csv('samples.csv')   # handles old and new header formats
    age, err = pipeline.run_mc(data.iloc[0])  # returns (mean_yr, 2sigma_yr)
"""
import multiprocessing as mp
import numpy as np
import pandas as pd

from isotope_constants import lambda_234, lambda_238, U238_atom, Th232_atom
from age_calc_Th230_U238_CS import age_calc_Th230_U238_CS
from m_Th230_U238_a_from_concs import m_Th230_U238_a_from_concs
from workers import workers, Measurement


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

def load_csv(path):
    """Load a U-Th mass-spec CSV (old flat or new hierarchical format).

    Returns a DataFrame with normalised column names and derived measurement
    columns ready for run_mc():
        Sample, U238_ppb, U238_rsd, Th230_pmolg, Th230_pmolg_err,
        Th232_pmol, Th232_err, delta_234U, delta_234U_err
    """
    probe = pd.read_csv(path, nrows=1)
    if '[238U](pmol_g)' in probe.columns:
        return _normalize_old(pd.read_csv(path))
    else:
        return _normalize_new(pd.read_csv(path, header=[0, 1]))


def _normalize_old(df):
    out = pd.DataFrame()
    out['Sample']          = df['Sample']
    out['U238_pmolg']      = df['[238U](pmol_g)']
    out['U238_pmolg_err']  = df['[238U]_er']
    out['Th230_pmolg']     = df['[230Th](pmol_g)']
    out['Th230_pmolg_err'] = df['[230Th]_er']
    out['Th232_Th230']     = df['232/230_Th']
    out['Th232_Th230_err'] = df['232/230_Ther']
    out['U234_U238']       = df['234/238_U']
    out['U234_U238_err']   = df['234/238_Uer']
    if 'age_estimation' in df.columns:
        out['age_estimation'] = df['age_estimation']
    return _add_derived(out)


def _normalize_new(df):
    # Flatten MultiIndex, forward-filling level-0 names over the unnamed interleaved columns.
    # The first column ('Unnamed: 0_level_0', 'Sample') becomes just 'Sample'.
    level0 = pd.Series(df.columns.get_level_values(0))
    level0 = level0.where(~level0.str.startswith('Unnamed:')).ffill()
    level1 = df.columns.get_level_values(1)
    df.columns = [str(b) if pd.isna(a) else f'{a}|{b}' for a, b in zip(level0, level1)]
    out = pd.DataFrame()
    out['Sample']          = df['Sample']
    out['U238_pmolg']      = df['[238U]s (pmol/g)|Average'].astype(float)
    out['U238_pmolg_err']  = df['[238U]s (pmol/g)|uc (k=1)'].astype(float)
    out['Th230_pmolg']     = df['[230Th]s (pmol/g)|Average'].astype(float)
    out['Th230_pmolg_err'] = df['[230Th]s (pmol/g)|uc (k=1)'].astype(float)
    out['Th232_Th230']     = df['232Th/230Th|Average'].astype(float)
    out['Th232_Th230_err'] = df['232Th/230Th|uc (k=1)'].astype(float)
    out['U234_U238']       = df['234U/238U|Average'].astype(float)
    out['U234_U238_err']   = df['234U/238U|uc (k=1)'].astype(float)
    return _add_derived(out)


def _add_derived(df):
    ratio_eq = lambda_238 / lambda_234

    df['U238_ppb']       = df['U238_pmolg'] * U238_atom / 1000
    df['U238_rsd']       = df['U238_pmolg_err'] / df['U238_pmolg']

    df['Th232_pmol']     = df['Th232_Th230'] * df['Th230_pmolg']
    df['Th232_err']      = df['Th232_Th230_err'] * df['Th230_pmolg']

    df['Th232_ppt']          = df['Th232_pmol'] * Th232_atom
    df['Th232_ppt_err']      = df['Th232_err']  * Th232_atom

    R = 1.0 / df['Th232_Th230']
    df['Th230_Th232_1e6']     = R * 1e6
    df['Th230_Th232_1e6_err'] = (df['Th232_Th230_err'] / df['Th232_Th230'] ** 2) * 1e6

    df['delta_234U']     = (df['U234_U238'] / ratio_eq - 1) * 1000
    df['delta_234U_err'] = (df['U234_U238_err'] / ratio_eq) * 1000

    return df


# ---------------------------------------------------------------------------
# MC orchestration
# ---------------------------------------------------------------------------

def run_mc(row, n_walks=10000, T_est=0):
    """Run Monte Carlo age calculation for one preprocessed sample row.

    Parameters
    ----------
    row     : pandas Series with U238_ppb, U238_rsd, Th230_pmolg,
              Th230_pmolg_err, Th232_pmol, Th232_err, delta_234U, delta_234U_err
    n_walks : number of MC iterations (default 10 000)
    T_est   : initial age guess in years; 0 triggers auto-estimate

    Returns
    -------
    (age_uncorr, err_uncorr, age_corr, err_corr) in years,
    or (None, None, None, None) if solver always exceeds MAX_AGE
    """
    u238  = Measurement(row['U238_ppb'],      row['U238_rsd'],        'rsd')
    th230 = Measurement(row['Th230_pmolg'],   row['Th230_pmolg_err'], 'absolute')
    th232 = Measurement(row['Th232_pmol'],    row['Th232_err'],       'absolute')
    d234u = Measurement(row['delta_234U'],    row['delta_234U_err'],  'absolute')

    if T_est == 0:
        act = m_Th230_U238_a_from_concs(row['Th230_pmolg'], row['U238_ppb'])
        T_est = age_calc_Th230_U238_CS(0, act, row['delta_234U'])
        if T_est is None:
            T_est = 500_000

    age_MC_uncorr, age_MC_corr = [], []

    def _cb(result):
        u, c = result
        age_MC_uncorr.append(u)
        age_MC_corr.append(c)

    pool = mp.Pool(mp.cpu_count())
    for _ in range(n_walks):
        pool.apply_async(
            workers,
            args=(u238, th230, th232, d234u, T_est),
            callback=_cb,
        )
    pool.close()
    pool.join()

    valid_u = [a for a in age_MC_uncorr if a is not None]
    valid_c = [a for a in age_MC_corr   if a is not None]
    if not valid_c:
        return None, None, None, None
    return (
        int(np.mean(valid_u)), int(np.std(valid_u) * 2),
        int(np.mean(valid_c)), int(np.std(valid_c) * 2),
    )
