# Decay constants — Cheng et al. (2013) Quat. Geochronol. 20, 142–158
# Correction factors applied to λ₂₃₄ (×0.9985) and λ₂₃₀ (×1.0014) per same reference.
lambda_230 = 0.0000091577 * 1.0014
lambda_232 = 0.000000000049475
lambda_234 = 0.0000028263 * 0.9985
lambda_238 = 0.000000000155125

# Atomic masses of U and Th isotopes (IUPAC 2019)
U_atom = 238.028913
U234_atom = 234.040950
U235_atom = 235.043928
U238_atom = 238.050785

Th_atom = 232.03774
Th230_atom = 230.033132
Th232_atom = 232.032

# Secular equilibrium 230/232 initial ratio (ref?)
Th230_Th232_ini = 0.0000044

# Solver configuration — numerical and scientific assumptions
MAX_AGE = 700_000          # age-solver upper limit in years (model breaks down beyond ~500 ka)
CONVERGENCE_DIFF = 1e-5    # activity-ratio convergence threshold for age solver
DETRITAL_CORRECTION_ITERS = 10   # detrital Th-230 correction iterations per MC walk
DETRITAL_INI_RATIO_RSD = 0.5     # 50% RSD on Th230_Th232_ini (ref?)
