"""
Monte Carlo uncertainty analysis for validated PEM electrolyzer model.
Perturbs calibrated parameters by ±5% and computes resulting variability
in cell voltage and specific energy consumption.

Used to support the generalization/uncertainty argument in Section 4.1.3
of the manuscript and the response to Reviewer 2, Comment 2.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Calibrated parameters (from Section 3.1.2 of manuscript)
# ---------------------------------------------------------------------------
alpha_0  = 0.448938          # transfer coefficient [-]
j0_0     = 5.428292e-2       # exchange current density [A/cm^2]
r_0      = 0.169475          # area-specific resistance [Ohm.cm^2]
kconc_0  = 4e-6              # concentration coefficient [V/(A/cm^2)^2]

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
R   = 8.314                  # gas constant [J/(mol.K)]
F   = 96485                  # Faraday constant [C/mol]
n_e = 2                      # electrons per H2 molecule

# ---------------------------------------------------------------------------
# Operating point at which uncertainty is evaluated
# ---------------------------------------------------------------------------
J_eval = 1.0                 # current density [A/cm^2]
T_C    = 75.0                # temperature [degrees C]
T_K    = T_C + 273.15

# ---------------------------------------------------------------------------
# Calibration RMSE (from Section 4.1)
# ---------------------------------------------------------------------------
RMSE = 0.01899               # V

# ---------------------------------------------------------------------------
# Scaling constant for specific energy
# f1 [kWh/kg] = (2F) / (M_H2 * 3.6e6) * V_cell
#              = 26.59 * V_cell
# ---------------------------------------------------------------------------
M_H2       = 2.016e-3        # molar mass H2 [kg/mol]
k_f1       = 2 * F / (M_H2 * 3.6e6)   # = 26.589 kWh/(kg.V)

# Pareto front width (from tab:pareto_stats, 48 L/min)
E_min_pareto = 40.23         # kWh/kg
E_max_pareto = 50.21         # kWh/kg
pareto_width = E_max_pareto - E_min_pareto

# ---------------------------------------------------------------------------
# Monte Carlo setup
# ---------------------------------------------------------------------------
N_MC = 1000
PERTURB_FRAC = 0.05          # ±5%
rng = np.random.default_rng(42)   # reproducible seed

alpha_s = alpha_0 * (1 + PERTURB_FRAC * (2 * rng.random(N_MC) - 1))
j0_s    = j0_0    * (1 + PERTURB_FRAC * (2 * rng.random(N_MC) - 1))
r_s     = r_0     * (1 + PERTURB_FRAC * (2 * rng.random(N_MC) - 1))
kconc_s = kconc_0 * (1 + PERTURB_FRAC * (2 * rng.random(N_MC) - 1))

# ---------------------------------------------------------------------------
# Evaluate V_cell for each perturbation
# ---------------------------------------------------------------------------
E_ref    = 1.42 - 0.0009 * (T_C - 75.0)
eta_act  = (R * T_K) / (alpha_s * n_e * F) * np.arcsinh(J_eval / (2 * j0_s))
eta_ohm  = r_s * J_eval
eta_conc = kconc_s * J_eval**2
V_cell   = E_ref + eta_act + eta_ohm + eta_conc

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
V_mean = np.mean(V_cell)
V_std  = np.std(V_cell, ddof=1)
V_rel  = 100.0 * V_std / V_mean

f1_mean = k_f1 * V_mean
f1_std  = k_f1 * V_std
f1_rel  = 100.0 * f1_std / f1_mean

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
print("=" * 70)
print("MONTE CARLO UNCERTAINTY ANALYSIS — PEM ELECTROLYZER MODEL")
print("=" * 70)
print(f"Samples                 : {N_MC}")
print(f"Parameter perturbation  : ±{PERTURB_FRAC*100:.0f}%")
print(f"Operating point         : J = {J_eval} A/cm², T = {T_C}°C")
print()
print("--- Cell voltage uncertainty ---")
print(f"V_cell (mean)           : {V_mean:.4f} V")
print(f"V_cell (std)            : {V_std:.4f} V")
print(f"V_cell (relative)       : ±{V_rel:.2f}%")
print(f"V_cell (95% CI)         : ±{1.96*V_std:.4f} V")
print()
print("--- Specific energy uncertainty ---")
print(f"f1 (mean)               : {f1_mean:.2f} kWh/kg")
print(f"f1 (std)                : {f1_std:.4f} kWh/kg")
print(f"f1 (relative)           : ±{f1_rel:.2f}%")
print(f"f1 (95% CI)             : ±{1.96*f1_std:.3f} kWh/kg")
print()
print("--- Comparison to Pareto front width ---")
print(f"Pareto front width      : {pareto_width:.2f} kWh/kg")
print(f"Ratio (f1_std / width)  : {100*f1_std/pareto_width:.2f}%")
print()
print("--- Analytical bound (based on calibration RMSE) ---")
print(f"RMSE-based f1 bound     : ±{k_f1 * RMSE:.3f} kWh/kg")
print(f"Ratio (RMSE bound/width): {100*k_f1*RMSE/pareto_width:.2f}%")
print()
print("=" * 70)
print("INTERPRETATION")
print("=" * 70)
print(f"Parametric uncertainty in f1 is ±{f1_std:.3f} kWh/kg "
      f"({100*f1_std/pareto_width:.1f}% of Pareto width).")
print("This is small compared to the trade-off range, confirming that")
print("optimization conclusions are robust to calibration uncertainty.")
print("=" * 70)
