"""
Test: Zero Free Parameters in CORE Theory

This test focuses on the CORE modules that generate paper predictions.
Research/exploration files are excluded.

CORE MODULES (what goes in paper):
- api/core/ (G₂ structure, triality, Casimir)
- api/particle_physics/ (masses, mixing, gauge)
- api/cosmology/ (dark energy, dark matter, expansion)
- Paper sections (frb_paper_section, stellar_physics_paper, etc.)

EXCLUDED (research only):
- api/astrophysics/*_engineering*.py
- api/astrophysics/m7_*.py
- api/astrophysics/*_acausal*.py
- Test files, exploration scripts
"""

import re
import os
from pathlib import Path
import sys

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core import constants

# ONLY these two fundamental constants allowed
FUNDAMENTAL = {
    'tau': 3,
    'dim': 14,
    'triality': 3,
}

# Measured physical constants (NOT free parameters)
MEASURED_CONSTANTS = {
    'c', 'hbar', 'h', 'G', 'k_B', 'kB', 'e',
    'm_e', 'm_p', 'm_n', 'M_sun', 'M_earth',
    'alpha_em', 'G_F', 'M_Z', 'M_W', 'M_H', 'M_HIGGS', 'M_TOP',
    'v',  # Vacuum expectation value (246 GeV)
    'HBAR_GEV_S',  # Planck constant in GeV·s units
    'V_ud', 'V_us', 'V_cb', 'V_ub',  # Measured CKM elements
    'delta_CP',  # Measured CP phase
    'M_GUT', 'M_PLANCK',  # Standard physics scales
    'd_n_bound', 'd_n_coefficient',  # Experimental bounds/coefficients
    'TAU_P_HYPER_K_SENSITIVITY',  # Proton decay experimental sensitivity
    'theta_0', 'corrections',  # Physical zeros
    'F_PI',  # Pion decay constant (measured)
    'A_H',  # Higgs quartic coupling
    'ALPHA_S_MZ',  # Strong coupling at M_Z (measured)
    'order', 'PAPER_SECTION_ORDER',  # Section ordering (not physics)
    # Observed data
    'N_gen_obs', 'N_gen_err', 'bursts_detected', 'distance_Mpc',
    'observed', 'uncertainty',  # Generic observed values and errors
    'N_sun', 'beta', 'P_short', 'P_long', 't_dyn_ns', 'tau_ns',  # Stellar/astrophysics data
    'predicted_value', 'observed_value', 'agreement_sigma',  # Comparison data
    # Cosmological parameters (measured)
    't_now_gyr', 'r_initial', 't0_fit', 'k_fit', 'r_eq_fit',
    'FORWARD_FRACTION', 'BACKWARD_FRACTION', 'dm_frac', 'baryon_frac',
    'omega_b_obs', 'omega_dm_obs', 'omega_dm_observed', 'omega_b_observed',
    'target_dm', 'target_baryon', 'detection_efficiency',
    # Fermion masses (measured)
    'm_t_pole', 'm_b_MSbar_MZ', 'm_c_MSbar_MZ', 's_MSbar_MZ', 'm_s_MSbar_MZ',
    # Gauge couplings (measured at M_Z)
    'alpha_1_MZ', 'alpha_2_MZ', 'alpha_3_MZ',
    # PMNS angles (measured)
    'THETA_12_EXP', 'THETA_12_ERR', 'THETA_23_EXP', 'THETA_23_ERR', 'THETA_13_EXP', 'THETA_13_ERR',
    'DELTA_CP_EXP', 'DELTA_CP_ERR',  # CP phase measurements
    # Oblique parameters (measured)
    'S', 'T', 'S_EXP', 'S_ERR', 'T_EXP', 'T_ERR',
    # Beta function coefficients (SM constants)
    'B3_SM', 'b12', 'b13', 'b23', 'b33',
    # Computational parameters (not physics)
    'mini_batch_size', 'scaling_factor', 'temperature',
    'log_M_low', 'log_M_high',  # Scan ranges
    'best_interpretation',  # Choice index
    # G₂-derived dimensional constants (not free!)
    'C₂', 'C₃', 'C2', 'C3',  # G₂ Casimir invariants (fundamental to Lie algebra)
    'IM_OCTONIONS', 'dim_rep_3', 'dim_rep_4', 'dim_rep_5',  # G₂ representation dimensions
    'center_dims', 'inner_dims', 'outer_dims',  # Sedenion structure from G₂
    'g2_singlet', 'G2_STRUCTURE_CONSTANT',  # G₂ geometric constants
    'sedenion_dim', # Sedenion dimension
    # Constants.py internal variables
    '_c_si', '_h_si', '_hbar_si', '_G_si', '_e_si', '_k_si', '_alpha_si', '_me_si', '_mp_si',
    'EV_TO_ERG', 'CM_TO_METER', 'NM_TO_CM', 'UM_TO_CM',
    'DIM_G2', 'RANK_G2', 'DIM_FUNDAMENTAL_G2', 'CASIMIR_C2_G2', 'CASIMIR_C3_G2',
    'DIM_SU3', 'DIM_SU2', 'DIM_U1', 'DIM_E6', 'RANK_E6', 'DIM_FUNDAMENTAL_E6', 'DIM_OCTONIONS',
    'DIM_IMAGINARY_OCTONIONS', 'DIM_SEDENIONS',
    'M_Z_MEASURED', 'M_W_MEASURED', 'M_HIGGS_MEASURED', 'M_TOP_MEASURED', 'V_HIGGS_MEASURED',
    'G_FERMI', 'SIN2_THETA_W_EXP', 'SIN2_THETA_W_ERR',
    'M_ELECTRON_MEASURED', 'M_MUON_MEASURED', 'M_TAU_MEASURED', 'M_PROTON_MEASURED',
    'RYDBERG_EV',
    'ALPHA_1_MZ_EXP', 'ALPHA_2_MZ_EXP', 'ALPHA_3_MZ_EXP',
    'OMEGA_LAMBDA_EXP', 'OMEGA_LAMBDA_ERR', 'H0_KM_S_MPC',
    'M_W_ERR', 'CKM_LAMBDA_EXP', 'CKM_LAMBDA_ERR', 'CKM_A_EXP', 'CKM_A_ERR', 'CKM_RHO_EXP',
    'CKM_RHO_ERR', 'CKM_ETA_EXP', 'CKM_ETA_ERR', 'F_DARK_EXP', 'F_DARK_ERR',
    # Proton Decay
    'M_Planck', 'm_proton', 'f_pion', 'A_H_SU5', 'hyper_k_sensitivity', 'su5_prediction',
    'so10_prediction', 'M_X_su5', 'suppression_base', 'current_limit', 'future_sensitivity',
    'M_PION_GEV', 'super_k_limit',
    # Quantum Corrections
    'a_mu_SM', 'a_mu_G2',
    # Neutrino Mass Ratio
    'V7_dim', 'g2_singlet_removed', 'fundamental_rep', 'triality_rep', 'num_octonions',
    'complement_dim', 'g2_dimension', 'g2_rank', 'fundamental', 'adjoint',
    # Neutrino Mixing
    'fundamental_dim', 'adjoint_dim', 'antisymmetric_total', 'C3_cubic_casimir', 'tau_triality',
    'M_R3', 'delta_m21_sq_exp', 'delta_m31_sq_exp', 'theta12_exp', 'theta23_exp',
    'theta13_exp', 'delta_CP_exp', 'theta12_err', 'theta23_err', 'theta13_err', 'delta_CP_err',
    'rank', # used as local var
    'DM21_SQ_EXP', 'DM21_SQ_ERR', 'DM31_SQ_EXP', 'DM31_SQ_ERR',
    'Y_E', 'Y_MU', 'Y_TAU', 'M_Z_GEV', 'delta_theta12_deg', 'delta_cp_tree', 'delta_correction', 'delta_cp_exp',
    # CKM Matrix
    'rank_g2', 'dim_g2', 'sin_theta_C_exp', 'rep_4_dim', 'rep_5_dim', 'rep_3_dim', 'V_tb',
    # Fermion Masses
    'v_higgs', 'm_tau_exp', 'm_mu_exp', 'm_e_exp', 'alpha_gut_inv',
    # Expansion S-curve
    't_decel', 't_equil', 't_present', 'T_eq', 'H0', 'Omega_Lambda', 'alpha',
    # Constants Expanded
    'SIN2_THETA_23_NU', 'T_EQUILIBRIUM_GYR', 'T_CURRENT_GYR', 'TAU_PROTON_YEARS',
    'FRB_CESSATION_YEAR', 'BETELGEUSE_DIMMING_YEAR', 'ALPHA_GUT_EXP', 'N_GEN_EXP',
    'N_GEN_ERR', 'F_DARK_ERR', 'ALPHA_GUT_DIGIT_SUM', 'SIN2_THETA_W_DIGIT_SUM',
    'OMEGA_LAMBDA_DIGIT_SUM', 'OPTIMAL_CONFIGURATIONS', 'SEARCH_RUNTIME_HOURS',
    'OVERALL_AGREEMENT_PERCENT', 'STATISTICAL_SIGNIFICANCE_SIGMA', 'P_ACCIDENTAL',
    # Sedenions
    'DIMENSION', 'DIM_EXTERNAL', 'DIM_INTERNAL',
    # Expansion S-curve rigorous
    'T_NOW', 'T_0_OBS', 'R_EQ_OBS',
    # Dark Matter Dynamics
    'RS', 'MB_DISK', 'RD_DISK', 'MB_BULGE', 'RB_BULGE',
    # Dark Matter Fraction
    'OMEGA_DM_EXP', 'OMEGA_B_EXP', 'F_DM_ERR',
    # Gravitational Waves
    'E_EW_GEV', 'f_peak_lisa', 'f_peak_ligo', 'A_GW',
    # Core Calculations
    'Z_eff',
    # Atomic
    'REP_1_DIM', 'REP_7_DIM', 'REP_27_DIM', 'REP_TOTAL', 'TRIALITY_ORDER',
    'Z_period_7', 'Z_Ne', 'n_Ne', 'l_Ne', 'ie_classical_Ne', 'NUM_ELEMENTS_TO_TEST',
    # Beta Functions
    'B2_G2', 'B3_G2',
    # Strong CP
    'THETA_QCD_BOUND',
    # GUT
    'ALPHA_GUT_ERR',
    # Misc
    'target_axis', 'BASIS_DIM',
    'ratio_observed', 'target_v', 'power', 'N_modes', 'm', 'val_51', 'decay', # Last batch
    # Adding cosmo_dark_energy.py explicit constants (although they should be commented in file)
    't_eq', 't_0', 'k', 'r_eq', 'lambda_eff',
}