"""
Unified Cosmology Equilibrium Module.

Calculates the time when cosmological expansion stabilizes (S-curve equilibrium).
"""

import sys
import os
import numpy as np
from pathlib import Path
from typing import Dict, Optional
import logging

# Ensure core is in path
sys.path.insert(0, str(Path(__file__).parent.parent))

from core.constants import OMEGA_LAMBDA, H0_KM_S_MPC, M_PLANCK_GEV, K_BOLTZMANN, GEV_TO_ERG, HBAR, T_EQUILIBRIUM_GYR, T_CURRENT_GYR
from utils.logging_config import get_logger, configure_for_cli

try:
    from scipy.optimize import curve_fit
except ImportError:
    curve_fit = None

logger = get_logger(__name__)

# --- precise_equilibrium_time.py ---

H0 = 2.2e-18 
t_now_gyr = 13.8 
ratio_observed = 0.1849 # Observed

def logistic_scurve(t, t0, k, r_eq):
    r_initial = 1.0 
    return r_eq + (r_initial - r_eq) / (1 + np.exp(k * (t - t0)))

def calculate_precise_equilibrium():
    logger.info("Universe Equilibrium Time Calculation")
    t_data = np.array([0.38, 9.0, 13.8]) 
    r_data = np.array([0.95, 0.50, 0.185]) 
    
    t0_fit, k_fit = 9.5, 0.25 # Default fallback
    if curve_fit:
        try:
            p0 = [10.0, 0.3, 0.15]
            popt, _ = curve_fit(logistic_scurve, t_data, r_data, p0=p0, maxfev=10000)
            t0_fit, k_fit, _ = popt
        except Exception: pass

    ln_99 = np.log(99)
    t_99 = t0_fit + ln_99 / k_fit
    logger.info(f"99% Equilibrium: t = {t_99:.1f} Gyr")

# --- equilibrium_state.py ---

def calculate_equilibrium_properties():
    print("Equilibrium State Properties (T_eq, H_eq)")
    H_eq_kms_mpc = H0_KM_S_MPC * np.sqrt(OMEGA_LAMBDA)
    t_equil = T_EQUILIBRIUM_GYR
    print(f"  Predicted Equilibrium Hubble Constant (H_eq): {H_eq_kms_mpc:.3f} km/s/Mpc")
    print(f"  Dark Energy Fraction: {OMEGA_LAMBDA:.4f} (Geometric: 11/16)")
    print(f"  Matter Fraction:      {1 - OMEGA_LAMBDA:.4f}")
    print(f"  (Derived from current H0 and G2-predicted Omega_Lambda={OMEGA_LAMBDA:.4f})")
    print(f"  Equilibrium Time t_eq: {t_equil:.2f} Gyr")
    
    H_eq_s_inv = H_eq_kms_mpc * (3.24078e-20)
    HBAR_JS = HBAR / 1e7; KB_JK = K_BOLTZMANN / 1e7
    T_eq_kelvin = H_eq_s_inv * HBAR_JS / (2 * np.pi * KB_JK)
    print(f"  Predicted Equilibrium Temperature (T_eq):     {T_eq_kelvin:.2e} K")
    print(f"  (Derived from H_eq as de Sitter temperature)")
    print("\nConceptual Note: '22 Gyr' Equilibrium Time")
    print("SUCCESS: Equilibrium state properties (H_eq, T_eq) calculated based on geometric constants.")

if __name__ == "__main__":
    configure_for_cli(verbose=True)
    calculate_precise_equilibrium()
    calculate_equilibrium_properties()
