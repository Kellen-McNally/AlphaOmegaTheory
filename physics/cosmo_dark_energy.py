"""
Unified Dark Energy and Expansion Module.

Combines calculations for Dark Energy density, S-curve expansion history,
and Vacuum Stability.
"""

import sys
import os
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from core.constants import (
    C_LIGHT, M_PLANCK_GEV, OMEGA_LAMBDA, DIM_G2, CASIMIR_C3_G2,
    CASIMIR_C2_G2, H0_KM_S_MPC, TRIALITY, RANK_G2,
    OMEGA_LAMBDA_EXP, OMEGA_LAMBDA_ERR,
    T_CURRENT_GYR, T_EQUILIBRIUM_GYR, M_HIGGS_MEASURED, V_HIGGS_MEASURED
)
from utils.logging_config import get_logger, configure_for_cli

logger = get_logger(__name__)

# --- dark_energy.py & dark_energy_rigorous.py ---

def vacuum_energy_from_casimir(cutoff_scale: float = 1e15) -> dict:
    """Calculates vacuum energy density from Cubic Casimir."""
    C3 = CASIMIR_C3_G2
    Lambda = cutoff_scale
    coefficient = C3 / (16 * np.pi**2)
    rho_vac = coefficient * (Lambda ** 4)
    return {'C3': C3, 'rho_vac_numerical': rho_vac}

def critical_density_from_dof(cutoff_scale: float = 1e15) -> dict:
    """Calculates critical density from total gauge degrees of freedom."""
    total_dof = DIM_G2 + RANK_G2
    Lambda = cutoff_scale
    coefficient = total_dof / (16 * np.pi**2)
    rho_crit = coefficient * (Lambda ** 4)
    return {'total_dof': total_dof, 'rho_crit_numerical': rho_crit}

def omega_lambda_rigorous() -> dict:
    """Derives Ω_Λ = ρ_vac/ρ_crit = 11/16."""
    vac = vacuum_energy_from_casimir()
    crit = critical_density_from_dof()
    omega_lambda_predicted = CASIMIR_C3_G2 / (DIM_G2 + RANK_G2)
    error_abs = abs(omega_lambda_predicted - OMEGA_LAMBDA_EXP)
    error_pct = (error_abs / OMEGA_LAMBDA_EXP) * 100
    return {
        'prediction': {'value': omega_lambda_predicted, 'fraction': '11/16'},
        'experiment': {'value': OMEGA_LAMBDA_EXP, 'error': OMEGA_LAMBDA_ERR},
        'comparison': {'error_percent': error_pct}
    }

def dark_energy_summary() -> dict:
    return omega_lambda_rigorous() # Simplified to use rigorous

# --- expansion_scurve.py & expansion_scurve_rigorous.py ---

class ExpansionSCurve:
    """Cosmological expansion equilibrium from G₂ thermodynamics."""
    def __init__(self):
        self.parameters = {
            't0': 8.7, 'k': 0.35, 'r_eq': 0.045, 'r_initial': 0.95
        }
        self.milestones = {
            'cmb_time': 0.38, 'transition': 9.0, 'present': T_CURRENT_GYR,
            'equilibrium_90': 15.0, 'equilibrium_99': 22.0
        }
        self.observations = {'present_ratio': 0.185}

    def matter_ratio_scurve(self, t: np.ndarray) -> np.ndarray:
        r_eq, r_init, k, t0 = self.parameters['r_eq'], self.parameters['r_initial'], self.parameters['k'], self.parameters['t0']
        return r_eq + (r_init - r_eq) / (1 + np.exp(k * (t - t0)))

    def hubble_parameter_scurve(self, t: np.ndarray, H0: float = H0_KM_S_MPC) -> np.ndarray:
        t_decel, t_equil = 9.0, T_EQUILIBRIUM_GYR
        H = np.zeros_like(t)
        mask1 = t < t_decel
        H[mask1] = H0 * (T_CURRENT_GYR / t[mask1])**(2/3)
        mask2 = (t >= t_decel) & (t < t_equil)
        t_norm = (t[mask2] - t_decel) / (t_equil - t_decel)
        H_max = H0 * 1.2
        H[mask2] = H0 + (H_max - H0) * t_norm / (1 + np.exp(5*(t_norm - 0.5)))
        mask3 = t >= t_equil
        H[mask3] = H0 * 0.8
        return H

    def scale_factor_evolution(self, t: np.ndarray, a0: float = 1.0) -> np.ndarray:
        t_present = T_CURRENT_GYR
        a = np.zeros_like(t)
        for i, time in enumerate(t):
            if time <= 9.0:
                a[i] = a0 * (time / t_present)**(2/3)
            elif time <= T_EQUILIBRIUM_GYR:
                H_eff = H0_KM_S_MPC * (1 + 0.1 * (time - 9.0) / 13.0)
                a[i] = a0 * np.exp(H_eff * (time - t_present) / 3.086e19)
            else:
                a[i] = a0 * (time / t_present)
        return a

    def temperature_evolution(self, t: np.ndarray, T0: float = 2.725) -> np.ndarray:
        t_present = T_CURRENT_GYR
        a_ratio = self.scale_factor_evolution(t) / self.scale_factor_evolution(np.array([t_present]))[0]
        T = T0 / a_ratio
        t_equil, T_eq = T_EQUILIBRIUM_GYR, 0.1
        mask_equil = t > t_equil
        T[mask_equil] = T_eq + (T[mask_equil] - T_eq) * np.exp(-(t[mask_equil] - t_equil) / 10.0)
        return T

    def compare_with_lcdm(self, t_max: float = 30.0) -> Dict:
        t = np.linspace(1.0, t_max, 1000)
        H_scurve = self.hubble_parameter_scurve(t)
        a_scurve = self.scale_factor_evolution(t)
        T_scurve = self.temperature_evolution(t)
        H0, Omega_Lambda, t_present = H0_KM_S_MPC, 0.685, T_CURRENT_GYR
        H_infty = H0 * np.sqrt(Omega_Lambda)
        H_lcdm = H_infty * np.ones_like(t)
        a_lcdm = np.exp(H_infty * (t - t_present) / 3.086e19)
        T_lcdm = 2.725 / a_lcdm
        return {
            'time': t,
            'scurve': {'hubble': H_scurve, 'scale_factor': a_scurve, 'temperature': T_scurve, 'expansion_fate': 'equilibrium', 'final_temperature': 'T_eq ≠ 0'},
            'lcdm': {'hubble': H_lcdm, 'scale_factor': a_lcdm, 'temperature': T_lcdm, 'expansion_fate': 'runaway acceleration', 'final_temperature': 'T → 0 (heat death)'}
        }

def expansion_scurve_rigorous() -> dict:
    """Derives S-curve parameters from G₂ geometry."""
    tau, dim, rank, C2, C3 = TRIALITY, DIM_G2, RANK_G2, CASIMIR_C2_G2, CASIMIR_C3_G2
    T_NOW = 13.8; T_0_OBS = 8.7; K_OBS = 0.35; R_EQ_OBS = 0.045
    t_eq = (dim + rank) / (C3 - 1) * T_NOW
    t_0 = (rank / (tau + rank)) * t_eq
    k = tau / t_0
    r_eq = (tau ** 2) / (dim ** 2 + C2)
    return {'parameters': {'t_eq': t_eq, 't_0': t_0, 'k': k, 'r_eq': r_eq}}

# --- vacuum_stability.py ---

FACTOR_42 = TRIALITY * DIM_G2

def higgs_quartic_coupling() -> float:
    return M_HIGGS_MEASURED**2 / (2 * V_HIGGS_MEASURED**2)

def higher_dimensional_operators() -> Dict[str, float]:
    c6_loop = 1.0 / (16.0 * np.pi**2)
    return {"c6_G2": c6_loop * FACTOR_42}

def stabilization_scale(lambda_ew: float = None) -> float:
    if lambda_ew is None: lambda_ew = higgs_quartic_coupling()
    c6_G2 = higher_dimensional_operators()["c6_G2"]
    return np.sqrt(abs(lambda_ew) * M_PLANCK_GEV**2 / c6_G2)

def vacuum_is_stable() -> Tuple[bool, str]:
    phi_stab = stabilization_scale()
    if phi_stab < M_PLANCK_GEV: return True, f"Stabilization at {phi_stab:.2e} GeV"
    return False, "Marginal stability"

if __name__ == "__main__":
    configure_for_cli(verbose=True)
    logger.info("Running Dark Energy & Expansion Analysis")
    omega_lambda_rigorous()
    ExpansionSCurve().compare_with_lcdm()
    expansion_scurve_rigorous()
    vacuum_is_stable()
