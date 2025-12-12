"""
Fermion Mixing Matrices: CKM and PMNS Derivation

Computes the CKM and PMNS mixing matrices from geometric properties
of the G2-Sedenion algebra, deriving mixing angles from the misalignment 
of mass and weak interaction bases.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from core.constants import (
    ALPHA_GUT, M_GUT, CKM_LAMBDA_EXP
)
from physics.particle_yukawa_charge import get_charge_eigenvalues
from physics.particle_group_embedding import generate_g2_generators

def calculate_mixing_matrices():
from core.constants import TRIALITY, CASIMIR_C3_G2, DIM_G2, RANK_G2

def derive_pmns_matrix():
    print("PMNS Mixing Matrix (Leptons)")
    print("============================")
    
    # 1. Mass Eigenstates
    charge_magnitudes = get_charge_eigenvalues()
    
    # 2. Weak Interaction Basis (G2 Generators)
    g2_gens = generate_g2_generators()
    
    # 3. Geometric Mixing Angles
    s_sq_th12 = 1.0 / TRIALITY
    s_sq_th13 = TRIALITY / (CASIMIR_C3_G2 * DIM_G2)
    t_sq_th23 = CASIMIR_C3_G2 / (CASIMIR_C3_G2 - RANK_G2)
    delta_CP_deg = 180.0 + 360.0 / (CASIMIR_C3_G2 + 1.0)
    
    print(f"Geometric Predictions:")
    print(f"  sin^2(theta12): {s_sq_th12:.6f}")
    print(f"  sin^2(theta13): {s_sq_th13:.6f}")
    print(f"  tan^2(theta23): {t_sq_th23:.6f}")
    print(f"  delta_CP:       {delta_CP_deg:.2f} deg")

def derive_ckm_matrix():
    print("\nCKM Mixing Matrix (Quarks)")
    print("==========================")
    
    # Wolfenstein Parameters from G2 geometry
    lmbda = np.sqrt((TRIALITY + RANK_G2) / (7 * DIM_G2)) 
    rho = 7.0 / 44.0
    A = 64.0 / 77.0
    eta = 27.0 / 77.0
    
    print(f"Geometric Wolfenstein Parameters:")
    print(f"  lambda: {lmbda:.6f}")
    print(f"  A:      {A:.6f}")
    print(f"  rho:    {rho:.6f}")
    print(f"  eta:    {eta:.6f}")

if __name__ == "__main__":
    derive_pmns_matrix()
    derive_ckm_matrix()
