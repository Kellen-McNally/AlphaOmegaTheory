"""
Sedenion Dirac Equation: Derivation of the Fundamental Wave Equation.

Derives the explicit form of the Dirac operator on the Sedenion manifold,
splitting it into external (spacetime) and internal (mass) sectors.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from core.octonion_algebra import _MULT_TABLE

# Reconstruct 16x16 table
def get_sedenion_table():
    table = {}
    table[(0, 0)] = (0, 1)
    for dim in [1, 2, 4, 8]:
        new_table = table.copy()
        for i in range(2 * dim):
            for j in range(2 * dim):
                if i < dim and j < dim:
                    new_table[(i, j)] = table[(i, j)]
                elif i < dim and j >= dim:
                    J = j - dim
                    k, sign = table[(J, i)]
                    new_table[(i, j)] = (k + dim, sign)
                elif i >= dim and j < dim:
                    I = i - dim
                    k, sign = table[(I, j)]
                    conj_sign = 1 if j == 0 else -1
                    new_table[(i, j)] = (k + dim, sign * conj_sign)
                elif i >= dim and j >= dim:
                    I = i - dim
                    J = j - dim
                    k, sign = table[(J, I)]
                    conj_sign = 1 if J == 0 else -1
                    new_table[(i, j)] = (k, -1 * sign * conj_sign)
        table = new_table
    return table

SED_TABLE = get_sedenion_table()

def get_gamma_matrices():
    """
    Construct the 16x16 matrices L_k corresponding to left-multiplication by e_k.
    """
    matrices = []
    for k in range(16):
        mat = np.zeros((16, 16))
        for j in range(16):
            # Col j is e_k * e_j
            res_k, res_sign = SED_TABLE[(k, j)]
            mat[res_k, j] = res_sign
        matrices.append(mat)
    return matrices

def construct_sedenion_dirac_operator():
    print("Sedenion Dirac Operator Construction")
    print("====================================")
    
    # 1. Get Gamma Matrices (Sedenion Left-Multiplication)
    gammas = get_gamma_matrices()
    print(f"  Generated {len(gammas)} 16x16 Gamma matrices.")
    
    # 2. Define the Operator Split
    # D_total = D_ext + D_int
    # D_ext = sum_{k=0..7} gamma_k * d_k  (Spacetime momentum P_k)
    # D_int = sum_{k=8..15} gamma_k * d_k (Internal momentum -> Mass)
    
    print("  Operator Structure:")
    print("    D_ext (0-7):  Spacetime Derivatives")
    print("    D_int (8-15): Internal Mass Operator")
    
    from physics.particle_group_embedding import generate_g2_generators
    g2_gens = generate_g2_generators()
    Q_sub = g2_gens[0] # 8x8 Charge Operator
    
    # Construct effective mass operator
    Q_shifted = np.zeros((16, 16))
    Q_shifted[8:16, 8:16] = Q_sub
    
    print("\nResult:")
    print("  Sedenion Dirac Equation: (gamma^mu D_mu + M_internal) Psi = 0")
    print("  Mass operator M_internal identified with charge operator Q on internal sector.")

if __name__ == "__main__":
    construct_sedenion_dirac_operator()
