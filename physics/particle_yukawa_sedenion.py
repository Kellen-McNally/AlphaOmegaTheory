"""
Yukawa Derivation: Sedenion Geometric Mass Eigenvalues

The Octonion spectrum was degenerate (all imaginary units had mass 96).
We now extend the calculation to the full 16D Sedenion algebra.
Hypothesis: The Sedenion non-associativity breaks the degeneracy between
the Internal and External sectors, and potentially within the Internal sector
due to the specific Cayley-Dickson construction.

We compute the "Associator Mass" for all 16 basis elements.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

# We need the full Sedenion multiplication table.
# Since we don't have a pre-computed 16x16 table in core/octonions/algebra.py (only 8x8),
# we will generate it here using the same Cayley-Dickson logic.

def _generate_sedenion_table():
    """Generate 16x16 multiplication table."""
    table = {}
    # Base: Re
    table[(0, 0)] = (0, 1)
    
    for dim in [1, 2, 4, 8]: # 1->2->4->8->16
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

SED_TABLE = _generate_sedenion_table()

def get_sedenion_associator_norm(idx_psi):
    """
    Compute Associator Mass for basis element e_idx in 16D Sedenion.
    Sum over all 16x16 interactions.
    """
    psi = idx_psi
    total_norm = 0
    
    # We only sum over the basis vectors to get the "geometric weight"
    for i in range(16):
        for j in range(16):
            # (e_i e_j) e_psi
            k_ij, s_ij = SED_TABLE[(i, j)]
            k_final1, s_final1 = SED_TABLE[(k_ij, psi)]
            sign1 = s_ij * s_final1
            
            # e_i (e_j e_psi)
            k_jpsi, s_jpsi = SED_TABLE[(j, psi)]
            k_final2, s_final2 = SED_TABLE[(i, k_jpsi)]
            sign2 = s_jpsi * s_final2
            
            # Check difference
            # In Sedenions, result is always a basis vector (up to sign)
            # or zero? Sedenions have zero divisors!
            # But basis vector product is always a basis vector in Cayley-Dickson.
            # Zero divisors happen for linear combinations.
            
            if k_final1 != k_final2:
                total_norm += 2
            else:
                if sign1 != sign2:
                    total_norm += 4
                    
    return total_norm

def run_derivation():
    print("Computing Sedenion Mass Spectrum (16D)...")
    
    masses = []
    for k in range(16):
        m = get_sedenion_associator_norm(k)
        masses.append((k, m))
        
    print("\nSedenion Associator Scores:")
    print(f"{ 'Idx':<4} {'Element':<15} {'Mass':<10} {'Normalized'}")
    print("-" * 45)
    
    # Normalize to the lowest non-zero mass
    non_zeros = [m for k,m in masses if m > 0]
    base_mass = min(non_zeros) if non_zeros else 1
    
    for k, m in masses:
        label = "Real" if k == 0 else ("Octonion" if k < 8 else "Sedenion")
        print(f"{k:<4} {label:<15} {m:<10} {m/base_mass:.2f}")
        
    # Determine if hierarchy emerges
    # We expect e_0 = 0
    # e_1...e_7 (Octonion) might be different from e_8...e_15 (Sedenion extension)
    
    print("\nConclusion:")
    oct_masses = [m for k,m in masses if 0 < k < 8]
    sed_masses = [m for k,m in masses if k >= 8]
    
    avg_oct = np.mean(oct_masses)
    avg_sed = np.mean(sed_masses)
    
    print(f"Internal Octonion Mass (avg): {avg_oct}")
    print(f"External Sedenion Mass (avg): {avg_sed}")
    
    if avg_oct != avg_sed:
        print("SUCCESS: Symmetry breaking between Internal and External sectors found!")
        print(f"Ratio External/Internal: {avg_sed/avg_oct:.4f}")
    else:
        print("FAIL: Still degenerate. Need higher order corrections or specific projection.")

if __name__ == "__main__":
    run_derivation()
