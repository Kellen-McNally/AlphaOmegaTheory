"""
Yukawa Derivation: Geometric Mass Eigenvalues

Derives fermion mass hierarchy scaling from the "Associator Density" of the 
Sedenion algebra. Calculates the associativity eigenvalues for the 
Internal Octonion basis states.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from core.octonion_algebra import _MULT_TABLE

BASIS_DIM = 8 # Internal Octonion

def get_associator_norm(idx_psi):
    """
    Compute the 'Associator Mass' for a basis element e_idx.
    Mass ~ Sum_{i,j} || (e_i e_j) e_psi - e_i (e_j e_psi) ||^2
    """
    psi = idx_psi 
    total_norm = 0
    
    for i in range(BASIS_DIM):
        for j in range(BASIS_DIM):
            # Term 1: (e_i e_j) e_psi
            k_ij, s_ij = _MULT_TABLE[(i, j)]
            k_final1, s_final1 = _MULT_TABLE[(k_ij, psi)]
            sign1 = s_ij * s_final1
            
            # Term 2: e_i (e_j e_psi)
            k_jpsi, s_jpsi = _MULT_TABLE[(j, psi)]
            k_final2, s_final2 = _MULT_TABLE[(i, k_jpsi)]
            sign2 = s_jpsi * s_final2
            
            if k_final1 != k_final2:
                total_norm += 2
            else:
                if sign1 != sign2:
                    total_norm += 4
                
    return total_norm

def run_derivation():
    print("Associator Mass Spectrum (Internal Octonion Basis)")
    print("================================================")
    
    masses = []
    for k in range(BASIS_DIM):
        m = get_associator_norm(k)
        masses.append((k, m))
        
    for k, m in masses:
        print(f"  State e_{k}: {m}")
        
    real_mass = masses[0][1]
    imag_masses = [m for k,m in masses if k > 0]
    avg_imag = np.mean(imag_masses)
    
    print("-" * 40)
    print(f"Real (e0) Mass:          {real_mass}")
    print(f"Imaginary (e1-e7) Mean:  {avg_imag}")
    print("-" * 40)
    
    if real_mass == 0 and avg_imag > 0:
        print("Symmetry Breaking: Real axis is massless, Imaginary axes acquire mass.")

if __name__ == "__main__":
    run_derivation()
