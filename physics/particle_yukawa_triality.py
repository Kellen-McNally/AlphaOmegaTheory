"""
Yukawa Derivation: Triality Eigenstate Mass Spectrum

Calculates the Associator Mass expectation values for the Triality Eigenspaces
(V0, V1, V2) to determine the fermion generation mass hierarchy.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from physics.particle_yukawa_sedenion import SED_TABLE, get_sedenion_associator_norm

def mul_sed(i, j):
    k, s = SED_TABLE[(i, j)]
    return k, s

def apply_triality(vector):
    """Apply Triality operator T to a 16D vector."""
    new_vec = np.zeros(16, dtype=complex)
    
    # Cycles definition
    cycle_A = [1, 2, 4]
    cycle_B = [3, 6, 5]
    fixed = [0, 7]
    
    # Apply to first octonion (0-7)
    for idx in fixed:
        new_vec[idx] = vector[idx]
        
    for i in range(3):
        src = cycle_A[i]
        dst = cycle_A[(i+1)%3]
        new_vec[dst] = vector[src]
        
        src = cycle_B[i]
        dst = cycle_B[(i+1)%3]
        new_vec[dst] = vector[src]
        
    # Apply to second octonion (8-15) - mirror symmetry
    for idx in fixed:
        new_vec[idx+8] = vector[idx+8]
        
    for i in range(3):
        src = cycle_A[i] + 8
        dst = cycle_A[(i+1)%3] + 8
        new_vec[dst] = vector[src]
        
        src = cycle_B[i] + 8
        dst = cycle_B[(i+1)%3] + 8
        new_vec[dst] = vector[src]
        
    return new_vec

def get_triality_matrix():
    """Construct 16x16 matrix for T."""
    mat = np.zeros((16, 16), dtype=complex)
    basis = np.eye(16)
    for i in range(16):
        mat[:, i] = apply_triality(basis[i])
    return mat

def compute_associator_mass_matrix():
    """
    Construct the Mass Operator M_{mn} = <e_m | Mass | e_n>
    M = Sum_{a,b} A_{a,b}^T A_{a,b} where A_{a,b} is the associator operator.
    """
    print("Constructing Mass Matrix M...")
    M = np.zeros((16, 16))
    
    for i in range(16):
        for j in range(16):
            A_ij = np.zeros((16, 16))
            k_ij, s_ij = mul_sed(i, j)
            
            for k in range(16):
                # Term 1: (e_i e_j) e_k
                k_res1, s_res1 = mul_sed(k_ij, k)
                sign1 = s_ij * s_res1
                
                # Term 2: e_i (e_j e_k)
                k_jk, s_jk = mul_sed(j, k)
                k_res2, s_res2 = mul_sed(i, k_jk)
                sign2 = s_jk * s_res2
                
                if k_res1 == k_res2:
                    val = sign1 - sign2
                    if val != 0:
                        A_ij[k_res1, k] += val
                else:
                    A_ij[k_res1, k] += sign1
                    A_ij[k_res2, k] -= sign2
            
            M += A_ij.T @ A_ij
            
    return M

def run_analysis():
    # 1. Get Triality Matrix and Eigenvectors
    T = get_triality_matrix()
    vals, vecs = np.linalg.eig(T)
    
    w = np.exp(2j * np.pi / 3)
    
    v0_indices = []
    v1_indices = [] 
    v2_indices = [] 
    
    print("Triality Eigenspace Decomposition:")
    for i, val in enumerate(vals):
        if np.isclose(val, 1):
            v0_indices.append(i)
        elif np.isclose(val, w):
            v1_indices.append(i)
        elif np.isclose(val, w.conjugate()):
            v2_indices.append(i)
            
    print(f"  V0 (Singlet): dim={len(v0_indices)}")
    print(f"  V1 (Gen 1):   dim={len(v1_indices)}")
    print(f"  V2 (Gen 2):   dim={len(v2_indices)}")
    
    # 2. Construct Mass Matrix
    M = compute_associator_mass_matrix()
    
    # 3. Project Mass Matrix onto Eigenspaces
    def analyze_subspace(indices, name):
        print(f"\nSubspace Analysis: {name}")
        if not indices:
            print("  Empty subspace.")
            return []
            
        sub_basis = vecs[:, indices]
        M_sub = sub_basis.conj().T @ M @ sub_basis
        m_vals = np.real(np.linalg.eigvals(M_sub))
        m_vals.sort()
        
        print(f"  Mass Eigenvalues: {np.round(m_vals, 2)}")
        return m_vals

    m0 = analyze_subspace(v0_indices, "V0 (Singlet)")
    m1 = analyze_subspace(v1_indices, "V1 (Gen 1)")
    m2 = analyze_subspace(v2_indices, "V2 (Gen 2)")
    
    print("\nFinal Mass Spectrum Summary:")
    unique_masses = set()
    for mlist in [m0, m1, m2]:
        for m in mlist:
            unique_masses.add(round(m, 2))
            
    print(f"  Distinct Values: {sorted(list(unique_masses))}")

if __name__ == "__main__":
    run_analysis()
