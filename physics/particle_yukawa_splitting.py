"""
Yukawa Derivation: Generation Splitting via External Interaction

Investigates the breaking of mass degeneracy by introducing the interaction 
with the External Octonion (Spacetime geometry). Constructs the Dirac Mass 
Matrix based on the interaction between internal and external sectors.
"""

import numpy as np
from core.constants import ALPHA_GUT, TRIALITY
from physics.particle_yukawa_sedenion import SED_TABLE
from physics.particle_yukawa_triality import get_triality_matrix

def calculate_splitting_factors():

def mul_sed(i, j):
    k, s = SED_TABLE[(i, j)]
    return k, s

def compute_interaction_matrix():
    """
    Compute the mixing matrix induced by External Octonion (indices 8-15).
    """
    print("Constructing External Interaction Scores...")
    
    scores = []
    for v_idx in range(16):
        score = 0
        # Sum over i in External (8-15), j in Internal (0-7)
        for i in range(8, 16):
            for j in range(0, 8):
                # [e_i, e_j, e_v]
                # Term 1: (e_i e_j) e_v
                k_ij, s_ij = mul_sed(i, j)
                k_1, s_1 = mul_sed(k_ij, v_idx)
                sign1 = s_ij * s_1
                
                # Term 2: e_i (e_j e_v)
                k_jv, s_jv = mul_sed(j, v_idx)
                k_2, s_2 = mul_sed(i, k_jv)
                sign2 = s_jv * s_2
                
                if k_1 == k_2:
                    if sign1 != sign2:
                        score += 4
                else:
                    score += 2
        scores.append(score)
        
    return scores

def run_splitting_analysis():
    # 1. Get Triality Eigenvectors
    T = get_triality_matrix()
    vals, vecs = np.linalg.eig(T)
    
    # 2. Compute External Interaction Scores
    scores = compute_interaction_matrix()
    
    print("\nExternal Interaction Scores (Basis):")
    print(scores)
    
    # 3. Analyze Scores by Triality Cycle
    cycle_A = [1, 2, 4] 
    cycle_B = [3, 6, 5] 
    
    score_A = [scores[i] for i in cycle_A]
    score_B = [scores[i] for i in cycle_B]
    
    print(f"Cycle A (1,2,4): {score_A}")
    print(f"Cycle B (3,6,5): {score_B}")
    
    mean_A = np.mean(score_A)
    mean_B = np.mean(score_B)
    
    print(f"\nMean Mass A: {mean_A}")
    print(f"Mean Mass B: {mean_B}")
    
    if mean_A != mean_B:
        print(f"Degeneracy Broken. Ratio: {mean_A/mean_B}")
    else:
        print("Degeneracy persists under this interaction.")

if __name__ == "__main__":
    run_splitting_analysis()
