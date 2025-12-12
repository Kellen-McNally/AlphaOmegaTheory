"""
Group Embedding Proof: Solving for G2 Generators

Solves the linear system for derivations D(xy) = D(x)y + xD(y) to find the 14 
generators of the G2 algebra within the Octonion basis.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from core.octonion_algebra import _MULT_TABLE

BASIS = np.eye(8)

def build_derivation_constraints():
    """
    Build matrix A such that A * vec(D) = 0 enforces the derivation property.
    Constraint: D(ei * ej) - D(ei)*ej - ei*D(ej) = 0
    """
    constraints = []
    
    for i in range(8):
        for j in range(8):
            target_k, sign_k = _MULT_TABLE[(i, j)]
            
            for out in range(8):
                row = np.zeros((8, 8)) 
                
                # Term 1: (D(ei * ej))_out
                row[out, target_k] += sign_k
                
                # Term 2: -(D(ei) * ej)_out
                for p in range(8):
                    res_idx, res_sgn = _MULT_TABLE[(p, j)]
                    if res_idx == out:
                        row[p, i] -= res_sgn
                        
                # Term 3: -(ei * D(ej))_out
                for q in range(8):
                    res_idx, res_sgn = _MULT_TABLE[(i, q)]
                    if res_idx == out:
                        row[q, j] -= res_sgn
                        
                constraints.append(row.flatten())
                
    return np.array(constraints)

def run_solver():
    print("G2 Generator Construction (Derivation Solver)")
    print("===========================================")
    
    A = build_derivation_constraints()
    ns = null_space(A)
    dim = ns.shape[1]
    
    print(f"Derivation Space Dimension: {dim}")
    
    if dim == 14:
        generators = []
        for i in range(dim):
            vec = ns[:, i]
            mat = vec.reshape((8, 8))
            generators.append(mat)
            
        print("Structure Verification:")
        g0 = generators[0]
        print(f"  Sample Generator (Row 0/Col 0 Real Constraint):")
        print(f"  Row 0: {g0[0, :]}")
        print(f"  Col 0: {g0[:, 0]}")
        
        su3_gens = []
        target_axis = 7
        for G in generators:
            if np.allclose(G[:, target_axis], 0, atol=1e-5):
                su3_gens.append(G)
        print(f"  Generators fixing e7 (SU(3) candidates): {len(su3_gens)}")
        
        if len(su3_gens) == 8:
            print("  [SUCCESS] SU(3) subgroup identified within G2.")
            print("  Result: G2 (14 dim) contains SU(3) (8 dim).")
            
    else:
        print(f"  [FAIL] Expected 14 generators, found {dim}.")

if __name__ == "__main__":
    run_solver()
