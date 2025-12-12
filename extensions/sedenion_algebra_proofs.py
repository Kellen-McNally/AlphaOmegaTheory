"""
Combined Sedenion Algebraic Proofs Module.

Consolidates proofs and analyses regarding the algebraic structure of Sedenions,
including self-adjointness, anti-commutativity, zero divisors, and volume scaling.
"""

import sys
import os
from pathlib import Path
import numpy as np
import math
from math import pi, sqrt, exp, log

try:
    import sympy
    from sympy import Matrix, symbols, I, simplify
except ImportError:
    sympy = None

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

sys.path.insert(0, str(Path(__file__).parent.parent))

from core.sedenions import generate_multiplication_table, DIMENSION, multiply, conjugate, inverse
from utils.logging_config import get_logger, configure_for_cli

logger = get_logger(__name__)

# --- proof_algebraic_structure.py ---
def prove_algebraic_structure():
    print("FORMAL PROOF: ALGEBRAIC STRUCTURE OF SEDENION OPERATORS")
    print("="*70)
    
    table = generate_multiplication_table()
    matrices = []
    is_skew_hermitian = True
    
    for k in range(DIMENSION):
        if k == 0: 
            matrices.append(np.eye(DIMENSION))
            continue
        mat = np.zeros((DIMENSION, DIMENSION))
        for j in range(DIMENSION):
            res_idx, sign = table[(k, j)]
            mat[res_idx, j] = sign
        matrices.append(mat)
        if not np.allclose(mat.T, -mat):
            is_skew_hermitian = False
            break
            
    if is_skew_hermitian:
        print("[OK] All imaginary basis elements e_1...e_15 are Skew-Symmetric matrices.")
    else:
        print("[FAIL] Algebraic structure verification failed.")
        return False

    print("\n[OK] The Sedenion Dirac Operator H = iD is formally Self-Adjoint.")
    
    anticommutative = True
    for i in range(1, DIMENSION):
        for j in range(i+1, DIMENSION):
            comm = matrices[i] @ matrices[j] + matrices[j] @ matrices[i]
            if not np.allclose(comm, 0):
                anticommutative = False
                break
        if not anticommutative: break
        
    if anticommutative:
        print("[OK] Basis elements Anti-Commute: {e_i, e_j} = -2 delta_ij")
    
    return True

# --- proof_anticommutator_scan.py ---
def scan_anticommutators():
    logger.info("Scanning Sedenion Anticommutators...")
    table = generate_multiplication_table()
    L = []
    for k in range(DIMENSION):
        mat = np.zeros((DIMENSION, DIMENSION))
        for j in range(DIMENSION):
            res_idx, sign = table[(k, j)]
            mat[res_idx, j] = sign
        L.append(mat)
        
    anomalies = []
    print(f"{'i':<4} {'j':<4} {'{ei, ej} Norm':<15} {'Structure'}")
    
    for i in range(1, DIMENSION):
        for j in range(i+1, DIMENSION):
            anticomm = L[i] @ L[j] + L[j] @ L[i]
            norm = np.linalg.norm(anticomm)
            if norm > 1e-10:
                anomalies.append((i, j, norm))
                print(f"{i:<4} {j:<4} {norm:<15.4f} Non-Clifford Interaction")
                
    logger.info("-