#
#                    Sedenion Algebra Foundation
#
#   Functional implementation of 16-dimensional sedenion
#   algebra. Provides multiplication tables and operations
#   on numpy arrays.
#

import numpy as np
from functools import lru_cache
from typing import Tuple, Dict

# Constants
DIMENSION = 16
DIM_EXTERNAL = 8
DIM_INTERNAL = 8

BASIS_NAMES = [
    "e0", "e1", "e2", "e3", "e4", "e5", "e6", "e7",  # External
    "i0", "i1", "i2", "i3", "i4", "i5", "i6", "i7"   # Internal (mapped to e8-e15)
]

@lru_cache(maxsize=1)
def generate_multiplication_table() -> Dict[Tuple[int, int], Tuple[int, int]]:
    """
    Generate 16x16 multiplication table using Cayley-Dickson doubling.
    Returns dict: (i, j) -> (k, sign)
    where e_i * e_j = sign * e_k
    """
    table = {}
    
    # Base case: Re (index 0)
    table[(0, 0)] = (0, 1)
    
    # Recursive doubling: 1->2->4->8->16
    for dim in [1, 2, 4, 8]: 
        new_table = table.copy()
        for i in range(2 * dim):
            for j in range(2 * dim):
                # Cayley-Dickson Formula: (a,b)(c,d) = (ac - d*b, da + bc*)
                
                if i < dim and j < dim:
                    # (e_i, 0)(e_j, 0) = (e_i e_j, 0)
                    new_table[(i, j)] = table[(i, j)]
                    
                elif i < dim and j >= dim:
                    # (e_i, 0)(0, e_J) = (0, e_J e_i)
                    J = j - dim
                    k, sign = table[(J, i)]
                    new_table[(i, j)] = (k + dim, sign)
                    
                elif i >= dim and j < dim:
                    # (0, e_I)(e_j, 0) = (0, e_I e_j*)
                    I = i - dim
                    k, sign = table[(I, j)]
                    conj_sign = 1 if j == 0 else -1
                    new_table[(i, j)] = (k + dim, sign * conj_sign)
                    
                elif i >= dim and j >= dim:
                    # (0, e_I)(0, e_J) = (-e_J* e_I, 0)
                    I = i - dim
                    J = j - dim
                    k, sign = table[(J, I)]
                    conj_sign = 1 if J == 0 else -1
                    new_table[(i, j)] = (k, -1 * sign * conj_sign)
        table = new_table

    return table

def multiply(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Multiply two 16D sedenion vectors.
    a, b: shape (16,) numpy arrays
    """
    if a.shape != (16,) or b.shape != (16,):
        raise ValueError("Inputs must be 16-dimensional vectors")
        
    result = np.zeros(16, dtype=float)
    table = generate_multiplication_table()
    
    # Optimize: Sparse multiplication
    idx_a = np.nonzero(a)[0]
    idx_b = np.nonzero(b)[0]
    
    for i in idx_a:
        val_a = a[i]
        for j in idx_b:
            val_b = b[j]
            k, sign = table[(i, j)]
            result[k] += sign * val_a * val_b
            
    return result

def conjugate(a: np.ndarray) -> np.ndarray:
    """Return conjugate (negate imaginary parts)."""
    res = a.copy()
    res[1:] *= -1
    return res

def norm_sq(a: np.ndarray) -> float:
    return np.sum(a**2)

def inverse(a: np.ndarray) -> np.ndarray:
    n2 = norm_sq(a)
    if n2 == 0:
        raise ValueError("Cannot invert zero sedenion")
    return conjugate(a) / n2

def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return multiply(a, b) - multiply(b, a)

def associator(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    return multiply(multiply(a, b), c) - multiply(a, multiply(b, c))

def split_parts(s: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (external, internal) octonion parts."""
    return s[:8], s[8:]

def create_basis(index: int) -> np.ndarray:
    v = np.zeros(16)
    v[index] = 1.0
    return v
