"""
Octonion Algebra Foundation

Functional implementation of 8-dimensional octonion algebra.
Provides multiplication tables and operations on numpy arrays.
"""

import numpy as np
from functools import lru_cache
from typing import Tuple, Dict

DIMENSION = 8

@lru_cache(maxsize=1)
def generate_multiplication_table() -> Dict[Tuple[int, int], Tuple[int, int]]:
    """
    Generate 8x8 multiplication table using Cayley-Dickson doubling.
    Returns dict: (i, j) -> (k, sign)
    """
    table = {}
    table[(0, 0)] = (0, 1)
    
    for dim in [1, 2, 4]: # 1->2, 2->4, 4->8
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

def multiply(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Multiply two 8D octonion vectors."""
    if a.shape != (8,) or b.shape != (8,):
        raise ValueError("Inputs must be 8-dimensional vectors")
        
    result = np.zeros(8, dtype=float)
    table = generate_multiplication_table()
    
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
    res = a.copy()
    res[1:] *= -1
    return res

def norm_sq(a: np.ndarray) -> float:
    return np.sum(a**2)

def inverse(a: np.ndarray) -> np.ndarray:
    n2 = norm_sq(a)
    if n2 == 0:
        raise ValueError("Cannot invert zero octonion")
    return conjugate(a) / n2

def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return multiply(a, b) - multiply(b, a)

def associator(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    return multiply(multiply(a, b), c) - multiply(a, multiply(b, c))
