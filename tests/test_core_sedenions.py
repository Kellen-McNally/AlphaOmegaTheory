import pytest
import numpy as np
from core.sedenions import (
    generate_multiplication_table, 
    multiply, 
    conjugate, 
    norm_sq, 
    inverse, 
    associator,
    create_basis,
    DIMENSION
)

def test_dimension():
    assert DIMENSION == 16

def test_table_generation():
    table = generate_multiplication_table()
    assert len(table) == 16 * 16
    assert table[(0, 0)] == (0, 1) # 1*1 = 1
    assert table[(1, 1)] == (0, -1) # e1*e1 = -1

def test_multiplication_basis():
    e1 = create_basis(1)
    e2 = create_basis(2)
    e3 = create_basis(3)
    
    # e1 * e2 = e3 (in standard Cayley-Dickson?)
    # Let's check what table gives
    prod = multiply(e1, e2)
    # Expect e3
    assert prod[3] == 1.0
    assert prod[0] == 0.0

    # e2 * e1 = -e3
    prod_rev = multiply(e2, e1)
    assert prod_rev[3] == -1.0

def test_associativity_failure():
    # Sedenions are non-associative
    # (e1 e2) e4 != e1 (e2 e4)
    e1 = create_basis(1)
    e2 = create_basis(2)
    e4 = create_basis(4)
    
    lhs = multiply(multiply(e1, e2), e4)
    rhs = multiply(e1, multiply(e2, e4))
    
    # Associator
    assoc = associator(e1, e2, e4)
    
    # Assert non-zero
    assert not np.allclose(lhs, rhs)
    assert not np.allclose(assoc, np.zeros(16))

def test_inverse():
    e1 = create_basis(1)
    inv_e1 = inverse(e1)
    # e1^-1 = -e1
    assert inv_e1[1] == -1.0
    
    prod = multiply(e1, inv_e1)
    assert np.isclose(prod[0], 1.0)
