"""
G₂ Casimir Operators and Geometric Predictions.

Computes Casimir invariants for the G₂ Lie algebra and derives physical
predictions including dark energy density and the weak mixing angle directly
from these geometric properties.
"""

import numpy as np
from typing import Dict, Any, Union
from .constants import (
    DIM_G2,
    RANK_G2,
    CASIMIR_C2_G2,
    CASIMIR_C3_G2,
    TRIALITY,
    OMEGA_LAMBDA,
    SIN2_THETA_W
)

def casimir_c2() -> float:
    """Returns the quadratic Casimir invariant C₂(G₂)."""
    return CASIMIR_C2_G2

def casimir_c3() -> float:
    """Returns the cubic Casimir invariant C₃(G₂)."""
    return CASIMIR_C3_G2

def triality_trace() -> float:
    """Returns the triality trace Tr_τ³(G₂), defined as C₃ + 6."""
    return CASIMIR_C3_G2 + 6

def dark_energy_density() -> float:
    """
    Computes dark energy density (Ω_Λ) from G₂ geometric invariants.
    
    Formula: Ω_Λ = C₃ / (dim(G₂) + rank(G₂))
    """
    return CASIMIR_C3_G2 / (DIM_G2 + RANK_G2)

def weak_mixing_angle() -> float:
    """
    Computes the weak mixing angle (sin²θ_W) from G₂ geometry.
    
    Formula: sin²θ_W = τ / (dim(G₂) - 1)
    """
    return TRIALITY / (DIM_G2 - 1)

def alpha_gut() -> float:
    """
    Computes GUT coupling (α_GUT) from G₂ geometry.
    
    Formula: α_GUT = 1 / (τ × dim(G₂)) = 1/42
    """
    return 1.0 / (TRIALITY * DIM_G2)

def casimir_invariants(representation: str = "adjoint") -> Dict[str, Any]:
    """
    Computes Casimir invariants for a specified G₂ representation.

    Args:
        representation: 'adjoint' (14D) or 'fundamental' (7D).

    Returns:
        Dict containing dimensions and Casimir values.
    
    Raises:
        ValueError: If an unknown representation is specified.
    """
    if representation in ["adjoint", "14"]:
        return {
            "representation": "adjoint",
            "dimension": DIM_G2,
            "C2": CASIMIR_C2_G2,
            "C3": CASIMIR_C3_G2,
            "triality_trace": triality_trace(),
            "physical_meaning": "Gauge bosons"
        }
    elif representation in ["fundamental", "7"]:
        return {
            "representation": "fundamental",
            "dimension": 7,
            "C2": 2.0,
            "C3": 0.0,
            "physical_meaning": "Fermion matter fields"
        }
    else:
        raise ValueError(f"Unknown G₂ representation: {representation}")

def g2_predictions() -> Dict[str, Any]:
    """Aggregates all physical predictions derived from G₂ Casimir operators."""
    
    omega_lambda_calc = dark_energy_density()
    sin2_theta_calc = weak_mixing_angle()
    
    return {
        "casimir_values": {
            "C2_G2": casimir_c2(),
            "C3_G2": casimir_c3(),
            "triality_trace": triality_trace()
        },
        "physical_predictions": {
            "dark_energy_density": {
                "value": omega_lambda_calc,
                "formula": "C₃/(dim + rank)",
                "experimental": OMEGA_LAMBDA,
                "agreement": f"{(1 - abs(omega_lambda_calc - OMEGA_LAMBDA)/OMEGA_LAMBDA)*100:.2f}%"
            },
            "weak_mixing_angle": {
                "value": sin2_theta_calc,
                "formula": "τ/(dim - 1)",
                "experimental": SIN2_THETA_W,
                "agreement": f"{(1 - abs(sin2_theta_calc - SIN2_THETA_W)/SIN2_THETA_W)*100:.2f}%"
            }
        },
        "mathematical_origin": {
            "source": "G₂ Lie algebra Casimir operators",
            "uniqueness": "C₃(G₂) exists only for G₂ among simple algebras",
            "constraint": "Determined by group theory"
        }
    }