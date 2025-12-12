#!/usr/bin/env python3
#
#   Maximum Atomic Number from G₂ Triality. Predicts Z_max =
#   172 from pure G₂ geometry, matching Fricke (1971)
#   relativistic Dirac-Fock calculations. Key prediction:
#   Period 8 contains 54 = 2τ³ elements
#

from core.constants import TRIALITY, DIM_G2, CASIMIR_C2_G2, CASIMIR_C3_G2
from utils.logging_config import get_logger

logger = get_logger(__name__)


#
#   Calculate maximum atomic number from G₂ invariants.
#   Z_max = C₂(G₂) × (τ × dim(G₂) + 1) = 4 × (3 × 14 + 1) =
#   4 × 43 = 172. Returns: dict (Maximum Z and period
#   structure)
#
def calculate_maximum_z():
    tau = TRIALITY
    dim = DIM_G2
    C2 = CASIMIR_C2_G2
    C3 = CASIMIR_C3_G2

    # Period 8 size from triality
    period_8_size = 2 * (tau ** 3)  # 2 × 27 = 54

    # Alternative: Period 8 = 32 + 2×C₃
    period_8_alt = 32 + 2 * C3  # 32 + 22 = 54

    # Last element of Period 7 (oganesson)
    Z_period_7 = 118

    # Maximum atomic number
    Z_max = C2 * (tau * dim + 1)  # 4 × 43 = 172

    # Verify: Z_max = 118 + 54
    Z_max_alt = Z_period_7 + period_8_size

    return {
        'Z_max': Z_max,
        'Z_max_alt': Z_max_alt,
        'period_8_size': period_8_size,
        'period_8_alt': period_8_alt,
        'Z_period_7': Z_period_7,
        'formula': f'C₂ × (τ × dim + 1) = {C2} × ({tau} × {dim} + 1) = {Z_max}',
        'triality_formula': f'2τ³ = 2 × {tau}³ = {period_8_size}',
    }
