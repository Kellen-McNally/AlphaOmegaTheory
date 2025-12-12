"""
Core mathematical foundations of the αΩ Framework.

This module contains the essential G₂ theory and sedenion algebra that
forms the mathematical foundation for all physics predictions.
"""

from .constants import *
from .casimir import casimir_c2, casimir_c3, dark_energy_density, weak_mixing_angle
from . import sedenions

__all__ = [
    # Constants
    "TRIALITY", "DIM_G2", "ALPHA_GUT", "OMEGA_LAMBDA", "SIN2_THETA_W",

    # Casimir
    "casimir_c2", "casimir_c3", "dark_energy_density", "weak_mixing_angle",

    # Sedenions
    "sedenions",
]
